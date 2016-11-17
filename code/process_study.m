%% Input-Output-Curve Analysis
% configuration
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
cls; %clc, clear all, close all, matlabrc, fclose('all'), addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','headmodel','setup','folder');

%% loading data
logfilename = [folder.code,'Logfile.log'];
logfileid   = fopen(logfilename,'wt');
D           = dir(folder.data.ioc);
D([1,2])    = [];
%define data matrix
latDATA     = []; %latenz data 
ampDATA     = []; %amplitude data
rawDATA     = []; %raw MEP waveform data
subAverage  = struct('ampDATA',[],'mepDATA',[],'latDATA',[],'rawDATA',[],'filtDATA',[],'msoDATA',[],'Design',[],'Condition',[]); %initialize for runwise averaged data
% load data from files and concatenate
for ss=1:length(setup.SUB.id), %every subject
    s = setup.SUB.id(ss);
    fprintf(logfileid,'\n Subject: %i \n',s);    
    for c=1:length(setup.IO.label.all), %all conditions
        filename = ['S',num2str(s),'C',num2str(c),'.mat'];        
        fprintf(logfileid,'Loading %s ',filename);
        %check whether file exists and load if it does, otherwise throw
        %error
        if sum(ismember({D.name},filename)==1) 
            load([folder.data.ioc,filename]);
        else
            disp(['File ',filename, 'does not exist'])
            fprintf(logfileid,'Error: File %s does not exist \n',filename);
        end
        
        % check whether complete dataset had to be rejected (due to
        % artifacts or other error
        if all(isnan(ioc.int))            
            fprintf(logfileid,'Skipped file %s\n',filename);
            continue;
        else
            fprintf(logfileid,'Processing %s\n',filename);
        end                        
        
        % LATENCY
        % concatenate data vectors after cleaning for artifacts (indicated
        % by negative values) and lack of MEP (indicated by zero latency)       
        cleaned             = ioc.lat;
        artifacted          = (ioc.amp<0) | (ioc.lat <0 | ioc.lat == 0);
        cleaned(artifacted) = NaN;       
        t_lat               = nanmean(reshape(cleaned,10,7),1)';
        subAverage.latDATA  = cat(1,subAverage.latDATA,t_lat);
        
        % MEP 
        % concatenate data vectors after cleaning for artifacts (indicated
        % by negative values)        
        cleaned             = ioc.amp;
        artifacted          = (ioc.amp<0) | (ioc.lat <0);
        cleaned(artifacted) = NaN;
        t_amp               = nanmean(reshape(cleaned,10,7),1)';
        t_mep               = nanmean(reshape(cleaned,10,7)>50,1)';
        subAverage.ampDATA  = cat(1,subAverage.ampDATA,t_amp);        
        subAverage.mepDATA  = cat(1,subAverage.mepDATA,t_mep); 
        
        % TIMECOURSE 
        % concatenate data vectors after cleaning for artifacts (indicated
        % by negative values)    
        cleaned             = reshape(ioc.raw,7*10,[]);        
        baselined           = utils.baseline(cleaned,1:100,1);
        filtered            = utils.filt_freqz(baselined,0,250,5000,2,2);
        artifacted          = (ioc.amp<0) | (ioc.lat <0);
        baselined(artifacted,:) = NaN;
        filtered(artifacted,:)  = NaN;        
        t_raw               = reshape(squeeze(nanmean(reshape(baselined,10,7,601),1)),7,[]);
        t_filt              = reshape(squeeze(nanmean(reshape(filtered,10,7,601),1)),7,[]);
        subAverage.rawDATA  = cat(1,subAverage.rawDATA,t_raw);
        subAverage.filtDATA  = cat(1,subAverage.filtDATA,t_filt);
        
        %concatenate factor vector for intensity in MSO and RMT
        mso                 = sort(unique(ioc.int));
        subAverage.msoDATA  = cat(1,subAverage.msoDATA,mso);
        
        %create Design Matrix: BI LM M1 SUB SI       
        subAverage.Design   = cat(1,subAverage.Design,cat(2,repmat(cat(2,setup.IO.BI(c)',setup.IO.LM(c)',setup.IO.M1(c)',s),7,1),[1:7]'));
    end
end
%create Design Matrix for each Condition
subAverage.Condition = cat(2,bi2de(subAverage.Design(:,1:3))+1,subAverage.Design(:,4:end)); 
fclose(logfileid);
%% PERMUTATION ANALYSIS

% GLOBAL SETTING
num_rep     = 1000;

% DEFINITION OF ANALYSIS FUNCTION
% output values are f-values from anova
% subject as random factor
% no interactions
% independently for each SI due to assumption of heteroscedasticity
do_test = @(x,y)anovan(x,y,'random',4,'display','off','varnames',{'BI','LM','M1','SUB'});

% CONSTRUCTION OF THE PERMUTATION ARRAY
% Output: A cell array containing for each stimulation intensity a matrix for later design matrix permutation 

PERM        = {};
for si=1:7,
    perm_set    = [];
    idx         = find(subAverage.Design(:,end)==si);
    for rep=1:num_rep, 
        perm_set    = cat(1,perm_set,idx(randperm(length(idx)))');
    end
    PERM{si} = perm_set';
end

% CONSTRUCTION OF THE BOOTSTRAPPING VECTOR
% Output: A cell array containing for each stimulation intensity a factor-balanced bootstrap matrix
BOT         = {};
orderByte   = bin2dec(num2str(subAverage.Design(1:7:54,1:3)));

for si=1:7.    
    idx         = find(subAverage.Design(:,end)==si);        
    bot_set     = NaN(num_rep,size(idx,1));
    sel_design  = bin2dec(num2str(subAverage.Design(idx,1:3)));
    for i_cond = 1:8, %balancing for each factor, bootstrapping only subjects, not factor levels
        cond_idx = (ismember(sel_design,orderByte(i_cond)));
        cond_sel = idx(cond_idx);
        for rep=1:num_rep, 
            bot_set(rep,cond_idx) = datasample(cond_sel,length(cond_sel));
        end        
    end
    BOT{si} = bot_set';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORM STATISTICAL ANALYSIS FOR EACH VARIABLE OF INTEREST
% DEFINE VARIABLES OF INTEREST

% LOG
logfilename = [folder.code,'Logfile.log'];
if exist(logfilename,'file'),delete(logfilename); end
utils.do_log(logfilename,@()fprintf(logfileid,'Num Reps %i \n ',num_rep));
utils.do_log(logfilename,@()fprintf(logfileid,'Statistical Value is based on %s \n',func2str(do_test)));

clear STAT_RESULTS;
VOI = {'ampDATA','mepDATA','latDATA'};
for voi = 1:length(VOI),

    eval(['DATA = subAverage.',VOI{voi},';']);

    % ESTIMATES BOOTSTRAP 95% CI and PERMUTATION TEST P values 
    tic
    [bot_CI,bot_D,P,true_M] = utils.do_botandperm_IOC(subAverage,DATA,PERM,BOT,do_test,num_rep,setup);
    toc 
    
    % AGGREGATE RESULTS FOR EACH VOI
    STAT_RESULTS(voi).bot_CI        = bot_CI;
    STAT_RESULTS(voi).bot_D         = bot_D;
    STAT_RESULTS(voi).P             = P;
    STAT_RESULTS(voi).M             = true_M;
    
    
    utils.do_log(logfilename,@()fprintf(logfileid,'On %s finished %s \n',datetime('now'),VOI{voi}));
end
save([folder.results.stats,'ioc.mat'],'STAT_RESULTS')
utils.do_log(logfilename,@()fprintf(logfileid,'Variables saved \n'));

% PERFORMS STATISTICAL ANALYSIS FOR EACH TIMEPOINT OF RAW DATA

rawDATA = subAverage.rawDATA;
delete(gcp('nocreate'))
parpool(4)
parfor t=1:size(rawDATA,2)
    
    DATA = rawDATA(:,t);

    % ESTIMATES BOOTSTRAP 95% CI and PERMUTATION TEST P values     
    % bot_CI Output: A matrix containing for each intensity the bootstrapped upper and lower  95% CI for each
    % factor and level (CiLow,CiUp,StimulationIntensity,Level,Factor)
    % P Output: A matrix containing for each stimulation intensity the p-value of
    % a two-sided hypothesis test
    tic
    [bot_CI,bot_D,P,true_M] = utils.do_botandperm_IOC(subAverage,DATA,PERM,BOT,do_test,num_rep,setup);
    clk_time = toc
    % AGGREGATE RESULTS FOR EACH VOI
    STAT_TIMECOURSE(t).bot_CI   = bot_CI;
    STAT_TIMECOURSE(t).bot_D    = bot_D;
    STAT_TIMECOURSE(t).P        = P;  
    STAT_TIMECOURSE(t).M        = true_M;
    utils.do_log(logfilename,@()fprintf(logfileid,'On %s finished rawDATA Timepoint %i after running for %i s \n',datetime('now'),t,int64(clk_time)))
end
delete(gcp('nocreate'))
save([folder.results.stats,'ioc_tc.mat'],'STAT_TIMECOURSE')
utils.do_log(logfilename,@()fprintf(logfileid,'Workers shut down & Timecourse saved \n'));
% ANALYSIS FINISHED FOR EACH VARIABLE OF INTEREST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
utils.do_log(logfilename,@()fprintf(logfileid,'Everything saved'));
%%
addpath('C:\Users\Robert Bauer\Documents\Matlab\other_toolboxes\gramm\')

figpos = [100 100 1200 400];
close all

x   = reshape(repmat(setup.IO.SI,1,2,3),1,[])';
lab = cat(2,reshape(repmat(setup.IO.label.BI,7,1),1,[]),reshape(repmat(setup.IO.label.LM,7,1),1,[]),reshape(repmat(setup.IO.label.M1,7,1),1,[]))';
org = reshape(repmat(setup.IO.label.ORG,14,1),1,[]);
ulab = unique(lab)
ulab = ulab([4,6,1,2,3,5]);
% PLOT LATENCY
y   = reshape(STAT_RESULTS(3).M,1,[])';
ylo = reshape(squeeze(STAT_RESULTS(3).bot_CI(1,:,:,:)),1,[])';
yup = reshape(squeeze(STAT_RESULTS(3).bot_CI(2,:,:,:)),1,[])';

clear g
g   = gramm('x',x,'y',y,'ymin',ylo,'ymax',yup,'color',lab); 
g.geom_interval('geom','area');
g.facet_grid([],org);
g.set_names('column','','x','Stimulation Intensity (RMT)','y','Latency (ms)','color','');
g.set_title('');
g.set_order_options('color',ulab);
g.axe_property('YLim',[25 28],'YTICK',[20:35]);

figure('Position',figpos);
g.draw();
print(gcf,[folder.results.figs,setup.IO.label.VOI{3},'.tif'],'-dtiff')

% PLOT MEP
y   = reshape(STAT_RESULTS(2).M,1,[])';
ylo = reshape(squeeze(STAT_RESULTS(2).bot_CI(1,:,:,:)),1,[])';
yup = reshape(squeeze(STAT_RESULTS(2).bot_CI(2,:,:,:)),1,[])';

clear g
g   = gramm('x',x,'y',y,'ymin',ylo,'ymax',yup,'color',lab); 
g.geom_interval('geom','area');
g.facet_grid([],org);
g.set_names('column','','x','Stimulation Intensity (RMT)','y','P (MEP) in %','color','');
g.set_title('');
g.set_order_options('color',ulab);
g.axe_property('YLim',[0 1],'YTICK',[0:0.1:1]);

figure('Position',figpos);
g.draw();
print(gcf,[folder.results.figs,setup.IO.label.VOI{2},'.tif'],'-dtiff')

% PLOT AMPLITUDE

y   = reshape(STAT_RESULTS(1).M,1,[])';
ylo = reshape(squeeze(STAT_RESULTS(1).bot_CI(1,:,:,:)),1,[])';
yup = reshape(squeeze(STAT_RESULTS(1).bot_CI(2,:,:,:)),1,[])';

clear g
g   = gramm('x',x,'y',y,'ymin',ylo,'ymax',yup,'color',lab); 
g.geom_interval('geom','area');
g.facet_grid([],org);
g.set_names('column','','x','Stimulation Intensity (RMT)','y','Amplitude (Vpp)','color','');
g.set_title('');
g.set_order_options('color',ulab);
g.axe_property('YLim',[0 700]);

figure('Position',figpos);
g.draw();
print(gcf,[folder.results.figs,setup.IO.label.VOI{1},'.tif'],'-dtiff')
%%

figpos = [100 100 1200 400];

x   = reshape(repmat(setup.IO.SI,1,3),1,[])';
lab = reshape(repmat(setup.IO.label.ORG,7,1),1,[])';
org = lab;
% PLOT LATENCY
y   = reshape(squeeze(diff(STAT_RESULTS(3).M,[],2)),1,[])';
ylo = reshape(squeeze(STAT_RESULTS(3).bot_D(1,:,:,:)),1,[])';
yup = reshape(squeeze(STAT_RESULTS(3).bot_D(2,:,:,:)),1,[])';

clear g
g   = gramm('x',x,'y',y,'ymin',ylo,'ymax',yup); 
g.geom_interval('geom','area');
g.facet_grid([],org);
g.set_title('');
g.set_names('column','','x','Stimulation Intensity (RMT)','y','Latency (ms)','color','Parameter');
g.axe_property('YLim',[-2 2],'YTICK',[-2:2]);

figure('Position',figpos);
g.draw();
print(gcf,[folder.results.figs,setup.IO.label.VOI{3},'_delta.tif'],'-dtiff')

% PLOT MEP
y   = reshape(squeeze(diff(STAT_RESULTS(2).M,[],2)),1,[])';
ylo = reshape(squeeze(STAT_RESULTS(2).bot_D(1,:,:,:)),1,[])';
yup = reshape(squeeze(STAT_RESULTS(2).bot_D(2,:,:,:)),1,[])';

clear g
g   = gramm('x',x,'y',y,'ymin',ylo,'ymax',yup); 
g.geom_interval('geom','area');
g.facet_grid([],org);
g.set_title('');
g.set_names('column','','x','Stimulation Intensity (RMT)','y','P (MEP) in %','color','');
g.axe_property('YLim',[-.5 .5],'YTICK',[-.5:.1:.5]);

figure('Position',figpos);
g.draw();
print(gcf,[folder.results.figs,setup.IO.label.VOI{2},'_delta.tif'],'-dtiff')

% PLOT AMPLITUDE

y   = reshape(squeeze(diff(STAT_RESULTS(1).M,[],2)),1,[])';
ylo = reshape(squeeze(STAT_RESULTS(1).bot_D(1,:,:,:)),1,[])';
yup = reshape(squeeze(STAT_RESULTS(1).bot_D(2,:,:,:)),1,[])';

clear g
g   = gramm('x',x,'y',y,'ymin',ylo,'ymax',yup); 
g.geom_interval('geom','area');
g.facet_grid([],org);
g.set_title('');
g.set_names('column','','x','Stimulation Intensity (RMT)','y','Amplitude (Vpp)','color','Parameter');
g.axe_property('YLim',[-300 300]);

figure('Position',figpos);
g.draw();
print(gcf,[folder.results.figs,setup.IO.label.VOI{1},'_delta.tif'],'-dtiff')
%%
M       = cat(4,STAT_TIMECOURSE.M);
CI      = cat(5,STAT_TIMECOURSE.bot_CI);
figpos  = [100 100 1200 400];
close all

raw_lab = cat(1,setup.IO.label.BI,setup.IO.label.LM,setup.IO.label.M1)';
x       = [];
y       = [];
ylo     = [];
yup     = [];
lab     = {};
org     = {};
clear g

for fctr = 1:size(M,3)
    x       = [];
    y       = [];
    ylo     = [];
    yup     = [];
    lab     = {};
    org     = {};
    for si=1:size(M,1)
        for lvl = 1:size(M,2)       
            for t=1:size(M,4)
                x       = cat(1,x,t);
                lab     = cat(1,lab,raw_lab(lvl,fctr));
                org     = cat(1,org,num2str(setup.IO.SI(si)));
                y       = cat(1,y,M(si,lvl,fctr,t));
                ylo     = cat(1,ylo,CI(1,si,lvl,fctr,t));
                yup     = cat(1,yup,CI(2,si,lvl,fctr,t));
            end
        end
    end
    % PLOT LATENCY
    
    clear g
    g       = gramm('x',x,'y',y,'ymin',ylo,'ymax',yup,'color',lab); 
    g.geom_interval('geom','area');
    g.facet_grid([],org);
    g.set_names('column','','x','Time in ms','y','Amplitude in µV','color','');
    g.set_title('Stimulation Intensity (RMT%)');
    g.set_order_options('color',raw_lab(:,fctr)');
    %g.axe_property('YLim',[25 28],'YTICK',[20:35]);

    figure('Position',figpos);
    g.draw();
    print(gcf,[folder.results.figs,setup.IO.label.ORG{fctr},'_TimeCourse.tif'],'-dtiff')    
end

