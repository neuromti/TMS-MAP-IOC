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
num_rep     = 10000;

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
% Output: A cell array containing for each stimulation intensity a matrix for later design matrix bootstrapping 
BOT        = {};
for si=1:7.
    bot_set    = [];
    idx         = find(subAverage.Design(:,end)==si);
    for rep=1:num_rep, 
        bot_set    = cat(1,bot_set,datasample(idx,length(idx))');
    end
    BOT{si} = bot_set';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORM STATISTICAL ANALYSIS FOR EACH VARIABLE OF INTEREST
% DEFINE VARIABLES OF INTEREST

% LOG
logfilename = [folder.code,'Logfile.log'];
logfileid   = fopen(logfilename,'wt');     
fprintf(logfileid,'Num Reps %i \n ',num_rep);
fprintf(logfileid,'do_test %s \n',func2str(do_test));

clear STAT_RESULTS;
VOI = {'ampDATA','mepDATA','latDATA'};
for voi = 1:length(VOI),

    eval(['DATA = subAverage.',VOI{voi},';']);

    % ESTIMATES BOOTSTRAP 95% CI and PERMUTATION TEST P values 
    [bot_CI,P] = utils.do_botandperm_IOC(subAverage,DATA,PERM,BOT,do_test,num_rep);
    
    % AGGREGATE RESULTS FOR EACH VOI
    STAT_RESULTS(voi).bot_CI        = bot_CI;
    STAT_RESULTS(voi).P             = P;

    
    fprintf(logfileid,'Finished %s \n',VOI{voi});
end

% PERFORMS STATISTICAL ANALYSIS FOR EACH TIMEPOINT OF RAW DATA
voi = size(STAT_RESULTS,2)+1;
for t=1:size(subAverage.rawDATA,2)
    
    DATA = subAverage.rawDATA(:,t);

    % ESTIMATES BOOTSTRAP 95% CI and PERMUTATION TEST P values     
    % bot_CI Output: A matrix containing for each intensity the bootstrapped upper and lower  95% CI for each
    % factor and level (CiLow,CiUp,StimulationIntensity,Level,Factor)
    % P Output: A matrix containing for each stimulation intensity the p-value of
    % a two-sided hypothesis test

    [bot_CI,P] = utils.do_botandperm_IOC(subAverage,DATA,PERM,BOT,do_test,num_rep);
    % AGGREGATE RESULTS FOR EACH VOI
    STAT_RESULTS(voi).bot_CI(:,:,:,:,t)     = bot_CI;
    STAT_RESULTS(voi).P(:,:,t)              = P;
    
    fprintf(logfileid,'Finished rawDATA Timepoint %i \n',t);
end

% ANALYSIS FINISHED FOR EACH VARIABLE OF INTEREST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([folder.stat_results,'ioc.mat'],'STAT_RESULTS')
fprintf(logfileid,'Everything saved');

fclose(logfileid);
%%





















%% TO DO PROPER VISUALIZATION
close all
figure
plot(STAT_RESULTS(2).bot_CI(:,:,1,1)','k')
hold on
plot(STAT_RESULTS(2).bot_CI(:,:,2,1)','r')
%%

close all
for k=1:3,
    h_fig = figure;
    h_axe = gca;
    set(gcf,'Position',[100 100 800 600],'paperpositionmode','auto')
    hold on

    switch int8(k)
        case 1, 
            h_plot = plot(BI)
            h_lgnd = legend(setup.IO.label.BI);
        case 2, 
            h_plot = plot(LM);
            h_lgnd = legend(setup.IO.label.LM);
        case 3, 
            h_plot = plot(M1);
            h_lgnd = legend(setup.IO.label.M1);
        otherwise
            disp('Error, impossible case');
    end
    p_any = sum([P(:,k)<.05,P(:,k)<.01,P(:,k)<.001],2);    
    for si=1:7,
        if any(p_any(si,:))
            plot(si,26,'k*');
        end
    end
    
    set(h_axe,'XLIM',[0 8])
    set(h_axe,'YLIM',[25 30])
    set(h_lgnd,'Location','northwest')
end    