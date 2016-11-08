%% Input-Output-Curve Analysis
% configuration
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
cls; %clc, clear all, close all, matlabrc, addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','headmodel','setup','folder');
%% loading data
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
    for c=1:length(setup.IO.label.all), %all conditions
        filename = ['S',num2str(s),'C',num2str(c),'.mat'];        
        %check whether file exists and load if it does, otherwise throw
        %error
        if sum(ismember({D.name},filename)==1) 
            load([folder.data.ioc,filename]);
        else
            disp(['File ',filename, 'does not exist'])
        end
        
        % check whether complete dataset had to be rejected (due to
        % artifacts or other error
        if all(isnan(ioc.int))
            continue;
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
%% PERMUTATION ANALYSIS
% TODO implement for all variables of interest
perm_DATA = subAverage.ampDATA;
perm_DATA = subAverage.latDATA;
% Design of the permutation vector
num_rep     = 1000;
PERM        = {};
for si=1:7,
    perm_set    = [];
    idx         = find(subAverage.Design(:,end)==si);
    for rep=1:num_rep, 
        perm_set    = cat(1,perm_set,idx(randperm(length(idx)))');
    end
    PERM{si} = perm_set';
end

%
% performing the permutation, based on F-Values from ANOVA
perm_F = [];

for si=1:7,    
    t_sel               = subAverage.Design(:,end)==si;
    t_design            = subAverage.Design(sel,1:4);        
    for rep=1:num_rep,
        t_idx               = PERM{si}(:,rep);
        t_values            = perm_DATA(t_idx);                         
        [p,tab,stats]       = anovan(t_values,t_design,'random',4,'display','off','varnames',{'BI','LM','M1','SUB'});
        perm_F(rep,si,:)    = [tab{2:5,6}];
    end   
end
% calculating the maesured (true) F-Values using ANOVA
true_F = [];
for si=1:7,
    t_sel               = (subAverage.Design(:,end)==si);
    t_design            = subAverage.Design(sel,1:4);      
    t_idx               = (subAverage.Design(:,end)==si);
    t_values            = perm_DATA(t_idx,:);
    [p,tab,stats]       = anovan(t_values,t_design,'random',4,'display','off','varnames',{'BI','LM','M1','SUB'});
    true_F(si,:)        = [tab{2:5,6}];
end
% estimating the p-values
P = [];
for si=1:7,
    matrix_P    = repmat(true_F(si,:),num_rep,1)>=squeeze(perm_F(:,si,:));
    t_p         = mean(matrix_P);
    t_p         = min(t_p,1-t_p).*2;
    P(si,:)     = t_p;
end

%calculation of averages
BI = [];
LM = [];
M1 = [];
for si=1:7
    for k=1:3,
    a           = nanmean(perm_DATA(subAverage.Design(:,k)==1 & subAverage.Design(:,5)==si));
    b           = nanmean(perm_DATA(subAverage.Design(:,k)==0 & subAverage.Design(:,5)==si));
    switch int8(k)
        case 1, 
            BI(si,:)    = [a,b];    
        case 2, 
            LM(si,:)    = [a,b];
        case 3, 
            M1(si,:)    = [a,b];
        otherwise
        disp('Error, impossible case');
    end
    end    
end


%% TO DO 

%visualization
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


%% BOOTSTRAP ANALYSIS (TODO)

% Design of the permutation vector
num_rep     = 1000;
PERM        = {};
for si=1:7.
    perm_set    = [];
    idx         = find(subAverage.Design(:,end)==si);
    for rep=1:num_rep, 
        perm_set    = cat(1,perm_set,idx(randperm(length(idx)))');
    end
    PERM{si} = perm_set';
end

%
% performing the permutation, based on F-Values from ANOVA
perm_F = [];

for si=1:7,    
    t_sel               = subAverage.Design(:,end)==si;
    t_design            = subAverage.Design(sel,1:4);        
    for rep=1:num_rep,
        t_idx               = PERM{si}(:,rep);
        t_values            = subAverage.ampDATA(t_idx);                         
        [p,tab,stats]       = anovan(t_values,t_design,'random',4,'display','off','varnames',{'BI','LM','M1','SUB'});
        perm_F(rep,si,:)    = [tab{2:5,6}];
    end   
end
% calculating the maesured (true) F-Values using ANOVA
true_F = [];
for si=1:7,
    t_sel               = (subAverage.Design(:,end)==si);
    t_design            = subAverage.Design(sel,1:4);      
    t_idx               = (subAverage.Design(:,end)==si);
    t_values            = subAverage.ampDATA(t_idx,:);
    [p,tab,stats]       = anovan(t_values,t_design,'random',4,'display','off','varnames',{'BI','LM','M1','SUB'});
    true_F(si,:)        = [tab{2:5,6}];
end
% estimating the p-values 
P = [];
for si=1:7,
    matrix_P    = repmat(true_F(si,:),num_rep,1)>=squeeze(perm_F(:,si,:));
    t_p         = mean(matrix_P);
    t_p         = min(t_p,1-t_p).*2;
    P(si,:)     = t_p;
end







