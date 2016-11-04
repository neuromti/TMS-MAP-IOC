%% Input-Output-Curve Analysis
% configuration
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
cls; %clc, clear all, close all, matlabrc, addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','headmodel','setup','folder');
%% loading data
D           = dir(folder.data.ioc);
D([1,2])    = [];
%define mixing matrix
SUB         = [];
KOND        = []; 
latDATA     = []; %latenz data 
ampDATA     = []; %amplitude data
rawDATA     = []; %raw MEP waveform data
msoINT      = []; %raw stimulation intensity in maximal stimulator output percent (MSO)
rmtINT      = []; %normalized stimulation intensity in resting motor threshold precent (RMT)
RMT         = NaN(length(setup.SUB.id),length(setup.IO.label.all)); %RMT in MSO
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
        
        
        %concatenate design vector for repeatedly measured subject 
        SUB         = cat(1,SUB,s*ones(size(ioc.int))); 
        
        %concatenate design vector for repeatedly measured condition
        KOND        = cat(1,KOND,c*ones(size(ioc.int))); 
        
        %concatenate factor vector for intensity in MSO and RMT
        mso         = sort(unique(ioc.int));
        msoINT      = cat(1,msoINT,ioc.int);
                
        [~, rmt]    = ismember(ioc.int,mso);
        rmtINT      = cat(1,rmtINT,setup.IO.SI(rmt)');

        %concatenate data vectors after cleaning for artifacts (indicated
        %by negative values)
        artifacted          = (ioc.amp<0) | (ioc.lat <0);
        cleaned             = ioc.amp;
        cleaned(artifacted) = NaN;
        ampDATA             = cat(1,ampDATA,cleaned);
        
        cleaned             = ioc.lat;
        cleaned(artifacted) = NaN;
        latDATA             = cat(1,latDATA,cleaned);                
        
        cleaned             = reshape(ioc.raw,7*10,[]);
        cleaned(artifacted,:) = NaN;
        rawDATA             = cat(1,rawDATA,cleaned);
        %estimate RMT in MSO 
        RMT(s,c)            = mso(2); %raw value
        %RMT(s,c)           = mean(diag(mso*(1./[.9,1,1.1,1.2,1.3,1.4,1.5]))); %inverse weighted average based on all seven intensities
                
    end
end
DESIGN = cat(2,setup.IO.BI(KOND)',setup.IO.LM(KOND)',setup.IO.M1(KOND)');
%% calculation of descriptive statistics for MEP parameters (Peak-to-Peak amplitude, Latency) and MEP waveform (raw and shifted to origin at latency)
M   = NaN(7,8,13);
E   = [];
BI  = [];
LM  = [];
M1  = [];
RW  = [];
shiftMEP = [];
for ssii=1:length(setup.IO.SI)
    
    si              = setup.IO.SI(ssii);    
    t_sub           = SUB(ismember(rmtINT,si));
    u_sub           = unique(t_sub);    
    
    t_data          = ampDATA(ismember(rmtINT,si));      
    t_raw           = rawDATA(ismember(rmtINT,si),:);      
    t_kond          = KOND(ismember(rmtINT,si));
    
    for ss=1:length(u_sub), %for all subjects
        s           = u_sub(ss);
        d           = t_data(ismember(t_sub,s),:);
        k           = t_kond(ismember(t_sub,s),:);
        %r           = t_raw(ismember(t_sub,s),:,:);
        r           = utils.filt_mova(t_raw(ismember(t_sub,s),:,:)',50,9)'; %raw MEP is filtered with a Gaussian Kernel (third argument==9), and a bandwidth of 50 samples (i.e. 10ms) second argument)        
        for knd = min(k):max(k),
           M(ssii,knd,ss)       = mean(d(k==knd));
           RW(ssii,knd,ss,:)    = mean(utils.baseline(r(k==knd,:),1:100,1)); %trials are baselined to the average of the 20ms prior to stimulation
           t_shifted            = utils.baseline(r(k==knd,:),1:100,1); %trials are baselined to the average of the 20ms prior to stimulation
           for e=1:size(t_shifted,1),
               % we calculate the latency as the average between negative
               % and positive peak, and shift the MEP accordingly
               [~,idx]  = sort(t_shifted(e,150:400));
               t_lat    = int32(150+mean([idx(1),idx(end)]));
               shiftMEP(ssii,knd,ss,:) = (t_shifted(e,t_lat-100:t_lat+99));
           end
        end
    end

end
%% Visual inspection of descriptive statistics
%% mep parameters
BI  = cat(2,nanmean(nanmean(M(:,setup.IO.BI==1,:),3),2),nanmean(nanmean(M(:,setup.IO.BI==0,:),3),2));
LM  = cat(2,nanmean(nanmean(M(:,setup.IO.LM==1,:),3),2),nanmean(nanmean(M(:,setup.IO.LM==0,:),3),2));
M1  = cat(2,nanmean(nanmean(M(:,setup.IO.M1==1,:),3),2),nanmean(nanmean(M(:,setup.IO.M1==0,:),3),2));
plot(M1)
plot(BI)
plot(LM)

%% raw waveform
rmtOI = 4;
BI  = squeeze(cat(2,nanmean(nanmean(RW(:,setup.IO.BI==1,:,:),3),2),nanmean(nanmean(RW(:,setup.IO.BI==0,:,:),3),2)));
LM  = squeeze(cat(2,nanmean(nanmean(RW(:,setup.IO.LM==1,:,:),3),2),nanmean(nanmean(RW(:,setup.IO.LM==0,:,:),3),2)));
M1  = squeeze(cat(2,nanmean(nanmean(RW(:,setup.IO.M1==1,:,:),3),2),nanmean(nanmean(RW(:,setup.IO.M1==0,:,:),3),2)));
plot(squeeze(M1(rmtOI,:,:))')
plot(squeeze(BI(rmtOI,:,:))')
plot(squeeze(LM(rmtOI,:,:))')

%% latency shifted
rmtOI = 4;
BI  = squeeze(cat(2,nanmean(nanmean(shiftMEP(:,setup.IO.BI==1,:,:),3),2),nanmean(nanmean(shiftMEP(:,setup.IO.BI==0,:,:),3),2)));
LM  = squeeze(cat(2,nanmean(nanmean(shiftMEP(:,setup.IO.LM==1,:,:),3),2),nanmean(nanmean(shiftMEP(:,setup.IO.LM==0,:,:),3),2)));
M1  = squeeze(cat(2,nanmean(nanmean(shiftMEP(:,setup.IO.M1==1,:,:),3),2),nanmean(nanmean(shiftMEP(:,setup.IO.M1==0,:,:),3),2)));
plot(squeeze(M1(rmtOI,:,:))')
plot(squeeze(BI(rmtOI,:,:))')
plot(squeeze(LM(rmtOI,:,:))')
%% 





















%% statistical analysis
%[p,tab,stats] = anovan(ampDATA,cat(2,rmtINT,KOND,SUB),'random',3,'display','off','varnames',{'SI','Condition','Subject'},'model',[1 0 0;0 1 0;0 0 1;1 1 0]);
%[p,tab,stats] = anovan(ampDATA,cat(2,rmtINT,DESIGN,SUB),'random',5,'display','off','varnames',{'SI','Biphasic','M1','90°','Subject'},'model',[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;1 1 1 1 0]);

M   = [];
E   = [];
BI  = [];
LM  = [];
M1  = [];
for ssii=1:length(setup.IO.SI)
    si              = setup.IO.SI(ssii);    
    [p,tab,stats]   = anovan(ampDATA(ismember(rmtINT,si)),cat(2,DESIGN(ismember(rmtINT,si),:),SUB(ismember(rmtINT,si))),'display','off','varnames',{'Biphasic','90°','M1','Subject'},'model',[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;1 1 1 0]);
    [c,m,h,nms]     = multcompare(stats,'dim',[1,2,3],'display','off');
    E               = cat(2,E,m(:,2));    
    
    [c,m,h,nms]     = multcompare(stats,'dim',[1],'display','off');
    BI              = cat(2,BI,m(:,1));
        
    [c,m,h,nms]     = multcompare(stats,'dim',[2],'display','off');
    LM              = cat(2,LM,m(:,1));

    [c,m,h,nms]     = multcompare(stats,'dim',[3],'display','off');
    M1              = cat(2,M1,m(:,1));

end


for c=1:length(setup.IO.label.all),
    
    
end







for ssii=1:length(setup.IO.SI)
    si          = setup.IO.SI(ssii);
    datapool    = ampDATA(ismember(rmtINT,si));
    subpool     = SUB(ismember(rmtINT,si));    
    
    
    
end



