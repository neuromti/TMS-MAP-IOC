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
        
        %concatenate data vectors after cleaning for artifacts (indicated
        %by negative values) and lack of MEP (indicated by zero latency)       
        cleaned             = ioc.lat;
        artifacted          = (ioc.amp<0) | (ioc.lat <0 | ioc.lat == 0);
        cleaned(artifacted) = NaN;       
        t_lat               = nanmean(reshape(cleaned,10,7),1)';
        subAverage.latDATA  = cat(1,subAverage.latDATA,t_lat);
        
        %concatenate data vectors after cleaning for artifacts (indicated
        %by negative values)        
        cleaned             = ioc.amp;
        artifacted          = (ioc.amp<0) | (ioc.lat <0);
        cleaned(artifacted) = NaN;
        t_amp               = nanmean(reshape(cleaned,10,7),1)';
        t_mep               = nanmean(reshape(cleaned,10,7)>50,1)';
        subAverage.ampDATA  = cat(1,subAverage.ampDATA,t_amp);        
        subAverage.mepDATA  = cat(1,subAverage.mepDATA,t_mep); 
        
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
subAverage.Condition = bi2de(subAverage.Design(:,1:3))+1; 
%%













































