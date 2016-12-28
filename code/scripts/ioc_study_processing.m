%% Input-Output-Curve Analysis
% configuration;
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
cls; %clc, clear, close all, matlabrc, fclose('all'), addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','setup','folder','Label');
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');

%% loading data
D           = dir([folder.data.ioc,'*.mat']);
clc
% Preparing data structure
SUB         = struct('AmplitudeAcrossIntensity',[],'MepAcrossIntensity',[],'LatencyAcrossIntensity',[],'WaveformAcrossIntensity',[],'SUBID',[],'CONDITION',[]);
% load data from files and concatenate
for idx_dataset = 1 : length(D)

    % Load dataset
    filename    = D(idx_dataset).name;
    load([folder.data.ioc,filename]);

    % Get Values
    SUB(idx_dataset).LatencyAcrossIntensity      = utils.get_LatencyAcrossIntensity(ioc);
    SUB(idx_dataset).AmplitudeAcrossIntensity    = utils.get_AmplitudeAcrossIntensity(ioc);
    SUB(idx_dataset).MepAcrossIntensity          = utils.get_MepAcrossIntensity(ioc);
    SUB(idx_dataset).TimeCourseAcrossIntensity   = utils.get_TimeCourseAcrossIntensity(ioc);

    % Get Condition
    T                                           = utils.scan_DataFileName(filename,setup);
    SUB(idx_dataset).SUBID                      = repmat(T.subID,7,1);
    tmp_DesignMatrix                            = [setup.IO.BI(T.data_cond),setup.IO.LM(T.data_cond),setup.IO.M1(T.data_cond)];
    SUB(idx_dataset).DESIGN                     = repmat(tmp_DesignMatrix,7,1);
    SUB(idx_dataset).CONDITION                  = bin2dec(num2str(tmp_DesignMatrix));
    SUB(idx_dataset).STIMINTENSITYinMSO         = sort(unique(ioc.int),'ascend');
    SUB(idx_dataset).STIMINTENSITYinSTEPS       = [1:7]';

end
save([folder.results.stats,'ioc_data.mat'],'SUB')
%% PERMUTATION ANALYSIS

load([folder.results.stats,'ioc_data.mat'],'SUB')
DESIGN  = cat(1,SUB.DESIGN);
SUBID   = cat(1,SUB.SUBID);
STIM    = cat(1,SUB.STIMINTENSITYinSTEPS);

FieldList = {'AmplitudeAcrossIntensity','MepAcrossIntensity','LatencyAcrossIntensity','TimeCourseAcrossIntensity'};
for i_d = 1:length(FieldList)
     % Initalize for Logging and Saving
    savefile            = [folder.results.stats,'ioc\',Label.ioc_Field{i_d},'_stats.mat'];
    if ~exist(fileparts(savefile),'dir'), mkdir(fileparts(savefile)); end
    log_betreff   = [Label.ioc_Field{i_d},': Statistics '];

    if exist (savefile,'file') %check whether already processed before
        notify_me([log_betreff,'already analyzed'],'');
    else
        % Read Data from Structure
        eval(['DATA    = cat(1,SUB.',FieldList{i_d},');'])
        betreff = [log_betreff,'Started'];
        disp(betreff)

        % Run Analysis and save results
        ProcessingTime = utils.measure_ProcessingTime();

        if strcmpi(Label.ioc_Field{i_d},'TimeCourse')
            TestResults  = utils.TimeCourse2statistics(DATA,DESIGN,SUBID,STIM);
            save(savefile,'TestResults');
        else
            [TestResults,ClusterResults] = utils.iocdata2statistics(DATA,DESIGN,SUBID,STIM);
            save(savefile,'TestResults','ClusterResults');
        end

        notify_me([log_betreff,'Finished'],utils.measure_ProcessingTime(ProcessingTime));
    end
end
%%
