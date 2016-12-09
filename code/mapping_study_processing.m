%% Configuration
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
utils.cls; %clc, clear, close all, matlabrc, fclose('all'), addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','gridmodel','headmodel','setup','folder');
%% Loading data
profile on -timer 'performance'

% List the .mat files int the data directory
D           = dir([folder.data.map,'\*.mat']);

% Find artifacted datasets 
IsArtifacted        = utils.find_ArtifactedSubjects(folder.data.map);
D(IsArtifacted)     = [];

logfilename = [folder.code,'Logfile.log'];
logfileid   = fopen(logfilename,'wt');
fprintf(logfileid,'started script on %s \n',datetime('now'));
utils.progressBar('[')
for idx_dataset = 1:length(D),    
    % Load dataset
    filename    = D(idx_dataset).name;
    load([folder.data.map,filename]);     
    % Arrange and set index variables, connect conditons and location with data     
    sub                             = utils.intialize_dataset(filename,setup);
    
    % Remove trials with artifacts, calculate MEP+ and center latency,
    mapping                         = utils.prepare_data(mapping);      
    
    % Estimate or Get Points of Interest for the Subject 
    sub.CoG                         = utils.calculate_CoG(mapping);
    sub.HS                          = utils.get_HS(setup,sub);
    sub.ANT                         = utils.get_ANT(setup,sub);
    sub.GridOrigin                  = utils.get_GridOrigin(mapping);
    SUB(idx_dataset)                = sub;
      
    % Put on Predefined Rectangular grid  
    GRD(idx_dataset)                = utils.project_on_Grid(mapping);

    % Log
    fprintf(logfileid,'At %s finished dataset %i \n',datetime('now'),idx_dataset);
    utils.progressBar('.')
end
utils.progressBar('1')
fclose(logfileid);
clear logfileid;

% Sort data according to subject ID
[~,sort_idx]    = sort(cat(1,SUB.subID));
SUB = SUB(sort_idx);
GRD = GRD(sort_idx);

% Save the processed data
save([folder.results.stats,'map_subject_data.mat'],'SUB')
save([folder.results.stats,'map_gridded.mat'],'GRD')

profile off
profile viewer

%% STATISTICAL ANALYSIS
% GROUP PARAMETERS
clc, clear, close all
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','gridmodel','headmodel','setup','folder','Label');
load([folder.results.stats,'map_subject_data.mat'],'SUB')
load([folder.results.stats,'map_gridded.mat'],'GRD')

DESIGN          = cat(1,SUB.DesignMatrix);
SUBID           = cat(1,SUB.subID);

for i_d = 1:length(Label.Dataset(1:3))
    % Where to save the files
    savefile                    = [folder.results.stats,Label.Dataset{i_d},'\stats.mat'];
    if ~exist(fileparts(savefile),'dir'), mkdir(fileparts(savefile)); end
    log_betreff   = [Label.Field{(ismember({'AMP','MEP','LAT'},Label.Dataset{i_d}(end-2:end)))},' ',Label.Weight{(ismember({'qu','th'},Label.Dataset{i_d}(1:2)))},': Statistics '];
    
    if exist (savefile,'file'), %check whether already processed before
        notify_me([log_betreff,'already analyzed'],'');
    else
        disp([log_betreff,'Started']);
        eval(['DATA = cat(2,GRD.',Label.Dataset{i_d},');']);
        
        ProcessingTime = utils.measure_ProcessingTime();
        [TestResults,ClusterResults] = utils.mappingdata2statistics(DATA,SUBID,DESIGN);     
        save(savefile,'TestResults','ClusterResults');      
        notify_me([log_betreff,'Finished'],utils.measure_ProcessingTime(ProcessingTime));
    end
end



