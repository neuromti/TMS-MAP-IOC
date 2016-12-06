%% Configuration
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
cls; %clc, clear, close all, matlabrc, fclose('all'), addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','gridmodel','headmodel','setup','folder');
%% Loading data
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
%% STATISTICAL ANALYSIS
% GROUP PARAMETERS
clc, clear, close all
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','gridmodel','headmodel','setup','folder','Label');
load([folder.results.stats,'map_subject_data.mat'],'SUB')
load([folder.results.stats,'map_gridded.mat'],'GRD')

DESIGN          = cat(1,SUB.DesignMatrix);
SUBID           = cat(1,SUB.subID);

delete(gcp('nocreate'));
parpool(3)
parfor i_d = 1:length(length(Label.Dataset(1:3)))
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
delete(gcp('nocreate'));
%% VISUALIZATION
addpath('C:\Users\Robert Bauer\Documents\MATLAB\other_toolboxes\CETperceptual_MATLAB'); %folder with colormaps
cmap = diverging_bwr_40_95_c42_n256;
set(groot,'DefaultFigureColormap',cmap)
clc

% Show Effects of Factors on the different measures on the grid 
% and project them on a headmodel
delete([folder.results.figs,Label.Weight{1},'\*.tif'])
delete([folder.results.figs,Label.Weight{2},'\*.tif'])

for i_d = 1:length(Label.Dataset(1:3))
    load([folder.results.stats,'\',Label.Dataset{i_d},'\stats.mat']);            
    i_field     = find(ismember({'AMP','MEP','LAT'},Label.Dataset{i_d}(end-2:end)));
    i_weight    = find(ismember({'qu','th'},Label.Dataset{i_d}(1:2)));
  
    close all
    
    for k=1:size(TestResults.Pval,2)        

        PlotLabel   = [Label.Title{k},' ->',Label.Field{i_field}];
         
        P           = -log10(TestResults.PermutationPval(:,k));
        S           = sign(TestResults.Coeffs(:,k));
        V           = P.*S;
        
        utils.plot_gridmodel(V,2)
        title(PlotLabel)             
        print(gcf,[folder.results.figs,Label.Weight{i_weight},'\GRID_',Label.Save{k},'-',Label.Field{i_field},'.tif'],'-dtiff')        
        
        utils.plot_headmodel(headmodel,V,2)
        annotation('textbox','Position',[0 0 1 1],'String',PlotLabel)
        print(gcf,[folder.results.figs,Label.Weight{i_weight},'\HEAD_',Label.Save{k},'-',Label.Field{i_field},'.tif'],'-dtiff')

        if ~isempty(ClusterResults(k).PermPval)
            for clus_idx = 1:length(ClusterResults(k).PermPval),
                Grid    = utils.get_DesignGrid;
                sel     = utils.mesh2vec(ClusterResults(k).MemberShip(:,:,1)==clus_idx);
                xyz     = double(Grid(sel,:));
                V       = utils.Target2Grid(xyz);

                %utils.plot_headmodel(headmodel,V,1)
                ClusterAddendum     = sprintf('(Cluster %i: p=%.2g )',clus_idx,ClusterResults(k).PermPval(clus_idx));
                ClusterLabel        = [PlotLabel,' ',ClusterAddendum];

                utils.plot_gridmodel(V,2)
                title(ClusterLabel)   
                colorbar off
                print(gcf,[folder.results.figs,Label.Weight{i_weight},'\CLUSTER_GRID_',Label.Save{k},'-',Label.Field{i_field},'-',num2str(clus_idx),'.tif'],'-dtiff')

                utils.plot_headmodel(headmodel,V,2)
                annotation('textbox','Position',[0 0 1 1],'String',ClusterLabel)
                print(gcf,[folder.results.figs,Label.Weight{i_weight},'\CLUSTER_HEAD_',Label.Save{k},'-',Label.Field{i_field},'-',num2str(clus_idx),'.tif'],'-dtiff')

            end                
        end
    end
    
    
end

% Plot Position of Anterior Target
close all
ANT         = utils.get_GroupAnt(SUB);
V           = utils.Target2Grid(ANT);
utils.plot_headmodel(headmodel,V,1)
annotation('textbox','Position',[0 0 1 1],'String','Distribution of Anterior Target')
print(gcf,[folder.results.figs,'\HEAD_PDF_Anterior_Target.tif'],'-dtiff')

utils.plot_gridmodel(V,2)
title('Distribution of Anterior Target')
print(gcf,[folder.results.figs,'\GRID_PDF_Anterior_Target.tif'],'-dtiff')
%
