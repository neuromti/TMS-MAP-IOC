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

for i_d = 1:length(Label.Dataset)
    % Where to save the files
    savefile                    = [folder.results.stats,Label.Dataset{i_d},'\stats.mat'];
    [savepath,filename,filext]  = fileparts(savefile);
    
    % Initalize for Logging
    log_betreff   = [Label.Field{(ismember({'AMP','MEP','LAT'},Label.Dataset{i_d}(end-2:end)))},' ',Label.Weight{(ismember({'qu','th'},Label.Dataset{i_d}(1:2)))},': Statistics '];
    
    if ~exist (savefile,'file'), %check whether already processed before 
        betreff = [log_betreff,'Started'];
        disp(betreff)
        eval(['DATA = cat(2,GRD.',Label.Dataset{i_d},');']);
        
        A = datetime('now');
        [TestResults,ClusterResults] = utils.mappingdata2statistics(DATA,SUBID,DESIGN);     
        O = datetime('now')-A;
        O = sprintf('It ran for %.1g years %.1g months %.1g days %.1g hours %.2g minutes %.2g seconds',datevec(O));
        if ~exist(savepath,'dir'), mkdir(savepath); end
        save(savefile,'TestResults','ClusterResults');   
        
        betreff = [log_betreff,'Finished'];
        notify_me(betreff,O);
    else
        betreff = [log_betreff,'already analyzed'];
        notify_me(betreff,'');

    end
end
%% VISUALIZATION
addpath('C:\Users\Robert Bauer\Documents\MATLAB\other_toolboxes\CETperceptual_MATLAB'); %folder with colormaps
cmap = diverging_bwr_40_95_c42_n256;
set(groot,'DefaultFigureColormap',cmap)
clc

% Show Effects of Factors on the different measures on the grid 
% and project them on a headmodel

for i_d = 1:length(Label.Dataset)
    load([folder.results.stats,'\',Label.Dataset{i_d},'\stats.mat']);            
    i_field     = find(ismember({'AMP','MEP','LAT'},Label.Dataset{i_d}(end-2:end)));
    i_weight    = find(ismember({'qu','th'},Label.Dataset{i_d}(1:2)));
  
    close all
    delete([folder.results.figs,Label.Weight{i_w},'\*.tif'])
    for k=1:size(TestResults.Pval,2)
        

        PlotLabel   = [Label.Field{i_field},' ',Label.Weight{i_weight},' ',Label.Title{k}];
        
        % P           = -log10(TestResults.Pval(:,k));
        
        P           = -log10(TestResults.PermutationPval(:,k));
        S           = sign(TestResults.Coeffs(:,k));
        V           = P.*S;
        
        utils.plot_gridmodel(V,2)
        title(Label.Title{k})             
        print(gcf,[folder.results.figs,Label.Weight{i_w},'\GRID_',Label.Save{k},'-',Label.Field{i_field},'.tif'],'-dtiff')        
        
        utils.plot_headmodel(headmodel,V,2)
        annotation('textbox','Position',[0 0 1 1],'String',PlotLabel)
        print(gcf,[folder.results.figs,Label.Weight{i_w},'\HEAD_',Label.Save{k},'-',Label.Field{i_field},'.tif'],'-dtiff')

        
        for clus_idx = 1:length(ClusterResults(k).PermPval),
            Grid    = utils.get_DesignGrid;
            sel     = utils.mesh2vec(ClusterResults(k).MemberShip(:,:,1)==clus_idx);
            xyz     = double(Grid(sel,:));
            V       = utils.Target2Grid(xyz);

            %utils.plot_headmodel(headmodel,V,1)

            utils.plot_gridmodel(V,2)
            title(PlotLabel)             
            colorbar off
            print(gcf,[folder.results.figs,Label.Weight{i_w},'\GRID_',Label.Save{k},'-',Label.Field{i_field},'-Cluster_',num2str(clus_idx),'.tif'],'-dtiff')
            
            utils.plot_headmodel(headmodel,V,2)
            annotation('textbox','Position',[0 0 1 1],'String',PlotLabel)
            print(gcf,[folder.results.figs,Label.Weight{i_w},'\HEAD_',Label.Save{k},'-',Label.Field{i_field},'-Cluster_',num2str(clus_idx),'.tif'],'-dtiff')
            
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
