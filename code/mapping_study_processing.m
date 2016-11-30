%% Configuration
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
cls; %clc, clear all, close all, matlabrc, fclose('all'), addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils

addpath('C:\Users\Robert Bauer\Documents\MATLAB\other_toolboxes\CETperceptual_MATLAB'); %folder with colormaps
cmap        = linear_kry_5_95_c72_n256;
set(groot,'DefaultFigureColormap',cmap)

load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','gridmodel','headmodel','setup','folder');

%% Loading data
% List the .mat files int the data directory
D           = dir([folder.data.map,'\*.mat']);

% Find artifacted datasets 
IsArtifacted        = utils.find_ArtifactedSubjects(folder.data.map);
D(IsArtifacted)     = [];
%%
clc
clear SUB GRD
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
utils.progressBar(']')
fclose(logfileid);
clear logfileid;
save([folder.results.stats,'map_subject_data.mat'],'SUB')
save([folder.results.stats,'map_gridded.mat'],'GRD')
%%
% GROUP PARAMETERS
load([folder.results.stats,'map_subject_data.mat'],'SUB')
load([folder.results.stats,'map_gridded.mat'],'GRD')

SUBID   = cat(1,SUB.subID);
[~,sort_idx] = sort(SUBID);
DESIGN      = cat(1,SUB(sort_idx).DesignMatrix);
SUBID       = SUBID(sort_idx);

cmap                        = diverging_bwr_40_95_c42_n256;
%cmap                        = calc_colormap(2,1.3/2,'white',256);
%cmap                        = calc_colormap(2,0,'white',24);
set(groot,'DefaultFigureColormap',cmap)
anno_label = {[setup.MAP.label.BI{1},' > ',setup.MAP.label.BI{2}]...
            [setup.MAP.label.LM{1},' > ',setup.MAP.label.LM{2}],...
            [setup.MAP.label.BI{1},' ',setup.MAP.label.LM{1},' & ',setup.MAP.label.BI{2},' ',setup.MAP.label.LM{2}],...
            };
        
SaveLabel  =  {'Waveform', 'Orientation', 'W-O-Interaction'};
FieldLabel = {'Amplitude','MEP','Latency'};
WeightLabel = {'Quadratic','Spherical'};   
DataList = {'quad_AMP'    'quad_MEP'    'quad_LAT'    'th_AMP'    'th_MEP'    'th_LAT'};
FieldList = fieldnames(GRD);

for i_f = 1:length(FieldList)
    
    eval(['DATA = cat(3,GRD(sort_idx).',FieldList{i_f},');']);
    tic
    [TestResults,ClusterResults] = utils.mappingdata2statistics(DATA,SUBID,DESIGN);     
    toc
    if ~exist([folder.results.stats,'\',FieldList{i_f},'\'],'dir'),mkdir([folder.results.stats,'\',FieldList{i_f},'\']); end
    save([folder.results.stats,'\',FieldList{i_f},'\stats.mat'],'TestResults','ClusterResults');        
    
    notify_me([FieldLabel{find(ismember({'AMP','MEP','LAT'},DataList{i_f}(end-2:end)))},' ',WeightLabel{find(ismember({'qu','th'},DataList{i_f}(1:2)))},'Statistics Finished'],'Grid Interpolation');
end
%%
clc
for i_d = 1:length(DataList)
    load([folder.results.stats,'\',FieldList{i_d},'\stats.mat']);            
    i_f = find(ismember({'AMP','MEP','LAT'},DataList{i_d}(end-2:end)));
    i_w = find(ismember({'qu','th'},DataList{i_d}(1:2)));
    TitleLab = FieldLabel{i_f};
    
    for k=1:length(anno_label),
        
        PlotLabel   = [TitleLab,' ',anno_label{k}];

        %P           = -log10(TestResults.Pval(:,k));
        P           = -log10(TestResults.PermutationPval(:,k));
        S           = sign(TestResults.Coeffs(:,k));
        V           = P.*S;
        
        utils.plot_gridmodel(V,2)
        title(PlotLabel)             
        print(gcf,[folder.results.figs,WeightLabel{i_w},'\GRID_',SaveLabel{k},'-',TitleLab,'.tif'],'-dtiff')
        
        
        utils.plot_headmodel(headmodel,V,2)
        annotation('textbox','Position',[0 0 1 1],'String',PlotLabel)
        print(gcf,[folder.results.figs,WeightLabel{i_w},'\HEAD_',SaveLabel{k},'-',TitleLab,'.tif'],'-dtiff')
        close all
    end
       
end
