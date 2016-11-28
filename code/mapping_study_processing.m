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
clear SUB GRD ALGN
fprintf('[')
logfilename = [folder.code,'Logfile.log'];
logfileid   = fopen(logfilename,'wt');
fprintf(logfileid,'started script on %s \n',datetime('now'));
for idx_dataset = 1:length(D),    
          
    % Load dataset
    filename    = D(idx_dataset).name;
    load([folder.data.map,filename]);     
    fprintf('-')    
    
    % Arrange and set index variables, connect conditons and location with data     
    sub                             = utils.intialize_dataset(filename,setup);
    
    % Remove trials with artifacts, calculate MEP+ and center latency,
    [mapping,sub]                    = utils.prepare_data(mapping,sub);      
    
    % Estimate or Get Points of Interest for the Subject 
    sub.CoG                         = utils.calculate_CoG(sub);
    sub.HS                          = utils.get_HS(setup,sub);
    sub.ANT                         = utils.get_ANT(setup,sub);
    sub.GridOrigin                  = utils.get_GridOrigin(mapping);
    SUB(idx_dataset)                = sub;
    fprintf('\b%s','|')
         
    % Put on Predefined Rectangular grid  
    GRD(idx_dataset)                        = utils.project_on_Grid(mapping);
    fprintf('\b%s','-')       
        
    % Log
    fprintf('\b%s','.')
    fprintf(logfileid,'At %s finished dataset %i \n',datetime('now'),idx_dataset);
end
fprintf(']\n')
fclose(logfileid);
clear logfileid;
save([folder.results.stats,'map_subject_data.mat'],'SUB')
save([folder.results.stats,'map_gridded.mat'],'GRD')
%%
% GROUP PARAMETERS
load([folder.results.stats,'map_gridded.mat'],'GRD')
load([folder.results.stats,'map_subject_data.mat'],'SUB')


SUBID   = cat(1,SUB.subID);
[~,sort_idx] = sort(SUBID);
DESIGN      = cat(1,SUB(sort_idx).DesignMatrix);
SUBID       = SUBID(sort_idx);

FieldList = fieldnames(GRD);
for i_f = 1:length(FieldList)
    eval(['DATA = cat(3,GRD(sort_idx).',FieldList{i_f},');']);
    [true_M,true_P,perm_P,bot_CI,bot_P] = utils.do_botandperm_MAP(DATA,SUBID,DESIGN);
    if ~exist([folder.results.stats,'\',FieldList{i_f},'\'],'dir'),mkdir([folder.results.stats,'\',FieldList{i_f},'\']); end
    save([folder.results.stats,'\',FieldList{i_f},'\nonparametric.mat'],'true_M','true_P','perm_P','bot_CI','bot_P');    
    
end
%%

%     FaceColor = utils.grid2headmodel(headmodel,-log10(true_P(:,1))); 
%     figure
%     ft_plot_mesh(headmodel,'edgealpha',0,'facecolor',FaceColor)
%     viewpoint = [-88 60];
%     view(viewpoint)
%    % annotation('textbox','Position',[0 0 1 1],'String',anno_label{k})
%     hnd_cb = colorbar;
%     yticklab                    = fix(100*linspace(-2,2,11))./100;
%     set(hnd_cb,'YLIM',[0 1],'YTICK',[0:.1:1],'YTICKLABEL',yticklab)
% 
% % PLOT ON GYRI HEADMODEL
% close all
% for k=1:length(anno_label),
%     
% 
%  
% 
%     %yticklab                    = fix(100*linspace(-tmp_max_abs,tmp_max_abs,11))./100;
%     
%     light('Position',-diag(rotx(viewpoint(1))*roty(viewpoint(2))))
%     material dull       	
% %    lighting gouraud
% end
% % 
% % 
% cmap                        = diverging_bwr_40_95_c42_n256;
% %cmap                        = calc_colormap(2,1.3/2,'white',256);
% set(groot,'DefaultFigureColormap',cmap)
%       
% 
% % anno_label = {'BI=0,LM=0'    'BI=1,LM=0'    'BI=0,LM=1'    'BI=1,LM=1'}
% % anno_label = {[setup.MAP.label.BI{1},' > ',setup.MAP.label.BI{2}]...
% %             [setup.MAP.label.LM{1},' > ',setup.MAP.label.LM{2}],...
% %             [setup.IO.label.ORG{1},' X ',setup.IO.label.ORG{2}],...
% %             }
% % 
% % 
