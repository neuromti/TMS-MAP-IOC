%% Input-Output-Curve Analysis
% configuration
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
cls; %clc, clear all, close all, matlabrc, fclose('all'), addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils

addpath('C:\Users\Robert Bauer\Documents\MATLAB\other_toolboxes\CETperceptual_MATLAB'); %folder with colormaps
cmap        = linear_kry_5_95_c72_n256;
set(groot,'DefaultFigureColormap',cmap)

load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','gridmodel','headmodel','setup','folder');

%% loading data
logfilename = [folder.code,'Logfile.log'];
logfileid   = fopen(logfilename,'wt');
fprintf(logfileid,'started script on %s \n',datetime('now'));

% List the .mat files int the data directory
D           = dir([folder.data.map,'\*.mat']);

% Find artifacted datasets 
IsArtifacted        = utils.check_artifacts(folder.data.map);
D(IsArtifacted)     = [];

%%
clear IPOL SHFT ALGN
for idx_dataset = 1:length(D),    
       
    % Load dataset
    load([folder.data.map,D(idx_dataset).name]);
   
    sub = struct();
    % Arrange and set index variables    
    % Glue conditons and location to data     
    sub                             = utils.glue_dataset(D(idx_dataset).name,setup,sub);
    
    % Remove trials with artifacts, calculate MEP+ and center latency,
    sub                             = utils.prepare_data(mapping,sub);   

    % Estimate or Get Points of Interest for the Subject 
    sub.CoG                         = utils.calculate_CoG(sub);
    sub.HS                          = utils.get_HS(setup,sub);
    sub.ANT                         = utils.get_ANT(setup,sub);
    
    SUB(idx_dataset)                = sub;
    
    % Interpolate individual subjects MEP+ and latency unto the surface based on squared distance    
    [quad_weights,threshed_weights]         = utils.calculate_distanceweights(sub.xyz,headmodel.pos);
    IPOL(idx_dataset)                       = utils.perform_weighting(sub,quad_weights,threshed_weights);    
    clear quad_weights threshed_weights
    
    % Interpolate as above, but based on CoG shifted to M1   
    sub                                     = utils.shift_by_CoG(sub);
    [quad_weights,threshed_weights]         = utils.calculate_distanceweights(sub.shifted_xyz,headmodel.pos);    
    SHFT(idx_dataset)                       = utils.perform_weighting(sub,quad_weights,threshed_weights);   
    clear quad_weights threshed_weights
    
    % Interpolate as above, but based on CoG shifted to M1 and main axis aligned towards anterior-posterioor     ´    
    sub                                     = utils.normalize_by_HsAnt(sub,setup);        
    [quad_weights,threshed_weights]         = utils.calculate_distanceweights(sub.aligned_xyz,gridmodel.pos,.5);         
    ALGN(idx_dataset)                       = utils.perform_weighting(sub,quad_weights,threshed_weights);
    clear quad_weights threshed_weights
    
    % Log
    fprintf(logfileid,'At %s finished dataset %i \n',datetime('now'),idx_dataset);
end
fclose(logfileid);
save([folder.results.stats,'map_subject_data.mat'],'SUB')
save([folder.results.stats,'map_interpolated.mat'],'IPOL')
save([folder.results.stats,'map_shift_interpolated.mat'],'SHFT')
save([folder.results.stats,'map_normalized_interpolated.mat'],'ALGN')
%% PERMUTATION ANALYSIS
% GROUP PARAMETERS
% mean(grpstats(cat(1,IPOL.ANT),SUB))-mean(grpstats(cat(1,IPOL.HS),SUB));

load([folder.results.stats,'map_interpolated.mat'],'IPOL')
SUB     = cat(1,IPOL.sub_idx);
[~,sort_idx] = sort(SUB);
DESIGN  = cat(1,IPOL(sort_idx).DESIGN);
AMP     = cat(2,IPOL(sort_idx).quad_AMP);
LAT     = cat(2,IPOL(sort_idx).quad_LAT);
MEP     = cat(2,IPOL(sort_idx).quad_MEP);
SUB     = SUB(sort_idx);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORM STATISTICAL ANALYSIS FOR EACH VARIABLE OF INTEREST
num_rep     = 1000;

field_list = fieldnames(IPOL(1));
f_idx = [];
for f=1:length(field_list),
    f_idx = [f_idx,~isempty([regexp(field_list{f},'MEP'),regexp(field_list{f},'AMP'),regexp(field_list{f},'LAT')])];
end
field_list = {field_list{find(f_idx)}};

for field_idx = 1:length(field_list),
    eval(['DATA = cat(2,IPOL(sort_idx).',field_list{field_idx},');'])
    utils.do_botandperm_MAP(DATA,SUB,DESIGN,[folder.results.stats,field_list{field_idx},'\'],num_rep);
end