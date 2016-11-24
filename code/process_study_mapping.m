%% Input-Output-Curve Analysis
% configuration
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
cls; %clc, clear all, close all, matlabrc, fclose('all'), addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils

addpath('C:\Users\Robert Bauer\Documents\MATLAB\other_toolboxes\CETperceptual_MATLAB'); %folder with colormaps
cmap        = linear_kry_5_95_c72_n256;
set(groot,'DefaultFigureColormap',cmap)

load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','headmodel','setup','folder');

%% loading data
logfilename = [folder.code,'Logfile.log'];
logfileid   = fopen(logfilename,'wt');
fprintf(logfileid,'started script on %s \n',datetime('now'));

% List the .mat files int the data directory
D           = dir([folder.data.map,'\*.mat']);

% Find artifacted datasets 
IsArtifacted        = utils.check_artifacts(folder.data.map);
D(IsArtifacted)     = [];

% Define the query grid based on surface mesh-grid 
sfc                 = load('C:\Users\Robert Bauer\Documents\Matlab\other_toolboxes\fieldtrip\template\anatomy\surface_white_left.mat');
clear queried
queried.pos         = sfc.bnd.pnt;
queried.tri         = sfc.bnd.tri;
%%
clear IPOL 
for idx_dataset = 1:length(D),    
       
    % Load dataset
    load([folder.data.map,D(idx_dataset).name]);
     
    % Arrange and set index variables
    % Glue conditons and location to data    
    t.s                     = regexp(D(idx_dataset).name,'S\w*');
    t.c                     = regexp(D(idx_dataset).name,'C\w*');
    t.d                     = regexp(D(idx_dataset).name,'.mat');
    t.data_sub              = int32(str2double(D(idx_dataset).name(t.s+1:t.c-1)));
    t.data_cond             = int32(str2double(D(idx_dataset).name(t.c+1:t.d-1)));
    t.sub                   = find(setup.SUB.id==t.data_sub);    
    t.sex                   = setup.SUB.sex(t.sub);
    t.age                   = setup.SUB.age(t.sub);
    
    IPOL(idx_dataset).HS        = [setup.SUB.pos{t.sub}{t.data_cond+1,2:4}];
    IPOL(idx_dataset).ANT       = [setup.SUB.pos{t.sub}{t.data_cond+1,5:7}];
    IPOL(idx_dataset).COND      = setup.MAP.label.all{t.data_cond};
    IPOL(idx_dataset).DESIGN    = logical([setup.MAP.BI(t.data_cond),setup.MAP.LM(t.data_cond)]);
    IPOL(idx_dataset).sub_sex   = t.sex;
    IPOL(idx_dataset).sub_age   = t.age;
    IPOL(idx_dataset).sub_idx   = t.data_sub;
    clear t

    
    % Remove trials with artifacts, calculate MEP+ and center latency,
    sub                         = utils.prepare_for_weighting(mapping);    
    % Interpolate individual subjects MEP+ and latency unto the surface based on squared distance    
    [quad_weights,threshold_weights]    = utils.calculate_distanceweights(sub.xyz,queried.pos);

    IPOL(idx_dataset).quad_AMP      = quad_weights*sub.amp;
    IPOL(idx_dataset).quad_MEP      = quad_weights*sub.mep;     
    IPOL(idx_dataset).quad_LAT      = quad_weights*sub.lat;     
    
    IPOL(idx_dataset).thresh_AMP    = threshold_weights*sub.amp;
    IPOL(idx_dataset).thresh_MEP    = threshold_weights*sub.mep;     
    IPOL(idx_dataset).thresh_LAT    = threshold_weights*sub.lat;     

    % Interpolate as above, but based on CoG shifted to M1   
    sub                         = utils.calculate_CoG(sub);
    sub                         = utils.shift_by_CoG(sub);
    [shifted_quad_weights,shifted_threshed_weights] = utils.calculate_distanceweights(sub.shifted_xyz,queried.pos);    

    IPOL(idx_dataset).CoG                   = sub.CoG;
    IPOL(idx_dataset).shifted_quad_AMP      = shifted_quad_weights*sub.amp;
    IPOL(idx_dataset).shifted_quad_MEP      = shifted_quad_weights*sub.mep;     
    IPOL(idx_dataset).shifted_quad_LAT      = shifted_quad_weights*sub.lat;   
    
    IPOL(idx_dataset).shifted_thresh_AMP    = shifted_threshed_weights*sub.amp;
    IPOL(idx_dataset).shifted_thresh_MEP    = shifted_threshed_weights*sub.mep;     
    IPOL(idx_dataset).shifted_thresh_LAT    = shifted_threshed_weights*sub.lat;     
    
    % Interpolate as above, but based on CoG shifted to M1 and main axis aligned towards anterior-posterioor     ´    
    aligned_quad_weights        = utils.calculate_distanceweights(sub.aligned_xyz,queried.pos,25);    

    IPOL(idx_dataset).CoG               = sub.CoG;
    IPOL(idx_dataset).aligned_quad_AMP  = aligned_quad_weights*sub.amp;
    IPOL(idx_dataset).aligned_quad_MEP  = aligned_quad_weights*sub.mep;     
    IPOL(idx_dataset).aligned_quad_LAT  = aligned_quad_weights*sub.lat;     
    
    fprintf(logfileid,'At %s finished dataset %i \n',datetime('now'),idx_dataset);
end
fclose(logfileid);
save([folder.results.stats,'map_interpolated.mat'],'IPOL')
% 
%% PERMUTATION ANALYSIS
% GROUP PARAMETERS
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