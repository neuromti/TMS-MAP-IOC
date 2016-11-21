%% Input-Output-Curve Analysis
% configuration
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
cls; %clc, clear all, close all, matlabrc, fclose('all'), addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','headmodel','setup','folder');

%% loading data
logfilename = [folder.code,'Logfile.log'];
logfileid   = fopen(logfilename,'wt');

% read in the datafiles int the directory
D           = dir(folder.data.map);
D([1,2])    = [];

% find artifacted datasets 
IsArtifacted = utils.check_artifacts(folder.data.map);

% Define the outlines of the group-wise grid for projection of individual
% maps
[grid_xyz,map_dimensions,X,Y,Z] = utils.make_groupgrid(folder.data.map,IsArtifacted);

%define some parameters and anonymous functions
do_norm         = @(x)x/norm(x,2);    
accepted_dist   = 15; %1.5 cm distance



clear IPOL 

for idx_dataset = 1:length(D),    
       
    % arrange and set index variables
    s               = regexp(D(idx_dataset).name,'S\w*');
    c               = regexp(D(idx_dataset).name,'C\w*');
    d               = regexp(D(idx_dataset).name,'.mat');
    idx_data_sub    = int32(str2double(D(idx_dataset).name(s+1:c-1)));
    idx_data_cond   = int32(str2double(D(idx_dataset).name(c+1:d-1)));
    idx_sub         = find(setup.SUB.id==idx_data_sub);
       
    % load dataset if it is not artifacted
    if artifacted{idx_dataset},
        continue;
    end
    load([folder.data.map,D(idx_dataset).name]);

    
    collect_XYZ     = NaN(size(mapping.xyz,1),size(grid_xyz,1));
    all_dist        = pdist2(mapping.xyz,grid_xyz);
    t_accepted      = ~(all_dist>accepted_dist);
    t_not_incl      = (all(~t_accepted,1));
    i_not_incl      = find(~t_not_incl);
    for idx_coordinate=i_not_incl,
             
        t_dist      = all_dist(:,idx_coordinate);        
        t_invdist   = 1./t_dist;
        t_remv      = ~t_accepted(:,idx_coordinate);
        t_invdist(t_remv)=[];
        t_weight    = do_norm((t_invdist));
        t_interpol  = t_weight'*mapping.amp(~t_remv);
        collect_XYZ(find(~t_remv),idx_coordinate) = t_interpol;
   
    end
    

    
    f_interpolated =  (nanmean(collect_XYZ,1));
    IPOL(idx_dataset).DATA  = f_interpolated;
    IPOL(idx_dataset).SUB   = idx_sub;
    IPOL(idx_dataset).COND  = idx_data_cond;
    % mesh(X(:,:,1),Y(:,:,1),inter_XYZ(:,:,1))
    
end   
    
 cond   = cat(1,IPOL.COND); 
 data   = cat(1,IPOL.DATA);
 sub    = cat(1,IPOL.SUB); 
    
mapped_data  = (reshape(nanmean(data(cond==1,:),1),map_dimensions(1),map_dimensions(2),map_dimensions(3)));
mesh(X(:,:,1),Y(:,:,1),mapped_data(:,:,1))