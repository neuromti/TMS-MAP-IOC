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

% read in the datafiles int the directory
D           = dir(folder.data.map);
D([1,2])    = [];

% find artifacted datasets 
IsArtifacted = utils.check_artifacts(folder.data.map);

% Define the outlines of the group-wise grid for projection of individual
% maps
[grid_xyz,map_dimensions,X,Y,Z] = utils.make_groupgrid(folder.data.map,IsArtifacted);
reshape_to_grid = @(x)reshape(x,map_dimensions(1),map_dimensions(2),map_dimensions(3));

clear IPOL 
for idx_dataset = 1:length(D),    
       
    % load dataset if it is not artifacted
    if IsArtifacted(idx_dataset),
        continue;
    end
    load([folder.data.map,D(idx_dataset).name]);
     
    % Arrange and set index variables
    % Glue conditons and location to data    
    % @return IPOL.HS and IPOL.ANT as HotSpot and Anterior Location in xyz  
    % @return IPOL.COND as Condition string
    % @return IPOL.DESIGN as Design Matrix (Biphasic / Lateromedial) yes/no
    t.s                     = regexp(D(idx_dataset).name,'S\w*');
    t.c                     = regexp(D(idx_dataset).name,'C\w*');
    t.d                     = regexp(D(idx_dataset).name,'.mat');
    i.data_sub              = int32(str2double(D(idx_dataset).name(t.s+1:t.c-1)));
    i.data_cond             = int32(str2double(D(idx_dataset).name(t.c+1:t.d-1)));
    i.sub                   = find(setup.SUB.id==i.data_sub);
    t.sex                   = setup.SUB.sex(i.sub);
    t.age                   = setup.SUB.age(i.sub);
    t.idx                   = i.data_sub;
    t.pos                   = setup.SUB.pos{i.sub};
    t.HS                    = [t.pos{i.data_cond+1,2:4}];
    t.ANT                   = [t.pos{i.data_cond+1,5:7}];
    
    IPOL(idx_dataset).HS        = t.HS;
    IPOL(idx_dataset).ANT       = t.ANT;
    IPOL(idx_dataset).COND      = setup.MAP.label.all{i.data_cond};
    IPOL(idx_dataset).DESIGN    = logical([setup.MAP.BI(i.data_cond),setup.MAP.LM(i.data_cond)]);
    IPOL(idx_dataset).sub_sex   = t.sex;
    IPOL(idx_dataset).sub_age   = t.age;
    IPOL(idx_dataset).sub_idx   = t.idx;
    
    clear t i

    % Interpolate individual subjects amplitude and MEP+ on the group grid based on squared distance
    % Interpolate individual subjects latency after centering on average
    quad_weights                =  utils.distanceweights_groupgrid(mapping.xyz,grid_xyz,25);

    sub_mep                     = double(mapping.amp>50 & mapping.lat>15);    
    sub_lat                     = mapping.lat;
    sub_lat(sub_lat>0)          = sub_lat(sub_lat>0)-mean(sub_lat(sub_lat>0));
    
    IPOL(idx_dataset).quad_AMP  = reshape_to_grid((quad_weights*mapping.amp));       
    IPOL(idx_dataset).quad_MEP  = reshape_to_grid((quad_weights*sub_mep));     
    IPOL(idx_dataset).quad_LAT  = reshape_to_grid((quad_weights*sub_lat));     
    
    fprintf(logfileid,'At %s finished dataset %i \n',datetime('now'),idx_dataset);
end
fclose(logfileid);
save([folder.results.stats,'map_interpolated.mat'],'IPOL')
%%
DESIGN  = cat(1,IPOL.DESIGN);
SUB     = cat(1,IPOL.sub_idx);
AMP     = cat(4,IPOL.quad_AMP);
LAT     = cat(4,IPOL.quad_LAT);
MEP     = cat(4,IPOL.quad_MEP);

%


















%%
F = [];
for i_x=1:map_dimensions(1),
    for i_y=1:map_dimensions(2),
        %for i_z=1:map_dimensions(3),
        i_z = 1;
            [p,tab,stats]       = anovan(squeeze(nanmean(MEP(i_x,i_y,:,:),3)),cat(2,DESIGN,SUB),'display','off','model',[1 0 0;0 1 0;1 1 0;0 0 1]);
            F(i_x,i_y,i_z,:)    = [tab{2:4,6}];
        %end
    end
end



close all
for k=1:3,
    figure
    set(gcf,'position',[100 100 [10,10].*map_dimensions(1:2)])
    contour(X(:,:,1),Y(:,:,1),F(:,:,k),25)
    hold on
    colorbar
    plot(-31,-31,'ko','markerfacecolor','k')
end

%%    Code Test by Visual Inspection / Plotting    
x = nanmean(LAT,4);
[h,p,ci,stats] = ttest(LAT,0,'dim',4);
[h,p,ci,stats] = ttest(AMP,0,'dim',4);
x = median(MEP,4);
%x = stats.tstat;
%x = p<0.05;

hold on
contour(X(:,:,1),Y(:,:,1),squeeze(nanmean(x,3)),25)
view([180,-90])

%colormap(diverging_bwr_40_95_c42_n256)