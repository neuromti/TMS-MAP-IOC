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

    
    % Remove trials with artifacts, calculate MEP+ and center latency,
    % calculate CoG
    sub                         = utils.prepare_for_weighting(mapping);
    
    % Interpolate individual subjects MEP+ and latency unto the surface based on squared distance    
    quad_weights                = utils.calculate_distanceweights(sub.xyz,queried.pos,25);
       
    IPOL(idx_dataset).quad_AMP  = quad_weights*sub.amp;
    IPOL(idx_dataset).quad_MEP  = quad_weights*sub.mep;     
    IPOL(idx_dataset).quad_LAT  = quad_weights*sub.lat;     

    % Interpolate as above, but based on CoG shifted to M1      ´    
    sub                         = utils.shift_by_CoG(sub);
    shifted_quad_weights        = utils.calculate_distanceweights(sub.shifted_xyz,queried.pos,25);    

    IPOL(idx_dataset).CoG               = sub.CoG;
    IPOL(idx_dataset).shifted_quad_AMP  = shifted_quad_weights*sub.amp;
    IPOL(idx_dataset).shifted_quad_MEP  = shifted_quad_weights*sub.mep;     
    IPOL(idx_dataset).shifted_quad_LAT  = shifted_quad_weights*sub.lat;     
    
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
DESIGN  = cat(1,IPOL.DESIGN);
AMP     = cat(2,IPOL.quad_AMP);
LAT     = cat(2,IPOL.quad_LAT);
MEP     = cat(2,IPOL.quad_MEP);

SUB     = SUB(sort_idx);
DESIGN  = DESIGN(sort_idx,:);
AMP     = AMP(:,:,:,sort_idx);
MEP     = MEP(:,:,:,sort_idx);
LAT     = LAT(:,:,:,sort_idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORM STATISTICAL ANALYSIS FOR EACH VARIABLE OF INTEREST
num_rep     = 1000;

utils.do_botandperm_MAP(AMP,SUB,DESIGN,folder.results.stats,num_rep)






%%

x = reshape_to_grid(mean(true_MM,2));
x = reshape_to_grid(mean(true_MM(:,1)-true_MM(:,2),2));
close all
figure
for z=4:4:64
    subplot(4,4,(z/4))
    %contour(X(:,:,1),Y(:,:,1),squeeze(nanmean(x(:,:,z),3)),[0:1:125],'fill','on');   
    %caxis([0 125])

    contour(X(:,:,1),Y(:,:,1),squeeze(nanmean(x(:,:,z),3)),[-50:50],'fill','on');   
    caxis([-50 50])
    hold on
    plot(-37,-21,'ko','markerfacecolor','k')
    view([180,-90])
end


%%

x = reshape_to_grid(true_STATVAL(:,2));
x = -log10(reshape_to_grid(1-fcdf(true_STATVAL(:,2),1,31)));
%x = -log10(reshape_to_grid(1-fcdf(true_STATVAL(:,1),1,31)));
close all
figure
%contour(X(:,:,1),Y(:,:,1),squeeze(nanmean(x,3)),25)
[~,hdl_c] = contour(X(:,:,1),Y(:,:,1),squeeze(nanmean(x,3)),1.3:.05:3,'fill','on');
%
caxis(gca,[1.3 3])
hold on
plot(-37,-21,'ko','markerfacecolor','k')
view([180,-90])



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