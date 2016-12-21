%% Configuration
matlabrc;
status = fclose('all');
close all;
clear;
clc;

addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %study package +utils
addpath('C:\Users\Robert Bauer\Documents\MATLAB\other_toolboxes\CETperceptual_MATLAB'); %folder with colormaps
addpath ('C:\Users\Robert Bauer\Documents\Matlab\other_toolboxes\fieldtrip');
ft_defaults;

load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','gridmodel','headmodel','setup','folder','Label');
set(groot,'DefaultFigureColormap',diverging_bwr_40_95_c42_n256())
%% SUB parameters
load([folder.results.stats,'map_subject_data.mat'])
SUBID = cat(1,SUB.subID);
sum(grpstats(cat(1,SUB.sex),SUBID))
mean(grpstats(cat(1,SUB.age),SUBID))
std(grpstats(cat(1,SUB.age),SUBID))
min(grpstats(cat(1,SUB.age),SUBID))
max(grpstats(cat(1,SUB.age),SUBID))
%% M1 Estimation
a               = (grpstats(cat(1,SUB.GridOrigin),SUBID))+[0,-10,0]; %grid origin is 10mm anterior to HotSpot
[h,p,ci,stats]  = ttest(a-utils.get_M1())
%% Get Mean and SD of RMT used for Mapping
[num,str,all]   = xlsread('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\zwischen\subject_data.xlsx',2);
% remember: mapping was performed at 110% of RMT
clc
close force all
DATA = num(:,1);
DSGN = cat(2,repmat([1,1,0,0],1,13)',repmat([1,0,1,0],1,13)',reshape(repmat(1:13,4,1),1,[])'); %BI vs MO, LM vs PLAM, SUB
[p,tab,stats]   = anovan(DATA,DSGN,'random',3,'model',[1 0 0;0 1 0;1 1 0;0 0 1],'display','off','Varnames',{'Biphasic','90°','Subject'});
cat(2,unique(all(:,4)),num2cell(grpstats(num(:,1),all(:,4),'mean')),num2cell(grpstats(num(:,1),all(:,4),'std')))
%% VISUALIZATION OF MAPPING RESULTS
% Show Influence of Factors on the different measure
% Visualize them on the grid and project them on a headmodel

% Delete already printed figures to prevent confusion about the correct version 
delete([folder.results.figs,Label.Weight{1},'\*.tif'])
delete([folder.results.figs,Label.Weight{2},'\*.tif'])

% Visualization of the cubic interpolation datasets for all three measures
% of interest (i.e. ('Amplitude' ; 'MEP' ; 'Latency')
for i_d = 1:length(Label.Dataset(1:3))
    load([folder.results.stats,Label.Dataset{i_d},'\stats.mat']);            
    i_field     = find(ismember({'AMP','MEP','LAT'},Label.Dataset{i_d}(end-2:end)));
    i_weight    = find(ismember({'qu','th'},Label.Dataset{i_d}(1:2)));
  
    close all    
    % for all modeled factors of interest (i.e. 'Biphasic > Monophasic' ; '90° > 45°' ; 'Biphasic 90° & Monophasic 45°' ; but not 'Subject')
    for k=1:size(TestResults.Pval,2)        

        % Prepare data for plotting
        PlotLabel   = [Label.Title{k},' ->',Label.Field{i_field}];         
        P           = -log10(TestResults.PermutationPval(:,k));
        S           = sign(TestResults.Coeffs(:,k));
        V           = P.*S;
        
        % Plot on Grid using the Design Grid
        utils.plot_gridmodel(V,2)
        title(PlotLabel)             
        print(gcf,[folder.results.figs,Label.Weight{i_weight},'\GRID_',Label.Save{k},'-',Label.Field{i_field},'.tif'],'-dtiff')        
        
        % Plot on Cortical Surface / Headmodel
        utils.plot_headmodel(headmodel,V,2)
        annotation('textbox','Position',[0 0 1 1],'String',PlotLabel)
        print(gcf,[folder.results.figs,Label.Weight{i_weight},'\HEAD_',Label.Save{k},'-',Label.Field{i_field},'.tif'],'-dtiff')

   
    end
       
end
close all
%% Plot Cluster Results of the cubic interpolation datasets for all three measures
% of interest (i.e. ('Amplitude' ; 'MEP' ; 'Latency')

COG             = {'Factor','Cluster ID','x','y','Area','Spread','PermPVal'}
for i_d = 1:length(Label.Dataset(1:3))
    load([folder.results.stats,Label.Dataset{i_d},'\stats.mat']);            
    i_field     = find(ismember({'AMP','MEP','LAT'},Label.Dataset{i_d}(end-2:end)));
    i_weight    = find(ismember({'qu','th'},Label.Dataset{i_d}(1:2)));
    close all     
    % for all modeled factors of interest (i.e. 'Biphasic > Monophasic' ; '90° > 45°' ; 'Biphasic 90° & Monophasic 45°' ; but not 'Subject')
    for k=1:size(TestResults.Pval,2)   
        PlotLabel   = [Label.Title{k},' ->',Label.Field{i_field}];         
        if ~isempty(ClusterResults(k).PermPval)
            for clus_idx = 1:length(ClusterResults(k).PermPval),
                
                % Prepare data for plotting
                ClusterLabel    = [PlotLabel,' ',sprintf('(Cluster %i: p=%.2g )',clus_idx,ClusterResults(k).PermPval(clus_idx))];
                ClusterMember   = ClusterResults(k).MemberShip(:,:,1)==clus_idx;
                [CoG,Spread,Area]      = utils.get_TopographicalParameters(ClusterMember);
                COG             = cat(1,COG,{PlotLabel,['Cluster #',num2str(k)],CoG(1),CoG(2),Area,Spread,ClusterResults(k).PermPval(clus_idx)});
                
                V               = utils.Clustermember2grid(ClusterMember);
                
                % Plot on Grid using the Design Grid
                utils.plot_gridmodel(V,2)
                title(ClusterLabel)   
                colorbar off
                print(gcf,[folder.results.figs,Label.Weight{i_weight},'\CLUSTER_GRID_',Label.Save{k},'-',Label.Field{i_field},'-',num2str(clus_idx),'.tif'],'-dtiff')
            
                % Plot on Cortical Surface / Headmodel
                utils.plot_headmodel(headmodel,V,2)
                annotation('textbox','Position',[0 0 1 1],'String',ClusterLabel)
                print(gcf,[folder.results.figs,Label.Weight{i_weight},'\CLUSTER_HEAD_',Label.Save{k},'-',Label.Field{i_field},'-',num2str(clus_idx),'.tif'],'-dtiff')

            end                
        end
    end
end
if exist([folder.results.stats,'cluster_topo.xlsx'],'file'), delete([folder.results.stats,'cluster_topo.xlsx']); end
xlswrite([folder.results.stats,'cluster_topo.xlsx'],COG);
close all

%% Plot Position of Anterior Target
close all
load([folder.results.stats,'map_subject_data.mat'],'SUB')
[xyz,names] = utils.get_landmark();
ANT         = utils.get_GroupAnt(SUB);
% closest region is PMd, but still significant anterior to it. 
P = []; CI = [];
for k = 1 : length(names)   
    [h,p,ci,stats] = ttest(ANT-xyz(k,:)); 
    P = cat(1,P,p); 
    CI = cat(3,CI,ci); 
end
table(names,P,squeeze(mean(CI,1))',squeeze(CI(1,:,:))',squeeze(CI(2,:,:))')


LandMarkDistance = pdist2(xyz(:,1:2),ANT(:,1:2));
figure
hold on
barh(mean(pdist2(xyz(:,1:2),ANT(:,1:2)),2),'w')
for lm_i = 1 : length(names)
    scatter(LandMarkDistance(lm_i,:),normrnd(lm_i,0.075,13,1),'ko','markerfacecolor','k')
end
set(gca,'YTICKLABEL',names,'YTICK',1 : length(names),'YLIM',[0 length(names)+1])
xlabel('Distance in mm')
colormap([.5 0 0])
print(gcf,[folder.results.figs,'\Landmark_to_Anterior_Target.tif'],'-dtiff')
%
ANT         = utils.get_GroupAnt(SUB);
V           = utils.Target2Grid(ANT);
utils.plot_headmodel(headmodel,V,1)
annotation('textbox','Position',[0 0 1 1],'String','Distribution of Anterior Target')
print(gcf,[folder.results.figs,'\HEAD_PDF_Anterior_Target.tif'],'-dtiff')

utils.plot_gridmodel(V,2)
title('Distribution of Anterior Target')
hnd_cb                      = colorbar;
set(hnd_cb,'YLIM',[0 1],'YTICK',[0:0.2:1],'YTICKLABEL',round(linspace(0,max(V)./sum(V),6)*1000)./1000)
print(gcf,[folder.results.figs,'\GRID_PDF_Anterior_Target.tif'],'-dtiff')

utils.plot_landmark2grid();
print(gcf,[folder.results.figs,'\GRID_PDF_Anterior_Target_w_landmarks.tif'],'-dtiff')
%% Visualize IOC for parameters (i.e. 'Amplitude' ; 'MEP' ; 'Latency')
labels_list = cat(1,setup.IO.label.BI,setup.IO.label.LM,setup.IO.label.M1);
DESIGN      = logical(cat(1,setup.IO.BI,setup.IO.LM,setup.IO.M1));
DESIGN_labels = cat(1,setup.IO.label.BI(~setup.IO.BI+1),setup.IO.label.LM(~setup.IO.LM+1),setup.IO.label.M1(~setup.IO.M1+1));
if exist([folder.results.figs,'IOC\']),delete([folder.results.figs,'IOC\*.*']), else mkdir([folder.results.figs,'IOC\']); end
for i_d = [1,3]
  
    loadfile            = [folder.results.stats,'ioc\',Label.ioc_Field{i_d},'_stats.mat'];   
    load(loadfile,'TestResults','ClusterResults');     
    
    close all
    for design_idx =1:3
        fig_DESIGN  = cat(1,DESIGN(design_idx,:),~DESIGN(design_idx,:));
        fig_LABEL   = [[unique(DESIGN_labels(design_idx,DESIGN(design_idx,:),:)),unique(DESIGN_labels(design_idx,~DESIGN(design_idx,:),:))],Label.ioc_Field{i_d}];
        Pval        = TestResults.PermutationPval(:,design_idx);
        utils.plot_ioc(TestResults.MargMeans,TestResults.CI,fig_DESIGN,fig_LABEL,Pval);
        print(gcf,[folder.results.figs,'IOC\',sprintf('%s_%s_%s',fig_LABEL{:}),'.tif'],'-dtiff')
    end

end
%% TIMECOURSE
loadfile            = [folder.results.stats,'ioc\',Label.ioc_Field{4},'_stats.mat'];   
load(loadfile,'TestResults');     
if exist([folder.results.figs,'TC\']),delete([folder.results.figs,'TC\*.*']), else mkdir([folder.results.figs,'TC\']); end

MargMeans   = (cat(3,TestResults.MargMeans));
Pval        = (cat(3,TestResults.PermutationPval));

for design_idx =1:3
    fig_DESIGN  = cat(1,DESIGN(design_idx,:),~DESIGN(design_idx,:));
    fig_LABEL   = [[unique(DESIGN_labels(design_idx,DESIGN(design_idx,:),:)),unique(DESIGN_labels(design_idx,~DESIGN(design_idx,:),:))],Label.ioc_Field{i_d}];
    
    close all
  	[h] = utils.plot_timecourse(MargMeans,fig_DESIGN,fig_LABEL,squeeze(Pval(:,design_idx,:)));
    for f=1:length(h)
        print(h(f),[folder.results.figs,'TC\',sprintf('%s_%s_%s_%i',fig_LABEL{1:2},'TC',f),'.tif'],'-dtiff')
    end
    
end





     