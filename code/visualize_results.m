%% Configuration
addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
utils.cls; %clc, clear, close all, matlabrc, fclose('all'), addpath of Fieldtrip and other toolboxes, ft_defaults
addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code'); %to access the package folder +utils
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\config.mat','gridmodel','headmodel','setup','folder','Label');
addpath('C:\Users\Robert Bauer\Documents\MATLAB\other_toolboxes\CETperceptual_MATLAB'); %folder with colormaps
cmap = diverging_bwr_40_95_c42_n256;
set(groot,'DefaultFigureColormap',cmap)
%% SUB parameters
load('map_subject_data.mat')
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

%% Plot Cluster Results of the cubic interpolation datasets for all three measures
% of interest (i.e. ('Amplitude' ; 'MEP' ; 'Latency')
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

%% Plot Position of Anterior Target
close all
load([folder.results.stats,'map_subject_data.mat'],'SUB')
ANT         = utils.get_GroupAnt(SUB);
V           = utils.Target2Grid(ANT);
utils.plot_headmodel(headmodel,V,1)
annotation('textbox','Position',[0 0 1 1],'String','Distribution of Anterior Target')
print(gcf,[folder.results.figs,'\HEAD_PDF_Anterior_Target.tif'],'-dtiff')

utils.plot_gridmodel(V,2)
title('Distribution of Anterior Target')
print(gcf,[folder.results.figs,'\GRID_PDF_Anterior_Target.tif'],'-dtiff')

%% 

% Visualize IOC for parameters (i.e. 'Amplitude' ; 'MEP' ; 'Latency')
labels_list = cat(1,setup.IO.label.BI,setup.IO.label.LM,setup.IO.label.M1);
DESIGN      = logical(cat(1,setup.IO.BI,setup.IO.LM,setup.IO.M1));
DESIGN_labels = cat(1,setup.IO.label.BI(~setup.IO.BI+1),setup.IO.label.LM(~setup.IO.LM+1),setup.IO.label.M1(~setup.IO.M1+1));
for i_d = 1:3
  
    loadfile            = [folder.results.stats,'ioc\',Label.ioc_Field{i_d},'_stats.mat'];   
    load(loadfile,'TestResults','ClusterResults');     
    
    close all
    for design_idx =1:3,

        fig_DESIGN  = cat(1,DESIGN(design_idx,:),~DESIGN(design_idx,:));
        fig_LABEL   = [[unique(DESIGN_labels(design_idx,DESIGN(design_idx,:),:)),unique(DESIGN_labels(design_idx,~DESIGN(design_idx,:),:))],Label.ioc_Field{i_d}];
        utils.plot_ioc(TestResults,fig_DESIGN,fig_LABEL)
 
    end    
    
   % close all
%     for design_idx =1:2,
%         a           = DESIGN(design_idx,:) & setup.IO.M1;
%         b           = ~DESIGN(design_idx,:) & setup.IO.M1;
%         fig_DESIGN  = cat(1,a,b);
%         trgt        = unique(cat(1,DESIGN_labels(3,a),DESIGN_labels(3,b)));
%         fig_LABEL   = [unique(DESIGN_labels(design_idx,a)),unique(DESIGN_labels(design_idx,b)),{[trgt{1},' ',Label.ioc_Field{i_d}]}];
%         utils.plot_ioc(TestResults,fig_DESIGN,fig_LABEL)
%         
%         a           = DESIGN(design_idx,:) & ~setup.IO.M1;
%         b           = ~DESIGN(design_idx,:) & ~setup.IO.M1;
%         fig_DESIGN  = cat(1,a,b);
%         trgt        = unique(cat(1,DESIGN_labels(3,a),DESIGN_labels(3,b)));
%         fig_LABEL   = [unique(DESIGN_labels(design_idx,a)),unique(DESIGN_labels(design_idx,b)),{[trgt{1},' ',Label.ioc_Field{i_d}]}];
%         utils.plot_ioc(TestResults,fig_DESIGN,fig_LABEL)
% 
%     end   
    
%     % Interaction
%     a           = ((setup.IO.BI & setup.IO.LM) | (~setup.IO.BI & ~setup.IO.LM)) & setup.IO.M1;
%     b           = ((setup.IO.BI & ~setup.IO.LM) | (~setup.IO.BI & setup.IO.LM)) & setup.IO.M1;
%     fig_DESIGN  = cat(1,a,b);
%     trgt        = unique(cat(1,DESIGN_labels(3,a),DESIGN_labels(3,b)));
%     fig_LABEL   = [{'Bi 90° | Mono 45°'},{'Bi 45° | Mono 90°'},{[trgt{1},' ',Label.ioc_Field{i_d}]}];
%     utils.plot_ioc(TestResults,fig_DESIGN,fig_LABEL)
%     
%     a           = ((setup.IO.BI & setup.IO.LM) | (~setup.IO.BI & ~setup.IO.LM)) & ~setup.IO.M1;
%     b           = ((setup.IO.BI & ~setup.IO.LM) | (~setup.IO.BI & setup.IO.LM)) & ~setup.IO.M1;
%     fig_DESIGN  = cat(1,a,b);
%     trgt        = unique(cat(1,DESIGN_labels(3,a),DESIGN_labels(3,b)));
%     fig_LABEL   = [{'Bi 90° | Mono 45°'},{'Bi 45° | Mono 90°'},{[trgt{1},' ',Label.ioc_Field{i_d}]}];
%     utils.plot_ioc(TestResults,fig_DESIGN,fig_LABEL)
%     
    % Full Model
    a           = ((setup.IO.BI & setup.IO.LM) | (~setup.IO.BI & ~setup.IO.LM)) & setup.IO.M1;
    b           = ((setup.IO.BI & setup.IO.LM) | (~setup.IO.BI & ~setup.IO.LM)) & ~setup.IO.M1;
    fig_DESIGN  = cat(1,a,b);
    fig_LABEL   = [{'Bi 90° | Mono 45° on M1'},{'Bi 90° | Mono 45° on NPMA'},Label.ioc_Field{i_d}];
    utils.plot_ioc(TestResults,fig_DESIGN,fig_LABEL)
        
    a           = ((setup.IO.BI & ~setup.IO.LM) | (~setup.IO.BI & setup.IO.LM)) & setup.IO.M1;    
    b           = ((setup.IO.BI & ~setup.IO.LM) | (~setup.IO.BI & setup.IO.LM)) & ~setup.IO.M1;
    fig_DESIGN  = cat(1,a,b);
    fig_LABEL   = [{'Bi 45° | Mono 90° on M1'},{'Bi 45° | Mono 90° on NPMA'},Label.ioc_Field{i_d}];
    utils.plot_ioc(TestResults,fig_DESIGN,fig_LABEL)

    
end
%%
    fig_DESIGN  = cat(1,ismember(setup.IO.label.all,{'Bi l-m M1','Mo pl-am M1'}),ismember(setup.IO.label.all,{'Bi pl-am M1','Mo l-m M1'}));
    fig_LABEL   = [{['Bi 90°',' & ','Mo 45 °']},{['Bi 45°',' & ','Mo 90°']},{[Label.ioc_Field{i_d}, ' M1']}];
    utils.plot_ioc(TestResults,fig_DESIGN,fig_LABEL)
    
    fig_DESIGN  = cat(1,(setup.IO.BI & setup.IO.M1),(setup.IO.BI & ~setup.IO.M1));
    fig_LABEL   = [{'M1'},{'NPMA'},{['Biphasic ',Label.ioc_Field{i_d}]}];
    utils.plot_ioc(TestResults,fig_DESIGN,fig_LABEL)
    
    fig_DESIGN  = cat(1,(~setup.IO.BI & setup.IO.M1),(~setup.IO.BI & ~setup.IO.M1));
    fig_LABEL   = [{'M1'},{'NPMA'},{['Monophasic ',Label.ioc_Field{i_d}]}];
    utils.plot_ioc(TestResults,fig_DESIGN,fig_LABEL)

     