addpath('C:\Users\Robert Bauer\Documents\Matlab\private_toolbox');
cls; %clc, clear all, close all, matlabrc, addpath of Fieldtrip and other toolboxes, ft_defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The dataset has been recorded by Cornelius Kuhn (CK) as part of his medical
% thesis work in the framework of an IZKF-stipend funding at the
% Translational Neurosurgery Research Group of the University Hopsital
% Tübingen. Supervision was performed by Robert Bauer (RB) and Alireza Gharabaghi (AG).
% Data was already visually inspected and organization streamlined by CK using a GUI developed by RB. 
% This scripts are written by RB with the purpose to analyze the input-output 
% curves and maps based on the visually inspected data; 
% RB follows a three-script approach. One for predefinitions and global
% settings (config.m, producing config.mat), one for signal processing
% (process.m) and one for statistical analysis and visualization
% (results.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining the folder with the visually inspected electrophysiological Data, and predefine processing output folder
folder.data.ioc         = 'C:\PROJECTS\Subject Studies\TMS-MAP-IOC\zwischen\ioc\';
folder.data.map         = 'C:\PROJECTS\Subject Studies\TMS-MAP-IOC\zwischen\mapping\';
folder.data.results     = 'C:\PROJECTS\Subject Studies\TMS-MAP-IOC\zwischen\results\';
folder.results.stats    = 'C:\PROJECTS\Subject Studies\TMS-MAP-IOC\results\stats\';
folder.results.figs     = 'C:\PROJECTS\Subject Studies\TMS-MAP-IOC\results\figs\';
folder.code             = 'C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\';
%% Data Organization and Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IOC Data is stored as .mat files folder.data.ioc according to the following schema
% S%C%.mat, with the % following S indicating anonymized subject ID; and
% the % following C indicating stimulation condition.
% Each number has a specific label. Additionally consider the condition
% mixing matrix defined by BI/LM/M1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setup.IO.label.all  = {'Bi l-m A','Bi l-m M1','Bi pl-am A','Bi pl-am M1','Mo l-m A','Mo l-m M1','Mo pl-am A','Mo pl-am M1'};
setup.IO.BI         = [1 1 1 1 0 0 0 0]; % 1: Biphasic; 0: Monophasic
setup.IO.LM         = [1 1 0 0 1 1 0 0]; % 1: Lateral to Medial Current flow (90°); 0: Posterior-lateral to anterior-medial current flow (45°)
setup.IO.M1         = [0 1 0 1 0 1 0 1]; % 1: Stimulation over M1 hotspot; 0: Stimulation over anterior hotspot (i.e. non-primary motor areas (NPMA))
setup.IO.SI         = [90,100,110,120,130,140,150]; % in percent of Resting motor threshold 
setup.IO.label.BI   = {'Biphasic','Monophasic'};
setup.IO.label.LM   = {'90°','45°'};
setup.IO.label.M1   = {'M1','NPMA'};
setup.IO.label.ORG  = {'Waveform','Orientation','Location'};
setup.IO.label.VOI  = {'Amplitude','MEP','Latency'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mapping Data is stored as .mat files in the folder.data.map according to the following schema
% S%C%.mat, with the % following S indicating anonymized subject ID; and
% the % following C indicating stimulation condition.
% Each number has a specific label. Additionally consider the condition
% mixing matrix defined by BI/LM/M1. Consider that mapping was performed at
% fixed intensity. Three exclamationmarks after the filename indicate that
% this condition was not recorded and/or had to be removed due to arficat
% contamination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setup.MAP.label.all = {'Bi_l-m'  'Bi_pl-am'  'Mono_l-m'  'Mono_pl-am'};
setup.MAP.BI        = [1 1 0 0]; % 1: Biphasic; 0: Monophasic
setup.MAP.LM        = [1 0 1 0]; % 1: Lateral to Medial Current flow (90°); 0: Posterior-lateral to anterior-medial current flow (45°)
setup.MAP.label.BI  = {'Biphasic','Monophasic'};
setup.MAP.label.LM  = {'90°','45°'};
%% Preloading subject parameters (age,sex,rmt)
% based on tabulated data collected by CK
load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\zwischen\IOC_positions.mat');
[num,str,all]   = xlsread('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\zwischen\subject_data.xlsx');
setup.SUB.id    = num(:,1); %corresponding to S% in datasets
setup.SUB.sex   = logical(num(:,2)); %male 1; female 0;
setup.SUB.age   = num(:,3); %in years

for idx_put=1:length(setup.SUB.id),
    idx_take                = setup.SUB.id(idx_put);    
    setup.SUB.pos{idx_put}  = COORD(idx_take).positions;
    setup.HS(idx_put,:)     = nanmean(cat(1,[COORD(idx_take).positions{2:end,2}],[COORD(idx_take).positions{2:end,3}],[COORD(idx_take).positions{2:end,4}]),2);
    setup.ANT(idx_put,:)    = nanmean(cat(1,[COORD(idx_take).positions{2:end,5}],[COORD(idx_take).positions{2:end,6}],[COORD(idx_take).positions{2:end,7}]),2);    
end
%%
% Define the query grid based on surface mesh-grid 
%sfc                 = load('C:\Users\Robert Bauer\Documents\Matlab\other_toolboxes\fieldtrip\template\anatomy\surface_white_left.mat');
sfc                 = load('C:\Users\Robert Bauer\Documents\Matlab\other_toolboxes\fieldtrip\template\anatomy\surface_pial_left.mat');
clear headmodel
headmodel.pos         = single(sfc.bnd.pnt);
headmodel.tri         = single(sfc.bnd.tri);
clear sfc
% Aligned Surface 
%gridmodel.resolution                    = [25,50,3];
%[gridmodel.X,gridmodel.Y,gridmodel.Z]   = (meshgrid(linspace(-35,35,gridmodel.resolution(1)),linspace(-65,65,gridmodel.resolution(2)),linspace(-1,1,gridmodel.resolution(3))));
%gridmodel.pos                           = single(cat(2,reshape(gridmodel.X,[],1),reshape(gridmodel.Y,[],1),reshape(gridmodel.Z,[],1)));
%
gridmodel.resolution                    = [50,100,1];
[gridmodel.X,gridmodel.Y,gridmodel.Z]   = (meshgrid(linspace(-25,25,gridmodel.resolution(1)),linspace(-50,50,gridmodel.resolution(2)),linspace(0,0,gridmodel.resolution(3))));
gridmodel.pos                           = single(cat(2,reshape(gridmodel.X,[],1),reshape(gridmodel.Y,[],1),reshape(gridmodel.Z,[],1)));
%% Define Labels and Titles for Plotting
Label.Title     = {[setup.MAP.label.BI{1},' > ',setup.MAP.label.BI{2}]...
            [setup.MAP.label.LM{1},' > ',setup.MAP.label.LM{2}],...
            [setup.MAP.label.BI{1},' ',setup.MAP.label.LM{1},' & ',setup.MAP.label.BI{2},' ',setup.MAP.label.LM{2}],...
            };
       
Label.Save      =  {'Waveform', 'Orientation', 'W-O-Interaction'};
Label.Field     = {'Amplitude','MEP','Latency'};
Label.Weight    = {'Quadratic','Spherical'};   
Label.Dataset   = {'quad_AMP'    'quad_MEP'    'quad_LAT'    'th_AMP'    'th_MEP'    'th_LAT'};
Label.ioc_Field = {'Amplitude','MEP','Latency','TimeCourse'};
%%
save([folder.code,'config.mat'],'headmodel','gridmodel','setup','folder','Label')