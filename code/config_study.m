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
folder.data.ioc     = 'C:\PROJECTS\Subject Studies\TMS-MAP-IOC\zwischen\ioc\';
folder.data.map     = 'C:\PROJECTS\Subject Studies\TMS-MAP-IOC\zwischen\mapping\';
folder.results      = 'C:\PROJECTS\Subject Studies\TMS-MAP-IOC\results\';
folder.code         = 'C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\';
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
% based on tabulated data
[num,str,all]   = xlsread('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\zwischen\subject_data.xlsx');
setup.SUB.id    = num(:,1); %corresponding to S% in datasets
setup.SUB.sex   = logical(num(:,1)); %male 1; female 0;
setup.SUB.age   = num(:,2); %in years
%% Preloading the headmodel for mapping visualization
load('C:\Users\Robert Bauer\Documents\MATLAB\private_toolbox\symetric_headmodel');
%%
save([folder.code,'config.mat'],'headmodel','setup','folder')