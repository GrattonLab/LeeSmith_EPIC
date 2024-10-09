%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
% Rank order consistency across tasks
%
% To run this script:
% You need mat-files named "FL_CEmat.mat","PP_CEmat.mat", and "ST_CEmat.mat"
% (You'll have these files once you run EPIC_preprocess_flanker.m, 
% EPIC_preprocess_primeprobe.m, and EPIC_preprocess_Stroop.m)
%
% What this script does:
% This script examines the rank order consistency of the congruency effect
% across the flanker, prime-probe, and Stroop tasks.
%
% What this script outputs:
% Matrices of rank orders for the three tasks regarding the congruency
% effect
% Supp. Fig. 5
%
% Created on 12/19/2023 by HJ Lee
% Last modified on 12/19/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')
%% Parameter settings
tmpSubjID = 1:12;
subjExcld = [1,2,11,9];  % excluded participants
subjID = setdiff(tmpSubjID, subjExcld);
nSubj = length(subjID);

nTask = 3;
palette = parula(nTask);
palette(3,:) = [1 0 0];  % red

%% Load processed data
load FL_CEmat.mat
load PP_CEmat.mat
load ST_CEmat.mat

%% 1. Rank order across tasks
% Preassigment
rnkRTmat = nan(nTask,nSubj,3);  % 3: congruency effect, congruent, incongruent
rnkACCmat = nan(nTask,nSubj,3);
%% Reaction time
% Flanker task
% Congruency effect
[~,p1] = sort(FL_CErtGrandMean,'descend');
r1 = 1:nSubj;
r1(p1) = r1;
rnkRTmat(1,:,1) = r1;

% Prime-probe task
% Congruency effect
[~,p2] = sort(PP_CErtGrandMean,'descend');
r2 = 1:nSubj;
r2(p2) = r2;
rnkRTmat(2,:,1) = r2;

% Stroop task
% Congruency effect
[~,p3] = sort(ST_CErtGrandMean,'descend');
r3 = 1:nSubj;
r3(p3) = r3;
rnkRTmat(3,:,1) = r3;

%% Accuracy
% Flanker task
% Congruency effect
[~,p1] = sort(FL_CEaccGrandMean,'descend');
r1 = 1:nSubj;
r1(p1) = r1;
rnkACCmat(1,:,1) = r1;

% Prime-probe task
% Congruency effect
[~,p2] = sort(PP_CEaccGrandMean,'descend');
r2 = 1:nSubj;
r2(p2) = r2;
rnkACCmat(2,:,1) = r2;

% Stroop task
% Congruency effect
[~,p3] = sort(ST_CEaccGrandMean,'descend');
r3 = 1:nSubj;
r3(p3) = r3;
rnkACCmat(3,:,1) = r3;

%% Plot matrices - Heatmap
% Reaction time
% Congruency effect
figure
xvalues = subjID;
yvalues = {'Flanker','Prime-Probe','Stroop'};
h1 = heatmap(xvalues,yvalues,rnkRTmat(:,:,1));
set(gca,'FontSize',18)
h1.XLabel = 'Participant ID';
%h1.YLabel = 'Task';
h1.FontSize = 20;
h1.Colormap = parula;
h1.Title = 'Congruency Effect Reaction Time';

% Accuracy
% Congruency effect
figure
h2 = heatmap(xvalues,yvalues,rnkACCmat(:,:,1));
set(gca,'FontSize',18)
h2.XLabel = 'Participant ID';
%h4.YLabel = 'Task';
h2.FontSize = 18;
h2.Colormap = parula;
h2.Title = 'Congruency Effect Accuracy';