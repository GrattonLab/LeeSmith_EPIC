%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
% Flanker, prime-probe, and Stroop tasks
%
% To run this script:
% You need mat-files named "FL_CEmat.mat","PP_CEmat.mat", and "ST_CEmat.mat"
% (You'll have these files once you run EPIC_preprocess_flanker.m, 
% EPIC_preprocess_primeprobe.m, and EPIC_preprocess_Stroop.m)
% You also need violinplot.m (download at
% https://github.com/bastibe/Violinplot-Matlab/blob/master/violinplot.m)
%
% What this script does:
% Plots the violin plots of 18 sessions and grand mean plots of all session
%
% What this script outputs:
% Fig. 2: congruency effect (CE) RT
% Supp. Fig. 1: CE accuracy
% Supp. Fig. 2: CE IES
% Supp. Fig. 3: Congruent/Incongruent RT
% Supp. Fig. 4: Congruent/Incongruent accuracy
% Supp. Fig. 5: Congruent/Incongruent IES
%
% Created on 12/10/2022 by HJ Lee
% Last modified on 06/29/2023
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

%% Violin plot and Grand mean plot together
% (1) CE
% RT
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CErtMat,[],'ViolinColor',palette(1,:))
xlabel('Participant ID')
ylabel('CE RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-20 180])
legend('1 session mean')
title('Flanker Task')
subplot(2,nTask,2)
violinplot(PP_CErtMat,[],'ViolinColor',palette(2,:))
xlabel('Participant ID')
ylabel('CE RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-20 180])
legend('1 session mean')
title('Prime-Probe Task')
subplot(2,nTask,3)
violinplot(ST_CErtMat,[],'ViolinColor',palette(3,:))
xlabel('Participant ID')
ylabel('CE RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-20 180])
legend('1 session mean')
title('Stroop Task')
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CErtGrandMean,FL_CErtConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(1,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('CE RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-20 180])
title('Flanker Task')
subplot(2,nTask,5)
a = errorbar(PP_CErtGrandMean,PP_CErtConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(2,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('CE RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-20 180])
title('Prime-Probe Task')
subplot(2,nTask,6)
a = errorbar(ST_CErtGrandMean,ST_CErtConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(3,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('CE RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-20 180])
title('Stroop Task')

% Acc
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CEaccMat,[],'ViolinColor',palette(1,:))
xlabel('Participant ID')
ylabel('CE accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
legend('1 session mean')
title('Flanker Task')
subplot(2,nTask,2)
violinplot(PP_CEaccMat,[],'ViolinColor',palette(2,:))
xlabel('Participant ID')
ylabel('CE accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
legend('1 session mean')
title('Prime-Probe Task')
subplot(2,nTask,3)
violinplot(ST_CEaccMat,[],'ViolinColor',palette(3,:))
xlabel('Participant ID')
ylabel('CE accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
legend('1 session mean')
title('Stroop Task')
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CEaccGrandMean,FL_CEaccConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(1,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('CE accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
title('Flanker Task')
subplot(2,nTask,5)
a = errorbar(PP_CEaccGrandMean,PP_CEaccConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(2,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('CE accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
title('Prime-Probe Task')
subplot(2,nTask,6)
a = errorbar(ST_CEaccGrandMean,ST_CEaccConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(3,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('CE accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
title('Stroop Task')
sgtitle('CE accuracy')

% IES
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CEiesMat,[],'ViolinColor',palette(1,:))
xlabel('Participant ID')
ylabel('CE IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-30 250])
legend('1 session mean')
title('Flanker Task')
subplot(2,nTask,2)
violinplot(PP_CEiesMat,[],'ViolinColor',palette(2,:))
xlabel('Participant ID')
ylabel('CE IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-30 250])
legend('1 session mean')
title('Prime-Probe Task')
subplot(2,nTask,3)
violinplot(ST_CEiesMat,[],'ViolinColor',palette(3,:))
xlabel('Participant ID')
ylabel('CE IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-30 250])
legend('1 session mean')
title('Stroop Task')
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CEiesGrandMean,FL_CEiesConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(1,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('CE IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-30 250])
title('Flanker Task')
subplot(2,nTask,5)
a = errorbar(PP_CEiesGrandMean,PP_CEiesConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(2,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('CE IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-30 250])
title('Prime-Probe Task')
subplot(2,nTask,6)
a = errorbar(ST_CEiesGrandMean,ST_CEiesConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(3,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('CE IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([-30 250])
title('Stroop Task')
sgtitle('CE IES')

%% (2) Congruent
% RT
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CrtMat,[],'ViolinColor',palette(1,:))
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([200 700])  % if you want the same length of range(400ms) across the three tasks
%legend('1 session mean')
title('Flanker Task')
subplot(2,nTask,2)
violinplot(PP_CrtMat,[],'ViolinColor',palette(2,:))
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(400ms) across the three tasks
%legend('1 session mean')
title('Prime-Probe Task')
subplot(2,nTask,3)
violinplot(ST_CrtMat,[],'ViolinColor',palette(3,:))
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(400ms) across the three tasks
%legend('1 session mean')
title('Stroop Task')
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CrtGrandMean,FL_CrtConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(1,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([200 700])  % if you want the same length of range(400ms) across the three tasks
title('Flanker Task')
subplot(2,nTask,5)
a = errorbar(PP_CrtGrandMean,PP_CrtConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(2,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(400ms) across the three tasks
title('Prime-Probe Task')
subplot(2,nTask,6)
a = errorbar(ST_CrtGrandMean,ST_CrtConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(3,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(400ms) across the three tasks
title('Stroop Task')
sgtitle('Congruent Trials RT')

% Acc
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CaccMat,[],'ViolinColor',palette(1,:))
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Flanker Task')
subplot(2,nTask,2)
violinplot(PP_CaccMat,[],'ViolinColor',palette(2,:))
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Prime-Probe Task')
subplot(2,nTask,3)
violinplot(ST_CaccMat,[],'ViolinColor',palette(3,:))
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Stroop Task')
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CaccGrandMean,FL_CaccConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(1,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
title('Flanker Task')
subplot(2,nTask,5)
a = errorbar(PP_CaccGrandMean,PP_CaccConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(2,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
title('Prime-Probe Task')
subplot(2,nTask,6)
a = errorbar(ST_CaccGrandMean,ST_CaccConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(3,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
title('Stroop Task')
sgtitle('Congruent Trials Accuracy')

% IES
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CiesMat,[],'ViolinColor',palette(1,:))
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Flanker Task')
subplot(2,nTask,2)
violinplot(PP_CiesMat,[],'ViolinColor',palette(2,:))
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Prime-Probe Task')
subplot(2,nTask,3)
violinplot(ST_CiesMat,[],'ViolinColor',palette(3,:))
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Stroop Task')
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CiesGrandMean,FL_CiesConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(1,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
title('Flanker Task')
subplot(2,nTask,5)
a = errorbar(PP_CiesGrandMean,PP_CiesConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(2,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
title('Prime-Probe Task')
subplot(2,nTask,6)
a = errorbar(ST_CiesGrandMean,ST_CiesConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(3,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
title('Stroop Task')
sgtitle('Congruent Trials IES')

%% (3) Incongruent
% RT
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_IrtMat,[],'ViolinColor',palette(1,:))
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([200 700])  % if you want the same length of range(500ms) across the three tasks
%legend('1 session mean')
title('Flanker Task')
subplot(2,nTask,2)
violinplot(PP_IrtMat,[],'ViolinColor',palette(2,:))
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(500ms) across the three tasks
%legend('1 session mean')
title('Prime-Probe Task')
subplot(2,nTask,3)
violinplot(ST_IrtMat,[],'ViolinColor',palette(3,:))
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(500ms) across the three tasks
%legend('1 session mean')
title('Stroop Task')
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_IrtGrandMean,FL_IrtConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(1,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([200 700])  % if you want the same length of range(500ms) across the three tasks
title('Flanker Task')
subplot(2,nTask,5)
a = errorbar(PP_IrtGrandMean,PP_IrtConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(2,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(500ms) across the three tasks
title('Prime-Probe Task')
subplot(2,nTask,6)
a = errorbar(ST_IrtGrandMean,ST_IrtConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(3,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('RT (ms)')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(500ms) across the three tasks
title('Stroop Task')
sgtitle('Incongruent Trials RT')

% Acc
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_IaccMat,[],'ViolinColor',palette(1,:))
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Flanker Task')
subplot(2,nTask,2)
violinplot(PP_IaccMat,[],'ViolinColor',palette(2,:))
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Prime-Probe Task')
subplot(2,nTask,3)
violinplot(ST_IaccMat,[],'ViolinColor',palette(3,:))
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Stroop Task')
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_IaccGrandMean,FL_IaccConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(1,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
title('Flanker Task')
subplot(2,nTask,5)
a = errorbar(PP_IaccGrandMean,PP_IaccConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(2,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
title('Prime-Probe Task')
subplot(2,nTask,6)
a = errorbar(ST_IaccGrandMean,ST_IaccConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(3,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('Accuracy')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([0.8 1])
title('Stroop Task')
sgtitle('Incongruent Trials Accuracy')

% IES
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_IiesMat,[],'ViolinColor',palette(1,:))
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Flanker Task')
subplot(2,nTask,2)
violinplot(PP_IiesMat,[],'ViolinColor',palette(2,:))
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Prime-Probe Task')
subplot(2,nTask,3)
violinplot(ST_IiesMat,[],'ViolinColor',palette(3,:))
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Stroop Task')
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_IiesGrandMean,FL_IiesConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(1,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
title('Flanker Task')
subplot(2,nTask,5)
a = errorbar(PP_IiesGrandMean,PP_IiesConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(2,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000]) 
title('Prime-Probe Task')
subplot(2,nTask,6)
a = errorbar(ST_IiesGrandMean,ST_IiesConHalf,'LineStyle','none');
a.Marker = 'o';
a.MarkerSize = 2;
a.Color = palette(3,:);
a.CapSize = 15;
xlabel('Participant ID')
ylabel('IES')
xlim([0 9])
xticks(1:9)
xticklabels(subjID(1:end))
ylim([300 1000])
title('Stroop Task')
sgtitle('Incongruent Trials IES')