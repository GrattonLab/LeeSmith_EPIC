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
% Fig. 2: Reaction time and accuracy for congruency effect (CE)
% Supp. Fig. 1: CE inverse efficiency score (IES)
% Supp. Fig. 2: Reaction time for congruent/incongruent trials
% Supp. Fig. 3: Accuracy for congruent/incongruent trials
% Supp. Fig. 4: IES for congruent/incongruent trials
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
subjExcld = [1,2,11];  % excluded participants
subjID = setdiff(tmpSubjID, subjExcld);
nSubj = length(subjID);

nTask = 3;
palette = parula(nTask);
palette(3,:) = [1 0 0];  % red

%% Load processed data
load FL_CEmat.mat
load PP_CEmat.mat
load ST_CEmat.mat

% Remove participant #9 if exists
if nSubj == 9
    rmIdx = find(subjID==9);
    % FL
    FL_CErtMat(:,rmIdx) = [];
    FL_CEaccMat(:,rmIdx) = [];
    FL_CEiesMat(:,rmIdx) = [];
    FL_CrtMat(:,rmIdx) = [];
    FL_IrtMat(:,rmIdx) = [];
    FL_CaccMat(:,rmIdx) = [];
    FL_IaccMat(:,rmIdx) = [];
    FL_CiesMat(:,rmIdx) = [];
    FL_IiesMat(:,rmIdx) = [];

    FL_CErtGrandMean(rmIdx) = [];
    FL_CEaccGrandMean(rmIdx) = [];
    FL_CEiesGrandMean(rmIdx) = [];
    FL_CrtGrandMean(rmIdx) = [];
    FL_IrtGrandMean(rmIdx) = [];
    FL_CaccGrandMean(rmIdx) = [];
    FL_IaccGrandMean(rmIdx) = [];
    FL_CiesGrandMean(rmIdx) = [];
    FL_IiesGrandMean(rmIdx) = [];

    FL_CErtConHalf(rmIdx) = [];
    FL_CEaccConHalf(rmIdx) = [];
    FL_CEiesConHalf(rmIdx) = [];
    FL_CrtConHalf(rmIdx) = [];
    FL_IrtConHalf(rmIdx) = [];
    FL_CaccConHalf(rmIdx) = [];
    FL_IaccConHalf(rmIdx) = [];
    FL_CiesConHalf(rmIdx) = [];
    FL_IiesConHalf(rmIdx) = [];

    % PP
    PP_CErtMat(:,rmIdx) = [];
    PP_CEaccMat(:,rmIdx) = [];
    PP_CEiesMat(:,rmIdx) = [];
    PP_CrtMat(:,rmIdx) = [];
    PP_IrtMat(:,rmIdx) = [];
    PP_CaccMat(:,rmIdx) = [];
    PP_IaccMat(:,rmIdx) = [];
    PP_CiesMat(:,rmIdx) = [];
    PP_IiesMat(:,rmIdx) = [];

    PP_CErtGrandMean(rmIdx) = [];
    PP_CEaccGrandMean(rmIdx) = [];
    PP_CEiesGrandMean(rmIdx) = [];
    PP_CrtGrandMean(rmIdx) = [];
    PP_IrtGrandMean(rmIdx) = [];
    PP_CaccGrandMean(rmIdx) = [];
    PP_IaccGrandMean(rmIdx) = [];
    PP_CiesGrandMean(rmIdx) = [];
    PP_IiesGrandMean(rmIdx) = [];

    PP_CErtConHalf(rmIdx) = [];
    PP_CEaccConHalf(rmIdx) = [];
    PP_CEiesConHalf(rmIdx) = [];
    PP_CrtConHalf(rmIdx) = [];
    PP_IrtConHalf(rmIdx) = [];
    PP_CaccConHalf(rmIdx) = [];
    PP_IaccConHalf(rmIdx) = [];
    PP_CiesConHalf(rmIdx) = [];
    PP_IiesConHalf(rmIdx) = [];

    % ST
    ST_CErtMat(:,rmIdx) = [];
    ST_CEaccMat(:,rmIdx) = [];
    ST_CEiesMat(:,rmIdx) = [];
    ST_CrtMat(:,rmIdx) = [];
    ST_IrtMat(:,rmIdx) = [];
    ST_CaccMat(:,rmIdx) = [];
    ST_IaccMat(:,rmIdx) = [];
    ST_CiesMat(:,rmIdx) = [];
    ST_IiesMat(:,rmIdx) = [];

    ST_CErtGrandMean(rmIdx) = [];
    ST_CEaccGrandMean(rmIdx) = [];
    ST_CEiesGrandMean(rmIdx) = [];
    ST_CrtGrandMean(rmIdx) = [];
    ST_IrtGrandMean(rmIdx) = [];
    ST_CaccGrandMean(rmIdx) = [];
    ST_IaccGrandMean(rmIdx) = [];
    ST_CiesGrandMean(rmIdx) = [];
    ST_IiesGrandMean(rmIdx) = [];

    ST_CErtConHalf(rmIdx) = [];
    ST_CEaccConHalf(rmIdx) = [];
    ST_CEiesConHalf(rmIdx) = [];
    ST_CrtConHalf(rmIdx) = [];
    ST_IrtConHalf(rmIdx) = [];
    ST_CaccConHalf(rmIdx) = [];
    ST_IaccConHalf(rmIdx) = [];
    ST_CiesConHalf(rmIdx) = [];
    ST_IiesConHalf(rmIdx) = [];
    
    subjID(rmIdx) = [];
    nSubj = nSubj-1;
end
%% Violin plot and Grand mean plot together
% (1) CE
% Reaction time
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CErtMat,[],'ViolinColor',palette(1,:))
set(gca,'FontSize',15,'FontWeight','bold')  % 12
%xlabel('Participant ID','FontSize',16)
ylabel({'Congruency effect'; 'reaction time (ms)'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-20 180])
%legend('1 session mean','FontSize',13,'Location','best')
title('Flanker Task','FontSize',18)
subplot(2,nTask,2)
violinplot(PP_CErtMat,[],'ViolinColor',palette(2,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'reaction time (ms)'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-20 180])
%legend('1 session mean','FontSize',13,'Location','best')
title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,3)
violinplot(ST_CErtMat,[],'ViolinColor',palette(3,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'reaction time (ms)'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-20 180])
%legend('1 session mean','FontSize',13,'Location','best')
title('Stroop Task','FontSize',18)
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CErtGrandMean,FL_CErtConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(1,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
ylabel({'Congruency effect'; 'reaction time (ms)'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-20 180])
%title('Flanker Task','FontSize',18)
subplot(2,nTask,5)
a = errorbar(PP_CErtGrandMean,PP_CErtConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 2
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(2,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'reaction time (ms)'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-20 180])
%title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,6)
a = errorbar(ST_CErtGrandMean,ST_CErtConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 2
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(3,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'reaction time (ms)'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-20 180])
%title('Stroop Task','FontSize',18)

%% Acc
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CEaccMat,[],'ViolinColor',palette(1,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
ylabel({'Congruency effect'; 'accuracy'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
%legend('1 session mean','FontSize',13,'Location','best')
title('Flanker Task','FontSize',18)
subplot(2,nTask,2)
violinplot(PP_CEaccMat,[],'ViolinColor',palette(2,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'accuracy'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
%legend('1 session mean','FontSize',13,'Location','best')
title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,3)
violinplot(ST_CEaccMat,[],'ViolinColor',palette(3,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'accuracy'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
%legend('1 session mean','FontSize',13,'Location','best')
title('Stroop Task','FontSize',18)
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CEaccGrandMean,FL_CEaccConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(1,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
ylabel({'Congruency effect'; 'accuracy'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
%title('Flanker Task','FontSize',18)
subplot(2,nTask,5)
a = errorbar(PP_CEaccGrandMean,PP_CEaccConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(2,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'accuracy'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
%title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,6)
a = errorbar(ST_CEaccGrandMean,ST_CEaccConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(3,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'accuracy'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-0.1 0.3])
%title('Stroop Task','FontSize',18)
%sgtitle('Congruency effect accuracy')

%% IES
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CEiesMat,[],'ViolinColor',palette(1,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
ylabel({'Congruency effect'; 'IES'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-30 250])
%legend('1 session mean','FontSize',13,'Location','northeast')
title('Flanker Task','FontSize',18)
subplot(2,nTask,2)
violinplot(PP_CEiesMat,[],'ViolinColor',palette(2,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'IES'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-30 250])
%legend('1 session mean','FontSize',13,'Location','northeast')
title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,3)
violinplot(ST_CEiesMat,[],'ViolinColor',palette(3,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'IES'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-30 250])
%legend('1 session mean','FontSize',13,'Location','northeast')
title('Stroop Task','FontSize',18)
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CEiesGrandMean,FL_CEiesConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(1,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
ylabel({'Congruency effect'; 'IES'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-30 250])
%title('Flanker Task','FontSize',18)
subplot(2,nTask,5)
a = errorbar(PP_CEiesGrandMean,PP_CEiesConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(2,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'IES'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-30 250])
%title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,6)
a = errorbar(ST_CEiesGrandMean,ST_CEiesConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(3,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel({'Congruency effect'; 'IES'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([-30 250])
%title('Stroop Task','FontSize',18)
%sgtitle('Congruency effect IES')

%% (2) Congruent
% Reaction time
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CrtMat,[],'ViolinColor',palette(1,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
ylabel('Reaction time (ms)','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([200 700])  % if you want the same length of range(400ms) across the three tasks
%legend('1 session mean')
title('Flanker Task','FontSize',18)
subplot(2,nTask,2)
violinplot(PP_CrtMat,[],'ViolinColor',palette(2,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel({'Congruent trial';'reaction time (ms)'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(400ms) across the three tasks
%legend('1 session mean')
title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,3)
violinplot(ST_CrtMat,[],'ViolinColor',palette(3,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel({'Congruent trial';'reaction time (ms)'},'FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(400ms) across the three tasks
%legend('1 session mean')
title('Stroop Task','FontSize',18)
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CrtGrandMean,FL_CrtConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(1,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
ylabel('Reaction time (ms)','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([200 700])  % if you want the same length of range(400ms) across the three tasks
%title('Flanker Task','FontSize',18)
subplot(2,nTask,5)
a = errorbar(PP_CrtGrandMean,PP_CrtConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(2,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('Reaction time (ms)','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(400ms) across the three tasks
%title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,6)
a = errorbar(ST_CrtGrandMean,ST_CrtConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(3,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('Reaction time (ms)','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([200 800])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(400ms) across the three tasks
%title('Stroop Task','FontSize',18)
%sgtitle('Congruent Trials Reaction time')

%% Acc
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CaccMat,[],'ViolinColor',palette(1,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Flanker Task','FontSize',18)
subplot(2,nTask,2)
violinplot(PP_CaccMat,[],'ViolinColor',palette(2,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,3)
violinplot(ST_CaccMat,[],'ViolinColor',palette(3,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Stroop Task','FontSize',18)
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CaccGrandMean,FL_CaccConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(1,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%title('Flanker Task','FontSize',18)
subplot(2,nTask,5)
a = errorbar(PP_CaccGrandMean,PP_CaccConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(2,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,6)
a = errorbar(ST_CaccGrandMean,ST_CaccConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(3,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%title('Stroop Task','FontSize',18)
%sgtitle('Congruent Trials Accuracy')

%% IES
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_CiesMat,[],'ViolinColor',palette(1,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Flanker Task','FontSize',18)
subplot(2,nTask,2)
violinplot(PP_CiesMat,[],'ViolinColor',palette(2,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,3)
violinplot(ST_CiesMat,[],'ViolinColor',palette(3,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Stroop Task','FontSize',18)
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_CiesGrandMean,FL_CiesConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(1,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%title('Flanker Task','FontSize',18)
subplot(2,nTask,5)
a = errorbar(PP_CiesGrandMean,PP_CiesConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(2,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,6)
a = errorbar(ST_CiesGrandMean,ST_CiesConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(3,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%title('Stroop Task','FontSize',18)
%sgtitle('Congruent Trials IES')

%% (3) Incongruent
% Reaction time
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_IrtMat,[],'ViolinColor',palette(1,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
ylabel('Reaction time (ms)','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([200 700])  % if you want the same length of range(500ms) across the three tasks
%legend('1 session mean')
title('Flanker Task','FontSize',18)
subplot(2,nTask,2)
violinplot(PP_IrtMat,[],'ViolinColor',palette(2,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel('Reaction time (ms)','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(500ms) across the three tasks
%legend('1 session mean')
title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,3)
violinplot(ST_IrtMat,[],'ViolinColor',palette(3,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel('Reaction time (ms)','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(500ms) across the three tasks
%legend('1 session mean')
title('Stroop Task','FontSize',18)
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_IrtGrandMean,FL_IrtConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(1,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
ylabel('Reaction time (ms)','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([200 700])  % if you want the same length of range(500ms) across the three tasks
%title('Flanker Task','FontSize',18)
subplot(2,nTask,5)
a = errorbar(PP_IrtGrandMean,PP_IrtConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(2,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('Reaction time (ms)','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(500ms) across the three tasks
%title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,6)
a = errorbar(ST_IrtGrandMean,ST_IrtConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(3,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('Reaction time (ms)','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
%ylim([300 900])  % if you want the same range across the three tasks
ylim([400 900])  % if you want the same length of range(500ms) across the three tasks
%title('Stroop Task','FontSize',18)
%sgtitle('Incongruent Trials Reaction time')

%% Acc
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_IaccMat,[],'ViolinColor',palette(1,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Flanker Task','FontSIze',18)
subplot(2,nTask,2)
violinplot(PP_IaccMat,[],'ViolinColor',palette(2,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,3)
violinplot(ST_IaccMat,[],'ViolinColor',palette(3,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%legend('1 session mean')
title('Stroop Task','FontSize',18)
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_IaccGrandMean,FL_IaccConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(1,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%title('Flanker Task','FontSize',18)
subplot(2,nTask,5)
a = errorbar(PP_IaccGrandMean,PP_IaccConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(2,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,6)
a = errorbar(ST_IaccGrandMean,ST_IaccConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(3,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('Accuracy','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([0.8 1])
%title('Stroop Task','FontSize',18)
%sgtitle('Incongruent Trials Accuracy')

%% IES
figure
% Violin plot
subplot(2,nTask,1)
violinplot(FL_IiesMat,[],'ViolinColor',palette(1,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Flanker Task','FontSize',18)
subplot(2,nTask,2)
violinplot(PP_IiesMat,[],'ViolinColor',palette(2,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,3)
violinplot(ST_IiesMat,[],'ViolinColor',palette(3,:))
set(gca,'FontSize',15,'FontWeight','bold')
%xlabel('Participant ID','FontSize',16)
%ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%legend('1 session mean')
title('Stroop Task','FontSize',18)
% Grand mean plot
subplot(2,nTask,4)
a = errorbar(FL_IiesGrandMean,FL_IiesConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(1,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%title('Flanker Task','FontSize',18)
subplot(2,nTask,5)
a = errorbar(PP_IiesGrandMean,PP_IiesConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(2,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000]) 
%title('Prime-Probe Task','FontSize',18)
subplot(2,nTask,6)
a = errorbar(ST_IiesGrandMean,ST_IiesConHalf,'LineStyle','none');
a.Marker = '.';
a.MarkerSize = 20;  % 12
a.Bar.LineWidth = 2;
a.Cap.LineWidth = 2;
a.Color = palette(3,:);
a.CapSize = 15;
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Participant ID','FontSize',16)
%ylabel('IES','FontSize',16)
xlim([0 nSubj+1])
xticks(1:nSubj+1)
xticklabels(subjID(1:end))
ylim([300 1000])
%title('Stroop Task','FontSize',18)
%sgtitle('Incongruent Trials IES')