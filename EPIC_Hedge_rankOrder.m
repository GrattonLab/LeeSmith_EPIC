%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
%
% To run this sciprt:
% You need Hedge et al.'s (2018) flanker and Stroop task data 
% (download at https://osf.io/cwzds/)
% You also need ICC.m (download at
% https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
% What this script does:
% Plots the rank correlation coefficients between CE and incongruent trial
%
% What this script outputs:
% Supp. Fig. 5B
%
% Created on 05/16/2023 by HJ Lee
% Last modified on 07/03/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')

%% Parameter settings
nCond = 2;  % number of experimental conditions: congruent/incongruent 
nStd = 3;  % STD from mean for excluding outliers

tmpSubjID1 = 1:50;
tmpSubjID2 = 1:62;
subjExcld1 = [6, 17, 37, 21, 33];  % based on Hedge et al.'s Readme.txt; also removed 21 33 based on accuracy
subjExcld2 = [25, 28, 34, 38, 56, 32];  % based on Readme.txt; also removed 32
subjID1 = setdiff(tmpSubjID1, subjExcld1);
subjID2 = setdiff(tmpSubjID2, subjExcld2);
nSubj1 = length(subjID1);
nSubj2 = length(subjID2);
nSubjT = nSubj1+nSubj2;
nTask = 1;  % just use the flanker task data
taskStrng1 = {'Flanker','Stroop'};
%nSession = 2;
nBlocks = 5;
nTrials = 144;  % per block
%nTrialsReal = 96;  % actual trials for calculating CE (excluding neutral)

drawLH = 4;  % ~nTdrawn; the dataset only has about 480 trials

%% 1. Load and process Hedge et al.'s (2018) data
GrtMatH = nan(nCond,nSubjT);  % grand mean RT
GaccMatH = nan(nCond,nSubjT);  % Acc
% Study 1 data
for i = 1:nSubj1
    for t = 1:nTask  % 1
        dir = ['HedgeData/Study1-' taskStrng1{t}];
        cd(dir)
        % session 1
        fileID  = ['Study1_P' num2str(subjID1(i)) taskStrng1{t} num2str(1) '.csv'];
        rawData = readmatrix(fileID);
        T1 = table((1:nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T1(T1.n0cong==1,:) = [];  % exclude neutral condition
        % session 2
        fileID  = ['Study1_P' num2str(subjID1(i)) taskStrng1{t} num2str(2) '.csv'];
        rawData = readmatrix(fileID);
        T2 = table((1+nBlocks*nTrials:2*nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T2(T2.n0cong==1,:) = [];  % exclude neutral condition
        T = [T1;T2];  % combine session 1&2
        cd('../../')  % home directory
        % Exclude outlier trials
        % T
        acc_1_ids = find(T.acc==1);
        if length(find(T.acc==1))/length(table2array(T))*100 < 70
            disp(['Study1 ' taskStrng1{t} ' task subject#' num2str(subjID1(i)) ' has below 70% accuracy'])
        end
        T_rm = T(acc_1_ids,:);  % accurate trials RT
        TGM = varfun(@mean,T_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS = varfun(@std,T_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        inbound = TGM.mean_rt-(TGS.std_rt*nStd);  % lower bound
        for i_inbound = 1:size(inbound,1)
            if inbound(i_inbound) < 0
                inbound(i_inbound) = 0;
            end
        end
        %inbound = 0.15;  % alternatively
        outbound = TGM.mean_rt+(TGS.std_rt*nStd);  % upper bound
        out_ids = [];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==0 & T_rm.rt<inbound(1))];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==2 & T_rm.rt<inbound(2))];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==0 & T_rm.rt>outbound(1))];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==2 & T_rm.rt>outbound(2))];
        T_rm.acc(ismember(T_rm.trialNum,out_ids)) = 99;
        T.acc(ismember(T.trialNum,out_ids)) = 99;
        % Calculate mean RT (secs)
        TGMrt = varfun(@mean,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGSrt = varfun(@std,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        GrtMatH(:,i) = TGMrt.mean_rt;
        % Calculate mean Acc
        TGM_total = varfun(@mean,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGM_corr = varfun(@mean,T(T.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS_total = varfun(@std,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        GaccMatH(1,i) = TGM_corr.GroupCount(1)/TGM_total.GroupCount(1);
        GaccMatH(2,i) = TGM_corr.GroupCount(2)/TGM_total.GroupCount(2);        
    end
end
% Study 2 data
for i = 1:nSubj2
    for t = 1:nTask
        dir = ['HedgeData/Study2-' taskStrng1{t}];
        cd(dir)
        % session 1
        fileID  = ['Study2_P' num2str(subjID2(i)) taskStrng1{t} num2str(1) '.csv'];
        rawData = readmatrix(fileID);
        T1 = table((1:nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T1(T1.n0cong==1,:) = [];  % exclude neutral condition
        % session 2
        fileID  = ['Study2_P' num2str(subjID2(i)) taskStrng1{t} num2str(2) '.csv'];
        rawData = readmatrix(fileID);
        T2 = table((1+nBlocks*nTrials:2*nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T2(T2.n0cong==1,:) = [];  % exclude neutral condition
        T = [T1;T2];
        cd('../../')  % home directory
        % Exclude outlier trials
        % T
        acc_1_ids = find(T.acc==1);
        if length(find(T.acc==1))/length(table2array(T))*100 < 70
            disp(['Study2 ' taskStrng1{t} ' task subject#' num2str(subjID2(i)) ' has below 70% accuracy'])
        end
        T_rm = T(acc_1_ids,:);
        TGM = varfun(@mean,T_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS = varfun(@std,T_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        inbound = TGM.mean_rt-(TGS.std_rt*nStd);  % lower bound
        for i_inbound = 1:size(inbound,1)
            if inbound(i_inbound) < 0
                inbound(i_inbound) = 0;
            end
        end
        %inbound = 0.15;  % alternatively
        outbound = TGM.mean_rt+(TGS.std_rt*nStd);  % higher bound
        out_ids = [];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==0 & T_rm.rt<inbound(1))];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==2 & T_rm.rt<inbound(2))];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==0 & T_rm.rt>outbound(1))];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==2 & T_rm.rt>outbound(2))];
        T_rm.acc(ismember(T_rm.trialNum,out_ids)) = 99;
        T.acc(ismember(T.trialNum,out_ids)) = 99;
        % Calculate mean RT (secs)
        TGMrt = varfun(@mean,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGSrt = varfun(@std,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        GrtMatH(:,nSubj1+i) = TGMrt.mean_rt;
        % Calculate mean Acc
        TGM_total = varfun(@mean,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGM_corr = varfun(@mean,T(T.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS_total = varfun(@std,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        GaccMatH(1,nSubj1+i) = TGM_corr.GroupCount(1)/TGM_total.GroupCount(1);
        GaccMatH(2,nSubj1+i) = TGM_corr.GroupCount(2)/TGM_total.GroupCount(2);
    end
end

%% Plot the correlation between CE and incongruent trials
% (1) RT
[~,p1] = sort(GrtMatH(2,:)-GrtMatH(1,:),'descend');
r1 = 1:nSubjT;
r1(p1) = r1;
[~,p2] = sort(GrtMatH(2,:),'descend');
r2 = 1:nSubjT;
r2(p2) = r2;
% Kendall's tau
tauH = corr(r1',r2','type','kendall');
iccH = ICC([r1',r2'],'A-k',0.05);
figure
scatter(GrtMatH(2,:),GrtMatH(2,:)-GrtMatH(1,:),'filled')
%lsline
refline
axis square
set(gca,'FontSize',13)
str = sprintf('  Tau = %1.2f\n  ICC = %1.2f', tauH, iccH);
rText = text(min(get(gca,'xlim')),max(get(gca,'ylim')), str);
set(rText, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
xlabel('Incongruent trial RT (sec)','FontSize',14)
ylabel('CE RT (sec)','FontSize',14)
title('Hedge et al.''s (2018) data','FontSize',15)
% (2) Accuracy
GpeMatH = 100.*(1-GaccMatH);
[~,p1] = sort(GpeMatH(2,:)-GpeMatH(1,:),'descend');
r1 = 1:nSubjT;
r1(p1) = r1;
[~,p2] = sort(GpeMatH(2,:),'descend');
r2 = 1:nSubjT;
r2(p2) = r2;
% Kendall's tau
tauH = corr(r1',r2','type','kendall');
iccH = ICC([r1',r2'],'A-k',0.05);
figure
scatter(GpeMatH(2,:),GpeMatH(2,:)-GpeMatH(1,:),'filled')
%lsline
refline
axis square
set(gca,'FontSize',13)
str = sprintf('  Tau = %1.2f\n  ICC = %1.2f', tauH, iccH);
rText = text(min(get(gca,'xlim')),max(get(gca,'ylim')), str);
set(rText, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
xlabel('Incongruent trial percent correct (%)','FontSize',14)
ylabel('CE percent correct (%)','FontSize',14)
title('Hedge et al.''s (2018) data','FontSize',15)