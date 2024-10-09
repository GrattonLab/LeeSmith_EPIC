%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
%
% To run this script:
% You need Hedge et al.'s (2018) flanker task data 
% (download at https://osf.io/cwzds/),
% a Matlab function named "ezdiffusion.m" (download with EPIC data), and
% ICC.m (download at
% https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
% What this script does:
% Simulates Hedge et al.'s (2018) flanker task data and performs
% EZ-diffusion modeling to examine how sampling size affects the
% reliability of the modeling parameters
% (RT needs to be in SECS to run EZ-diffusion modeling)
%
% What this script outputs:
% Supp. Fig. 20
% Plot of the split-half reliability of the parameters,
% drift rate, boundary separation, and nondecision time
%
% Created on 05/16/2023 by HJ Lee
% Last modified on 07/05/2023
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
subjExcld1 = [6, 17, 37, 21, 33];  % based on Hedge et al.'s Readme.txt; also removed 21 33
subjExcld2 = [25, 28, 34, 38, 56, 32];  % based on Readme.txt; also removed 32
subjID1 = setdiff(tmpSubjID1, subjExcld1);
subjID2 = setdiff(tmpSubjID2, subjExcld2);
nSubj1 = length(subjID1);  % study 1
nSubj2 = length(subjID2);  % study 2
nSubjT = nSubj1+nSubj2;
nTask = 1;  % just use the flanker task data for EZdiffusion; otherwise 2
taskStrng = {'Flanker','Stroop'};
%nSession = 2;
nBlocks = 5;
nTrials = 144;  % per block
%nTrialsCE = 96;  % number of trials for calculating CE (excluding neutral)

nTdrawn = [50, 100, 200, 400, 800, 1600 3200]/2;  % number of trials drawn (per condition)
drawLH = 4;  % number of trials drawn for Hedge et al.'s data; max 480 trials (cf. nTdrawn)
n_subs = 100;
bs_dist_n = 5000;
ws_dist_n = 10000;
numTest = 100;
nzsigr = 0.24;  % noise sigma; determined by evaluating SSE of ICCs; run noisesigmaEvaluation script; 0.2-0.3
nzsiga = 0.22;

%% 1. Load and preprocess Hedge et al.'s (2018) flanker task data
GrtMat = nan(nCond,nSubjT);  % grand mean RT; need this for simulation (copulas)
GaccMat = nan(nCond,nSubjT);  % accuracy
stdMatr = nan(nCond,nSubjT);  % within-subject STD; need this also for simulation
stdMata = nan(nCond,nSubjT);
HrtMat1 = nan(nCond,length(nTdrawn),nSubjT);  % mean RT in secs; split-half grand mean; need this to compare w/ simulated data
HrtMat2 = nan(nCond,length(nTdrawn),nSubjT);
HrtvMat1 = nan(nCond,length(nTdrawn),nSubjT);  % RT variance
HrtvMat2 = nan(nCond,length(nTdrawn),nSubjT);
HaccMat1 = nan(nCond,length(nTdrawn),nSubjT);  % mean accuracy
HaccMat2 = nan(nCond,length(nTdrawn),nSubjT);
% Study 1 data
for i = 1:nSubj1
    for t = 1:nTask  % 1
        dir = ['HedgeData/Study1-' taskStrng{t}];
        cd(dir)
        % session 1
        fileID  = ['Study1_P' num2str(subjID1(i)) taskStrng{t} num2str(1) '.csv'];
        rawData = readmatrix(fileID);
        T1 = table((1:nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T1(T1.n0cong==1,:) = [];  % exclude neutral condition
        % session 2
        fileID  = ['Study1_P' num2str(subjID1(i)) taskStrng{t} num2str(2) '.csv'];
        rawData = readmatrix(fileID);
        T2 = table((1+nBlocks*nTrials:2*nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T2(T2.n0cong==1,:) = [];  % exclude neutral condition
        T = [T1;T2];  % combine session 1&2
        cd('../../')  % home directory
        % Exclude outlier trials
        % T1
        T1_rm = T1(T1.acc==1,:);
        TGM1 = varfun(@mean,T1_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS1 = varfun(@std,T1_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        inbound = TGM1.mean_rt-(TGS1.std_rt*nStd);  % 3SD
        for i_inbound = 1:size(inbound,1)
            if inbound(i_inbound) < 0
                inbound(i_inbound) = 0;
            end
        end
        %inbound = 0.15;  % alternatively
        outbound = TGM1.mean_rt+(TGS1.std_rt*nStd);
        out_ids = [];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==0 & T1_rm.rt<inbound(1))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==2 & T1_rm.rt<inbound(2))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==0 & T1_rm.rt>outbound(1))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==2 & T1_rm.rt>outbound(2))];
        T1_rm.acc(ismember(T1_rm.trialNum,out_ids)) = 99;
        T1.acc(ismember(T1.trialNum,out_ids)) = 99;
        % T2
        T2_rm = T2(T2.acc==1,:);
        TGM2 = varfun(@mean,T2_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS2 = varfun(@std,T2_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        inbound = TGM2.mean_rt-(TGS2.std_rt*nStd);
        for i_inbound = 1:size(inbound,1)
            if inbound(i_inbound) < 0
                inbound(i_inbound) = 0;
            end
        end
        %inbound = 0.15;  % alternatively
        outbound = TGM2.mean_rt+(TGS2.std_rt*nStd);  % upper bound
        out_ids = [];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==0 & T2_rm.rt<inbound(1))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==2 & T2_rm.rt<inbound(2))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==0 & T2_rm.rt>outbound(1))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==2 & T2_rm.rt>outbound(2))];
        T2_rm.acc(ismember(T2_rm.trialNum,out_ids)) = 99;
        T2.acc(ismember(T2.trialNum,out_ids)) = 99;
        % T
        acc_1_ids = find(T.acc==1);
        if length(find(T.acc==1))/length(table2array(T))*100 < 70
            disp(['Study1 ' taskStrng{t} ' task subject#' num2str(subjID1(i)) ' has below 70% accuracy'])
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
        GrtMat(:,i) = TGMrt.mean_rt;
        stdMatr(:,i) = TGSrt.std_rt;
        % Calculate mean Acc
        TGM_total = varfun(@mean,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGM_corr = varfun(@mean,T(T.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS_total = varfun(@std,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        GaccMat(1,i) = TGM_corr.GroupCount(1)/TGM_total.GroupCount(1);
        GaccMat(2,i) = TGM_corr.GroupCount(2)/TGM_total.GroupCount(2);
        stdMata(:,i) = TGS_total.std_acc;

        % descriptive
        rtC1 = T1.rt(and(T1.n0cong==0,T1.acc==1));  % con
        rtI1 = T1.rt(and(T1.n0cong==2,T1.acc==1));  % inc
        rtC2 = T2.rt(and(T2.n0cong==0,T2.acc==1));
        rtI2 = T2.rt(and(T2.n0cong==2,T2.acc==1));
        accC1 = T1.acc(T1.n0cong==0);
        accI1 = T1.acc(T1.n0cong==2);
        accC2 = T2.acc(T2.n0cong==0);
        accI2 = T2.acc(T2.n0cong==2);
        accC1(find(accC1==99)) = [];  % removing outlier trials
        accI1(find(accI1==99)) = [];
        accC2(find(accC2==99)) = [];
        accI2(find(accI2==99)) = [];
        for k = 1:drawLH
            if length(rtC1) < nTdrawn(k)  % this happens when k == 4
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtC1))])
                HrtMat1(1,k,i) = mean(rtC1);  % don't need trimmean because outlier trials were marked as 99 above
                HrtvMat1(1,k,i) = var(rtC1);
                HaccMat1(1,k,i) = mean(accC1);
            else
                HrtMat1(1,k,i) = mean(rtC1(1:nTdrawn(k)));
                HrtvMat1(1,k,i) = var(rtC1(1:nTdrawn(k)));
                HaccMat1(1,k,i) = mean(accC1(1:nTdrawn(k)));
            end
            if length(rtI1) < nTdrawn(k)
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtI1))])
                HrtMat1(2,k,i) = mean(rtI1);
                HrtvMat1(2,k,i) = var(rtI1);
                HaccMat1(2,k,i) = mean(accI1);
            else
                HrtMat1(2,k,i) = mean(rtI1(1:nTdrawn(k)));
                HrtvMat1(2,k,i) = var(rtI1(1:nTdrawn(k)));
                HaccMat1(2,k,i) = mean(accI1(1:nTdrawn(k)));
            end
            if length(rtC2) < nTdrawn(k)
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtC2))])
                HrtMat2(1,k,i) = mean(rtC2);
                HrtvMat2(1,k,i) = var(rtC2);
                HaccMat2(1,k,i) = mean(accC2);
            else
                HrtMat2(1,k,i) = mean(rtC2(1:nTdrawn(k)));
                HrtvMat2(1,k,i) = var(rtC2(1:nTdrawn(k)));
                HaccMat2(1,k,i) = mean(accC2(1:nTdrawn(k)));
            end
            if length(rtI2) < nTdrawn(k)
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtI2))])
                HrtMat2(2,k,i) = mean(rtI2);
                HrtvMat2(2,k,i) = var(rtI2);
                HaccMat2(2,k,i) = mean(accI2);
            else
                HrtMat2(2,k,i) = mean(rtI2(1:nTdrawn(k)));
                HrtvMat2(2,k,i) = var(rtI2(1:nTdrawn(k)));
                HaccMat2(2,k,i) = mean(accI2(1:nTdrawn(k)));
            end
        end
    end
end
% Study 2 data
for i = 1:nSubj2
    for t = 1:nTask
        dir = ['HedgeData/Study2-' taskStrng{t}];
        cd(dir)
        % session 1
        fileID  = ['Study2_P' num2str(subjID2(i)) taskStrng{t} num2str(1) '.csv'];
        rawData = readmatrix(fileID);
        T1 = table((1:nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T1(T1.n0cong==1,:) = [];  % exclude neutral condition
        % session 2
        fileID  = ['Study2_P' num2str(subjID2(i)) taskStrng{t} num2str(2) '.csv'];
        rawData = readmatrix(fileID);
        T2 = table((1+nBlocks*nTrials:2*nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T2(T2.n0cong==1,:) = [];  % exclude neutral condition
        T = [T1;T2];
        cd('../../')  % home directory
        % Exclude outlier trials
        % T1
        T1_rm = T1(T1.acc==1,:);
        TGM1 = varfun(@mean,T1_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS1 = varfun(@std,T1_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        inbound = TGM1.mean_rt-(TGS1.std_rt*nStd);
        for i_inbound = 1:size(inbound,1)
            if inbound(i_inbound) < 0
                inbound(i_inbound) = 0;
            end
        end
        %inbound = 0.15;  % alternatively
        outbound = TGM1.mean_rt+(TGS1.std_rt*nStd);  % upper bound
        out_ids = [];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==0 & T1_rm.rt<inbound(1))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==2 & T1_rm.rt<inbound(2))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==0 & T1_rm.rt>outbound(1))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==2 & T1_rm.rt>outbound(2))];
        T1_rm.acc(ismember(T1_rm.trialNum,out_ids)) = 99;
        T1.acc(ismember(T1.trialNum,out_ids)) = 99;
        % T2
        T2_rm = T2(T2.acc==1,:);
        TGM2 = varfun(@mean,T2_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS2 = varfun(@std,T2_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        inbound = TGM2.mean_rt-(TGS2.std_rt*nStd);
        for i_inbound = 1:size(inbound,1)
            if inbound(i_inbound) < 0
                inbound(i_inbound) = 0;
            end
        end
        %inbound = 0.15;  % alternatively
        outbound = TGM2.mean_rt+(TGS2.std_rt*nStd);  % upper bound
        out_ids = [];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==0 & T2_rm.rt<inbound(1))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==2 & T2_rm.rt<inbound(2))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==0 & T2_rm.rt>outbound(1))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==2 & T2_rm.rt>outbound(2))];
        T2_rm.acc(ismember(T2_rm.trialNum,out_ids)) = 99;
        T2.acc(ismember(T2.trialNum,out_ids)) = 99;
        % T
        acc_1_ids = find(T.acc==1);
        if length(find(T.acc==1))/length(table2array(T))*100 < 70
            disp(['Study2 ' taskStrng{t} ' task subject#' num2str(subjID2(i)) ' has below 70% accuracy'])
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
        GrtMat(:,nSubj1+i) = TGMrt.mean_rt;
        stdMatr(:,nSubj1+i) = TGSrt.std_rt;
        % Calculate mean Acc
        TGM_total = varfun(@mean,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGM_corr = varfun(@mean,T(T.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS_total = varfun(@std,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        GaccMat(1,nSubj1+i) = TGM_corr.GroupCount(1)/TGM_total.GroupCount(1);
        GaccMat(2,nSubj1+i) = TGM_corr.GroupCount(2)/TGM_total.GroupCount(2);
        stdMata(:,nSubj1+i) = TGS_total.std_acc;

        % descriptive (cont.)
        rtC1 = T1.rt(and(T1.n0cong==0,T1.acc==1));
        rtI1 = T1.rt(and(T1.n0cong==2,T1.acc==1));
        rtC2 = T2.rt(and(T2.n0cong==0,T2.acc==1));
        rtI2 = T2.rt(and(T2.n0cong==2,T2.acc==1));
        accC1 = T1.acc(T1.n0cong==0);
        accI1 = T1.acc(T1.n0cong==2);
        accC2 = T2.acc(T2.n0cong==0);
        accI2 = T2.acc(T2.n0cong==2);
        accC1(find(accC1==99)) = [];
        accI1(find(accI1==99)) = [];
        accC2(find(accC2==99)) = [];
        accI2(find(accI2==99)) = [];
        for k = 1:drawLH
            if length(rtC1) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtC1))])
                HrtMat1(1,k,nSubj1+i) = mean(rtC1);  % don't need trimmean because outlier trials are marked as 99 above
                HrtvMat1(1,k,nSubj1+i) = var(rtC1);
                HaccMat1(1,k,nSubj1+i) = mean(accC1);
            else
                HrtMat1(1,k,nSubj1+i) = mean(rtC1(1:nTdrawn(k)));
                HrtvMat1(1,k,nSubj1+i) = var(rtC1(1:nTdrawn(k)));
                HaccMat1(1,k,nSubj1+i) = mean(accC1(1:nTdrawn(k)));
            end
            if length(rtI1) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtI1))])
                HrtMat1(2,k,nSubj1+i) = mean(rtI1);
                HrtvMat1(2,k,nSubj1+i) = var(rtI1);
                HaccMat1(2,k,nSubj1+i) = mean(accI1);
            else
                HrtMat1(2,k,nSubj1+i) = mean(rtI1(1:nTdrawn(k)));
                HrtvMat1(2,k,nSubj1+i) = var(rtI1(1:nTdrawn(k)));
                HaccMat1(2,k,nSubj1+i) = mean(accI1(1:nTdrawn(k)));
            end
            if length(rtC2) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtC2))])
                HrtMat2(1,k,nSubj1+i) = mean(rtC2);
                HrtvMat2(1,k,nSubj1+i) = var(rtC2);
                HaccMat2(1,k,nSubj1+i) = mean(accC2);
            else
                HrtMat2(1,k,nSubj1+i) = mean(rtC2(1:nTdrawn(k)));
                HrtvMat2(1,k,nSubj1+i) = var(rtC2(1:nTdrawn(k)));
                HaccMat2(1,k,nSubj1+i) = mean(accC2(1:nTdrawn(k)));
            end
            if length(rtI2) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtI2))])
                HrtMat2(2,k,nSubj1+i) = mean(rtI2);
                HrtvMat2(2,k,nSubj1+i) = var(rtI2);
                HaccMat2(2,k,nSubj1+i) = mean(accI2);
            else
                HrtMat2(2,k,nSubj1+i) = mean(rtI2(1:nTdrawn(k)));
                HrtvMat2(2,k,nSubj1+i) = var(rtI2(1:nTdrawn(k)));
                HaccMat2(2,k,nSubj1+i) = mean(accI2(1:nTdrawn(k)));
            end
        end
    end
end
% Within-subject standard deviation - for simulation
m_ws_std_RT(1,1) = mean(stdMatr(1,:));
m_ws_std_RT(2,1) = mean(stdMatr(2,:));
m_ws_std_Acc(1,1) = mean(stdMata(1,:));
m_ws_std_Acc(2,1) = mean(stdMata(2,:));

%% 2. Simulate data
%% Copulas: Generating dependent multivariate data to simulate the correlation between congruent and incongruent trials
% RT
% 1. Fit model to empirical data set
% [Fi1, xi1] = ecdf(GrtMat(1,:));
% Fi1_sm = ksdensity(GrtMat(1,:),xi1,'function','cdf','width',.01);  % kernel smoothing; adjust width until it matches ecdf
% [Fi2, xi2] = ecdf(GrtMat(2,:));
% Fi2_sm = ksdensity(GrtMat(2,:),xi2,'function','cdf','width',.01);
% figure
% subplot(1,nCond,1)
% stairs(xi1,Fi1,'b','LineWidth',2); hold on
% plot(xi1,Fi1_sm,'r-','LineWidth',1.5)
% xlabel('X1')
% ylabel('Cumulative Probability')
% legend('Empirical','Smoothed','Location','NW'); grid on
% title('Congruent')
% subplot(1,nCond,2)
% stairs(xi2,Fi2,'b','LineWidth',2); hold on
% plot(xi2,Fi2_sm,'r-','LineWidth',1.5)
% xlabel('X1')
% ylabel('Cumulative Probability')
% legend('Empirical','Smoothed','Location','NW'); grid on
% title('Incongruent')
% 2. Calculate correlation
tau_r = corr(GrtMat(1,:)',GrtMat(2,:)','type','kendall');  % Kendall's rank order
nu = 1;  % degree of freedom
rho_r = copulaparam('t',tau_r,nu,'type','kendall');  % corresponding linear correlation parameter for the t copula
% 3. Generate random values from the t copula
U_r = copularnd('t',[1 rho_r; rho_r 1], nu, bs_dist_n);
X1_r = ksdensity(GrtMat(1,:),U_r(:,1),'function','icdf','width',.01);
X2_r = ksdensity(GrtMat(2,:),U_r(:,2),'function','icdf','width',.01);
% 4. Plot and compare w/ empirical
% figure
% scatterhist(GrtMat(1,:),GrtMat(2,:),'Direction','out')
% title('Empirical Mean RT Group Distribution')
% figure
% scatterhist(X1_r,X2_r,'Direction','out')
% title('Simulated')

% Acc
% 1. Fit model to empirical data set
% [Fi1, xi1] = ecdf(GaccMat(1,:));
% Fi1_sm = ksdensity(GaccMat(1,:),xi1,'function','cdf','width',.008);  % kernel smoothing; adjust width until it matches ecdf
% [Fi2, xi2] = ecdf(GaccMat(2,:));
% Fi2_sm = ksdensity(GaccMat(2,:),xi2,'function','cdf','width',.01);
% figure
% subplot(1,2,1)
% stairs(xi1,Fi1,'b','LineWidth',2); hold on
% plot(xi1,Fi1_sm,'r-','LineWidth',1.5)
% xlabel('X1')
% ylabel('Cumulative Probability')
% legend('Empirical','Smoothed','Location','NW'); grid on
% title('Congruent')
% subplot(1,2,2)
% stairs(xi2,Fi2,'b','LineWidth',2); hold on
% plot(xi2,Fi2_sm,'r-','LineWidth',1.5)
% xlabel('X1')
% ylabel('Cumulative Probability')
% legend('Empirical','Smoothed','Location','NW'); grid on
% title('Incongruent')
% 2. Calculate correlation
tau_a = corr(GaccMat(1,:)',GaccMat(2,:)','type','kendall');  % Kendall's rank order
rho_a = copulaparam('t',tau_a,nu,'type','kendall');  % corresponding linear correlation parameter for the t copula
% 3. Generate random values from the t copula
U_a = copularnd('t',[1 rho_a; rho_a 1], nu, bs_dist_n);
X1_a = ksdensity(GaccMat(1,:),U_a(:,1),'function','icdf','width',.008);
X2_a = ksdensity(GaccMat(2,:),U_a(:,2),'function','icdf','width',.01);
% 4. Plot and compare w/ empirical
% figure
% scatterhist(GaccMat(1,:),GaccMat(2,:),'Direction','out')
% title('Empirical Mean Acc Group Distribution')
% figure
% scatterhist(X1_a,X2_a,'Direction','out')
% title('Simulated')

%% Simulate data to get split-half reliability across different sampling sizes
% Preassignment
% Subject X tested in day 1
simCErt_matc = nan(numTest,n_subs,length(nTdrawn));  % mean RT
simCErt_mati = nan(numTest,n_subs,length(nTdrawn));
simCEvar_matc = nan(numTest,n_subs,length(nTdrawn));  % RT variance
simCEvar_mati = nan(numTest,n_subs,length(nTdrawn));
simCEacc_matc = nan(numTest,n_subs,length(nTdrawn));  % accuracy
simCEacc_mati = nan(numTest,n_subs,length(nTdrawn));
vMc = nan(numTest,n_subs,length(nTdrawn));  % drift rate
vMi = nan(numTest,n_subs,length(nTdrawn));
aMc = nan(numTest,n_subs,length(nTdrawn));  % boundary separation
aMi = nan(numTest,n_subs,length(nTdrawn));
TerMc = nan(numTest,n_subs,length(nTdrawn));  % nondecision time
TerMi = nan(numTest,n_subs,length(nTdrawn));

% (same subject) tested in day 2
simCErt_matc2 = nan(numTest,n_subs,length(nTdrawn));
simCErt_mati2 = nan(numTest,n_subs,length(nTdrawn));
simCEvar_matc2 = nan(numTest,n_subs,length(nTdrawn));
simCEvar_mati2 = nan(numTest,n_subs,length(nTdrawn));
simCEacc_matc2 = nan(numTest,n_subs,length(nTdrawn));
simCEacc_mati2 = nan(numTest,n_subs,length(nTdrawn));
vMc2 = nan(numTest,n_subs,length(nTdrawn));
vMi2 = nan(numTest,n_subs,length(nTdrawn));
aMc2 = nan(numTest,n_subs,length(nTdrawn));
aMi2 = nan(numTest,n_subs,length(nTdrawn));
TerMc2 = nan(numTest,n_subs,length(nTdrawn));
TerMi2 = nan(numTest,n_subs,length(nTdrawn));

% Between-subject joint distribution
bs_dist_RT_c = X1_r;  % RT con
bs_dist_RT_i = X2_r;  % inc
bs_dist_Acc_c = X1_a;  % accuracy
bs_dist_Acc_i = X2_a;
bs_dist_Acc_c(find(bs_dist_Acc_c>1)) = [];
bs_dist_Acc_i(find(bs_dist_Acc_i>1)) = [];
minI = min(length(bs_dist_Acc_c),length(bs_dist_Acc_i));

for k = 1:length(nTdrawn)
    for l = 1:numTest
        for n = 1:n_subs
            rankr = randi(bs_dist_n,1,1);  % rank order in RT joint distribution
            ranka = randi(minI,1,1);  % rank order in accuracy joint distribution
            ws_mean_RT_c = bs_dist_RT_c(rankr);  % true mean
            ws_mean_RT_i = bs_dist_RT_i(rankr);
            ws_mean_Acc_c = bs_dist_Acc_c(ranka);
            ws_mean_Acc_i = bs_dist_Acc_i(ranka);

            % Simulated data - within-subject distribution in day 1
            ws_dist_RT_c = randn(1,ws_dist_n)*m_ws_std_RT(1)+ws_mean_RT_c;
            ws_dist_RT_i = randn(1,ws_dist_n)*m_ws_std_RT(2)+ws_mean_RT_i;
            ws_dist_Acc_c = randn(1,ws_dist_n)*m_ws_std_Acc(1)+ws_mean_Acc_c;
            ws_dist_Acc_i = randn(1,ws_dist_n)*m_ws_std_Acc(2)+ws_mean_Acc_i;
            ws_dist_Acc_c(find(ws_dist_Acc_c>1)) = [];  % accuracy should not exceed 1
            ws_dist_Acc_i(find(ws_dist_Acc_i>1)) = [];
            rtC = [];
            rtI = [];
            accC = [];
            accI = [];
            for j = 1:nTdrawn(k)  % number of trials
                disp(['[Sim' num2str(l) '] Sampling step:' num2str(k) ', Subject#' num2str(n) ', Simulation trial#' num2str(j)])
                rtC = [rtC; ws_dist_RT_c(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];  % add random noise
                rtI = [rtI; ws_dist_RT_i(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];
                accC = [accC; ws_dist_Acc_c(randi(length(ws_dist_Acc_c),1,1))+normrnd(0,nzsiga)];
                accI = [accI; ws_dist_Acc_i(randi(length(ws_dist_Acc_i),1,1))+normrnd(0,nzsiga)];
            end
            %rtC = rmoutliers(rtC,"mean");  % 3SD from mean
            %rtI = rmoutliers(rtI,"mean");
            simCErt_matc(l,n,k) = mean(rtC);
            simCErt_mati(l,n,k) = mean(rtI);
            simCEvar_matc(l,n,k) = var(rtC);
            simCEvar_mati(l,n,k) = var(rtI);
            simCEacc_matc(l,n,k) = mean(accC);
            simCEacc_mati(l,n,k) = mean(accI);

            % EZ-diffusion modeling
            [v, a, Ter] = ezdiffusion(simCEacc_matc(l,n,k),simCEvar_matc(l,n,k),simCErt_matc(l,n,k),nTdrawn(k));
            vMc(l,n,k) = v;
            aMc(l,n,k) = a;
            TerMc(l,n,k) = Ter; clear v a Ter
            [v, a, Ter] = ezdiffusion(simCEacc_mati(l,n,k),simCEvar_mati(l,n,k),simCErt_mati(l,n,k),nTdrawn(k));
            vMi(l,n,k) = v;
            aMi(l,n,k) = a;
            TerMi(l,n,k) = Ter; clear v a Ter

            % Simulated data within-subject distribution in day 2
            ws_dist_RT_c = randn(1,ws_dist_n)*m_ws_std_RT(1)+ws_mean_RT_c;
            ws_dist_RT_i = randn(1,ws_dist_n)*m_ws_std_RT(2)+ws_mean_RT_i;
            ws_dist_Acc_c = randn(1,ws_dist_n)*m_ws_std_Acc(1)+ws_mean_Acc_c;
            ws_dist_Acc_i = randn(1,ws_dist_n)*m_ws_std_Acc(2)+ws_mean_Acc_i;
            ws_dist_Acc_c(find(ws_dist_Acc_c>1)) = [];
            ws_dist_Acc_i(find(ws_dist_Acc_i>1)) = [];
            rtC = [];
            rtI = [];
            accC = [];
            accI = [];
            for j = 1:nTdrawn(k) % number of trials
                rtC = [rtC; ws_dist_RT_c(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];
                rtI = [rtI; ws_dist_RT_i(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];
                accC = [accC; ws_dist_Acc_c(randi(length(ws_dist_Acc_c),1,1))+normrnd(0,nzsiga)];
                accI = [accI; ws_dist_Acc_i(randi(length(ws_dist_Acc_i),1,1))+normrnd(0,nzsiga)];
            end
            %rtC = rmoutliers(rtC,"mean");
            %rtI = rmoutliers(rtI,"mean");
            simCErt_matc2(l,n,k) = mean(rtC);
            simCErt_mati2(l,n,k) = mean(rtI);
            simCEvar_matc2(l,n,k) = var(rtC);
            simCEvar_mati2(l,n,k) = var(rtI);
            simCEacc_matc2(l,n,k) = mean(accC);
            simCEacc_mati2(l,n,k) = mean(accI);
            [v, a, Ter] = ezdiffusion(simCEacc_matc2(l,n,k),simCEvar_matc2(l,n,k),simCErt_matc2(l,n,k),nTdrawn(k));
            vMc2(l,n,k) = v;
            aMc2(l,n,k) = a;
            TerMc2(l,n,k) = Ter; clear v a Ter
            [v, a, Ter] = ezdiffusion(simCEacc_mati2(l,n,k),simCEvar_mati2(l,n,k),simCErt_mati2(l,n,k),nTdrawn(k));
            vMi2(l,n,k) = v;
            aMi2(l,n,k) = a;
            TerMi2(l,n,k) = Ter; clear v a Ter
        end
    end
end
% Save simulation data
save('EZDsim_testretestResult_avg','simCErt_matc','simCErt_mati','simCEvar_matc','simCEvar_mati',...
    'simCEacc_matc','simCEacc_mati','vMc','vMi','aMc','aMi','TerMc','TerMi',...
    'simCErt_matc2','simCErt_mati2','simCEvar_matc2','simCEvar_mati2',...
    'simCEacc_matc2','simCEacc_mati2','vMc2','vMi2','aMc2','aMi2','TerMc2','TerMi2')
%load EZDsim_testretestResult_avg

%% Plot results
% Colormap
cmap1 = parula(3);
cmap2 = spring(3);
cmap3 = cool(3);
% difference scores
ceDiffe1 = squeeze(HrtMat1(2,:,:)-HrtMat1(1,:,:));  % RT CE empirical data; session 1
ceDiffe2 = squeeze(HrtMat2(2,:,:)-HrtMat2(1,:,:));  % session 2
ceDiffeA1 = squeeze(HaccMat1(1,:,:)-HaccMat1(2,:,:));  % accuracy
ceDiffeA2 = squeeze(HaccMat2(1,:,:)-HaccMat2(2,:,:));
ceDiff1 = squeeze(mean(simCErt_mati-simCErt_matc,1));  % RT CE simulated data; session 1
ceDiff2 = squeeze(mean(simCErt_mati2-simCErt_matc2,1));  % session 2
ceDiffA1 = squeeze(mean(simCEacc_matc-simCEacc_mati,1));  % accuracy
ceDiffA2 = squeeze(mean(simCEacc_matc2-simCEacc_mati2,1));
vDiff1 = squeeze(mean(vMc-vMi,1));  % drift rate; con-inc
vDiff2 = squeeze(mean(vMc2-vMi2,1));
aDiff1 = squeeze(mean(aMc-aMi,1));  % boundary separation
aDiff2 = squeeze(mean(aMc2-aMi2,1));
TerDiff1 = squeeze(mean(TerMi-TerMc,1));  % nondecision time; inc-con
TerDiff2 = squeeze(mean(TerMi2-TerMc2,1));

% Correlate day 1 and day 2 data
% ICC
figure
for k = 1:length(nTdrawn)
    % empirical data
    if k < 5
        subplot(length(nTdrawn),2,2*(k-1)+1)
        [Icc1] = ICC([ceDiffe1(k,:)',ceDiffe2(k,:)'],'A-k',0.05);
        scatter(ceDiffe1(k,:),ceDiffe2(k,:),8,'filled')
        %lsline
        str = sprintf(' ICC = %1.2f',Icc1);
        rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
        set(rText, 'fontsize', 7, 'verticalalignment', 'top', 'horizontalalignment','left');
        refline
        axis square
        ylabel(['N = ' num2str(2*nTdrawn(k))])
    end
    % simulated data
    subplot(length(nTdrawn),2,2*(k-1)+2)
    [Icc2] = ICC([ceDiff1(:,k),ceDiff2(:,k)],'A-k',0.05);
    scatter(ceDiff1(:,k),ceDiff2(:,k),8,'filled')
    %lsline
    str = sprintf(' ICC = %1.2f',Icc2);
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText, 'fontsize', 7, 'verticalalignment', 'top', 'horizontalalignment','left');
    refline
    axis square
    ylabel(['N = ' num2str(2*nTdrawn(k))])
end
subplot(length(nTdrawn),2,1)
title('Hedge et al. (2018)')
subplot(length(nTdrawn),2,2)
title('Simulated data')
sgtitle('CE - IntraClass Correlation (ICC)')

% DDM parameters (difference scores)
for k = 1:length(nTdrawn)
    % 1) RT CE
    figure
    scatter(ceDiff1(:,k),ceDiff2(:,k),120,cmap2(2,:),'filled')
    xlim([min([ceDiff1(:);ceDiff2(:)]) max([ceDiff1(:);ceDiff2(:)])])
    ylim([min([ceDiff1(:);ceDiff2(:)]) max([ceDiff1(:);ceDiff2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [ICc] = ICC([ceDiff1(:,k),ceDiff2(:,k)],'A-k',0.05);
    str = sprintf(' ICC = %1.2f',ICc); clear ICc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'fontweight','bold','fontsize',40,'verticalalignment','top','horizontalalignment','left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontWeight','bold','FontSize',9.5)

    % 2) Acc CE
    figure
    scatter(ceDiffA1(:,k),ceDiffA2(:,k),120,cmap2(2,:),'filled')
    xlim([min([ceDiffA1(:);ceDiffA2(:)]) max([ceDiffA1(:);ceDiffA2(:)])])
    ylim([min([ceDiffA1(:);ceDiffA2(:)]) max([ceDiffA1(:);ceDiffA2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [ICc] = ICC([ceDiffA1(:,k),ceDiffA2(:,k)],'A-k',0.05);
    str = sprintf(' ICC = %1.2f',ICc); clear ICc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'fontweight','bold','fontsize',40,'verticalalignment','top','horizontalalignment','left');
    refline
    axis square

    % 3) v
    figure
    scatter(vDiff1(:,k),vDiff2(:,k),120,cmap1(1,:),'filled')
    xlim([min([vDiff1(:);vDiff2(:)]) max([vDiff1(:);vDiff2(:)])])
    ylim([min([vDiff1(:);vDiff2(:)]) max([vDiff1(:);vDiff2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [ICc] = ICC([vDiff1(:,k),vDiff2(:,k)],'A-k',0.05);
    str = sprintf(' ICC = %1.2f',ICc); clear ICc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'fontweight','bold','fontsize',40,'verticalalignment','top','horizontalalignment','left');
    refline
    axis square

    % 4) a
    figure
    scatter(aDiff1(:,k),aDiff2(:,k),120,cmap1(2,:),'filled')
    xlim([real(min([aDiff1(:);aDiff2(:)])) real(max([aDiff1(:);aDiff2(:)]))])
    ylim([real(min([aDiff1(:);aDiff2(:)])) real(max([aDiff1(:);aDiff2(:)]))])
    %yticks([0.003 0.005 0.008 0.01])
    ax = gca;
    ax.FontSize = 36;
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    [ICc] = ICC([aDiff1(:,k),aDiff2(:,k)],'A-k',0.05);
    str = sprintf(' ICC = %1.2f',ICc); clear ICc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'fontweight','bold','fontsize',40,'verticalalignment','top','horizontalalignment','left');
    refline
    axis square

    % 5) Ter
    figure
    scatter(TerDiff1(:,k),TerDiff2(:,k),120,cmap1(3,:),'filled')
    xlim([min([TerDiff1(:);TerDiff2(:)]) max([TerDiff1(:);TerDiff2(:)])])
    ylim([min([TerDiff1(:);TerDiff2(:)]) max([TerDiff1(:);TerDiff2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [ICc] = ICC([TerDiff1(:,k),TerDiff2(:,k)],'A-k',0.05);
    str = sprintf(' ICC = %1.2f',ICc); clear ICc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'fontweight','bold','fontsize',40,'verticalalignment','top','horizontalalignment','left');
    refline
    axis square
end