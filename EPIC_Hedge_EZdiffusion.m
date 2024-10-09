%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
%
% To run this sciprt:
% You need Hedge et al.'s (2018) flanker task data
% (download at https://osf.io/cwzds/),
% a Matlab function file named "ezdiffusion.m" (download at
% https://osf.io/jk9nb), and ICC.m (download at
% https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
% What this script does:
% Plots the EZ-diffusion model parameters across varying numbers of trials
% (RT needs to be in SECS to run EZ-diffusion modeling)
%
% What this script outputs:
% Supp. Figs. 19A and 19C
%
% Created on 05/16/2023 by HJ Lee
% Last modified on 01/18/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')

%% Parameter settings
nCond = 2;  % number of experimental conditions: congruent, incongruent
nStd = 3;  % standard deviation from mean for excluding outliers
UpB = 97.5;  % 95% confidence interval
LwB = 2.5;
nBt = 1000;  % number of bootstrapping

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
nSession = 2;
nBlocks = 5;
nTrials = 144;  % per block
%nTrialsCE = 96;  % number of trials for calculating CE (excluding neutral)

nTdrawn = [50, 100, 200, 400, 800, 1600 3200]/2;  % number of trials drawn (per condition)
drawLH = 4;  % index for the max number of trials drawn for Hedge et al.'s data; max 480 trials (cf. nTdrawn)

%% 1. Load and preprocess Hedge et al.'s (2018) flanker task data
% Descriptive stats
HrtMat1 = nan(nCond,drawLH,nSubjT);  % mean RT (first half) in secs; session 1/first half
HrtMat2 = nan(nCond,drawLH,nSubjT);  % session 2/second half
HrtvMat1 = nan(nCond,drawLH,nSubjT);  % RT variance
HrtvMat2 = nan(nCond,drawLH,nSubjT);
HaccMat1 = nan(nCond,drawLH,nSubjT);  % mean accuracy
HaccMat2 = nan(nCond,drawLH,nSubjT);
% Bootstrapping
Hrt_bt1 = cell(nSubjT,1);  % mean RT (1st half) bootstrap for confidence interval
Hrt_bt2 = cell(nSubjT,1);
Hrtv_bt1 = cell(nSubjT,1);
Hrtv_bt2 = cell(nSubjT,1);
Hacc_bt1 = cell(nSubjT,1);
Hacc_bt2 = cell(nSubjT,1);
% Study 1 data
for i = 1:nSubj1
    for t = 1:nTask  % 1
        % preassignment for bootstrapping
        rt_bt1 = nan(nCond,drawLH,nBt);
        rt_bt2 = nan(nCond,drawLH,nBt);
        rtv_bt1 = nan(nCond,drawLH,nBt);
        rtv_bt2 = nan(nCond,drawLH,nBt);
        acc_bt1 = nan(nCond,drawLH,nBt);
        acc_bt2 = nan(nCond,drawLH,nBt);

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
        %T = [T1;T2];  % combine session 1&2
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
                HrtMat1(1,k,i) = mean(rtC1);  % don't need trimmean because outlier trials were marked as 99 and removed
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
            simrtC1 = [];
            simrtC2 = [];
            simrtI1 = [];
            simrtI2 = [];
            simrtvC1 = [];
            simrtvC2 = [];
            simrtvI1 = [];
            simrtvI2 = [];
            simaccC1 = [];
            simaccC2 = [];
            simaccI1 = [];
            simaccI2 = [];
            for i_b = 1:nBt
                tc1 = datasample(rtC1,nTdrawn(k));
                tc2 = datasample(rtC2,nTdrawn(k));
                ti1 = datasample(rtI1,nTdrawn(k));
                ti2 = datasample(rtI2,nTdrawn(k));
                simrtC1 = [simrtC1, mean(tc1)];
                simrtC2 = [simrtC2, mean(tc2)];
                simrtI1 = [simrtI1, mean(ti1)];
                simrtI2 = [simrtI2, mean(ti2)];
                simrtvC1 = [simrtvC1, var(tc1)];
                simrtvC2 = [simrtvC2, var(tc2)];
                simrtvI1 = [simrtvI1, var(ti1)];
                simrtvI2 = [simrtvI2, var(ti2)];
                simaccC1 = [simaccC1, mean(datasample(accC1,nTdrawn(k)))];
                simaccC2 = [simaccC2, mean(datasample(accC2,nTdrawn(k)))];
                simaccI1 = [simaccI1, mean(datasample(accI1,nTdrawn(k)))];
                simaccI2 = [simaccI2, mean(datasample(accI2,nTdrawn(k)))];
            end
            rt_bt1(1,k,:) = simrtC1;
            rt_bt2(1,k,:) = simrtC2;
            rt_bt1(2,k,:) = simrtI1;
            rt_bt2(2,k,:) = simrtI2;
            rtv_bt1(1,k,:) = simrtvC1;
            rtv_bt2(1,k,:) = simrtvC2;
            rtv_bt1(2,k,:) = simrtvI1;
            rtv_bt2(2,k,:) = simrtvI2;
            acc_bt1(1,k,:) = simaccC1;
            acc_bt2(1,k,:) = simaccC2;
            acc_bt1(2,k,:) = simaccI1;
            acc_bt2(2,k,:) = simaccI2;
        end
        Hrt_bt1{i,1} = rt_bt1;
        Hrt_bt2{i,1} = rt_bt2;
        Hrtv_bt1{i,1} = rtv_bt1;
        Hrtv_bt2{i,1} = rtv_bt2;
        Hacc_bt1{i,1} = acc_bt1;
        Hacc_bt2{i,1} = acc_bt2;
    end
end
% Study 2 data
for i = 1:nSubj2
    for t = 1:nTask
        % preassignment
        rt_bt1 = nan(nCond,drawLH,nBt);
        rt_bt2 = nan(nCond,drawLH,nBt);
        rtv_bt1 = nan(nCond,drawLH,nBt);
        rtv_bt2 = nan(nCond,drawLH,nBt);
        acc_bt1 = nan(nCond,drawLH,nBt);
        acc_bt2 = nan(nCond,drawLH,nBt);

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
        %T = [T1;T2];
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
                HrtMat1(1,k,nSubj1+i) = mean(rtC1);  % don't need trimmean because outlier trials are marked as 99 above and removed
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
            simrtC1 = [];
            simrtC2 = [];
            simrtI1 = [];
            simrtI2 = [];
            simrtvC1 = [];
            simrtvC2 = [];
            simrtvI1 = [];
            simrtvI2 = [];
            simaccC1 = [];
            simaccC2 = [];
            simaccI1 = [];
            simaccI2 = [];
            for i_b = 1:nBt
                tc1 = datasample(rtC1,nTdrawn(k));
                tc2 = datasample(rtC2,nTdrawn(k));
                ti1 = datasample(rtI1,nTdrawn(k));
                ti2 = datasample(rtI2,nTdrawn(k));
                simrtC1 = [simrtC1, mean(tc1)];
                simrtC2 = [simrtC2, mean(tc2)];
                simrtI1 = [simrtI1, mean(ti1)];
                simrtI2 = [simrtI2, mean(ti2)];
                simrtvC1 = [simrtvC1, var(tc1)];
                simrtvC2 = [simrtvC2, var(tc2)];
                simrtvI1 = [simrtvI1, var(ti1)];
                simrtvI2 = [simrtvI2, var(ti2)];
                simaccC1 = [simaccC1, mean(datasample(accC1,nTdrawn(k)))];
                simaccC2 = [simaccC2, mean(datasample(accC2,nTdrawn(k)))];
                simaccI1 = [simaccI1, mean(datasample(accI1,nTdrawn(k)))];
                simaccI2 = [simaccI2, mean(datasample(accI2,nTdrawn(k)))];
            end
            rt_bt1(1,k,:) = simrtC1;
            rt_bt2(1,k,:) = simrtC2;
            rt_bt1(2,k,:) = simrtI1;
            rt_bt2(2,k,:) = simrtI2;
            rtv_bt1(1,k,:) = simrtvC1;
            rtv_bt2(1,k,:) = simrtvC2;
            rtv_bt1(2,k,:) = simrtvI1;
            rtv_bt2(2,k,:) = simrtvI2;
            acc_bt1(1,k,:) = simaccC1;
            acc_bt2(1,k,:) = simaccC2;
            acc_bt1(2,k,:) = simaccI1;
            acc_bt2(2,k,:) = simaccI2;
        end
        Hrt_bt1{nSubj1+i,1} = rt_bt1;
        Hrt_bt2{nSubj1+i,1} = rt_bt2;
        Hrtv_bt1{nSubj1+i,1} = rtv_bt1;
        Hrtv_bt2{nSubj1+i,1} = rtv_bt2;
        Hacc_bt1{nSubj1+i,1} = acc_bt1;
        Hacc_bt2{nSubj1+i,1} = acc_bt2;
    end
end

%% 2. EZ-diffusion modeling
% EZ-diffusion model parameters
vMatH1 = nan(nCond,drawLH,nSubjT);  % drift rate; session 1
aMatH1 = nan(nCond,drawLH,nSubjT);  % boundary separation
TerMatH1 = nan(nCond,drawLH,nSubjT);  % nondecision time
vMatH2 = nan(nCond,drawLH,nSubjT);  % session 2
aMatH2 = nan(nCond,drawLH,nSubjT);
TerMatH2 = nan(nCond,drawLH,nSubjT);
% Standard error (95% confidence interval of the mean
vMatHse = nan(nSession,drawLH,nSubjT);  % nSession=2, sessions 1&2
aMatHse = nan(nSession,drawLH,nSubjT);
TerMatHse = nan(nSession,drawLH,nSubjT);
% ICC
vICC1 = nan(drawLH,nBt,nSubjT);  % for calculating 95% confidence interval of ICC
vICC2 = nan(drawLH,nBt,nSubjT);
aICC1 = nan(drawLH,nBt,nSubjT);
aICC2 = nan(drawLH,nBt,nSubjT);
TerICC1 = nan(drawLH,nBt,nSubjT);
TerICC2 = nan(drawLH,nBt,nSubjT);
for i = 1:nSubjT
    % Preassignment for bootstrapping
    RT_bt1 = Hrt_bt1{i,1};
    RT_bt2 = Hrt_bt2{i,1};
    RTV_bt1 = Hrtv_bt1{i,1};
    RTV_bt2 = Hrtv_bt2{i,1};
    ACC_bt1 = Hacc_bt1{i,1};
    ACC_bt2 = Hacc_bt2{i,1};

    vSE1 = nan(nCond,drawLH,nBt);  % 95% confidence interval
    vSE2 = nan(nCond,drawLH,nBt);
    aSE1 = nan(nCond,drawLH,nBt);
    aSE2 = nan(nCond,drawLH,nBt);
    TerSE1 = nan(nCond,drawLH,nBt);
    TerSE2 = nan(nCond,drawLH,nBt);

    for k = 1:drawLH
        [v1, a1, Ter1] = ezdiffusion(HaccMat1(1,k,i),HrtvMat1(1,k,i),HrtMat1(1,k,i),nTdrawn(k));  % con
        [v2, a2, Ter2] = ezdiffusion(HaccMat1(2,k,i),HrtvMat1(2,k,i),HrtMat1(2,k,i),nTdrawn(k));  % inc
        vMatH1(1,k,i) = v1;  % con
        vMatH1(2,k,i) = v2;  % inc
        aMatH1(1,k,i) = a1;
        aMatH1(2,k,i) = a2;
        TerMatH1(1,k,i) = Ter1;
        TerMatH1(2,k,i) = Ter2;
        clear v1 v2 a1 a2 Ter1 Ter2
        [v1, a1, Ter1] = ezdiffusion(HaccMat2(1,k,i),HrtvMat2(1,k,i),HrtMat2(1,k,i),nTdrawn(k));
        [v2, a2, Ter2] = ezdiffusion(HaccMat2(2,k,i),HrtvMat2(2,k,i),HrtMat2(2,k,i),nTdrawn(k));
        vMatH2(1,k,i) = v1;
        vMatH2(2,k,i) = v2;
        aMatH2(1,k,i) = a1;
        aMatH2(2,k,i) = a2;
        TerMatH2(1,k,i) = Ter1;
        TerMatH2(2,k,i) = Ter2;
        clear v1 v2 a1 a2 Ter1 Ter2

        % EZ-diffusion modeling of the bootstrapped data
        for ib = 1:nBt
            [vSE1(1,k,ib),aSE1(1,k,ib),TerSE1(1,k,ib)] = ezdiffusion(ACC_bt1(1,k,ib),RTV_bt1(1,k,ib),RT_bt1(1,k,ib),nTdrawn(k));
            [vSE2(1,k,ib),aSE2(1,k,ib),TerSE2(1,k,ib)] = ezdiffusion(ACC_bt2(1,k,ib),RTV_bt2(1,k,ib),RT_bt2(1,k,ib),nTdrawn(k));
            [vSE1(2,k,ib),aSE1(2,k,ib),TerSE1(2,k,ib)] = ezdiffusion(ACC_bt1(2,k,ib),RTV_bt1(2,k,ib),RT_bt1(2,k,ib),nTdrawn(k));
            [vSE2(2,k,ib),aSE2(2,k,ib),TerSE2(2,k,ib)] = ezdiffusion(ACC_bt2(2,k,ib),RTV_bt2(2,k,ib),RT_bt2(2,k,ib),nTdrawn(k));
        end
        % Standard error
        VceSE1 = squeeze(vSE1(1,k,:)-vSE1(2,k,:));  % con-inc
        vICC1(k,:,i) = VceSE1;  % to calculate 95% confidence interval of the ICC
        vMatHse(1,k,i) = (prctile(VceSE1,UpB)-prctile(VceSE1,LwB))/2;
        VceSE2 = squeeze(vSE2(1,k,:)-vSE2(2,k,:));
        vICC2(k,:,i) = VceSE2;
        vMatHse(2,k,i) = (prctile(VceSE2,UpB)-prctile(VceSE2,LwB))/2;
        AceSE1 = squeeze(aSE1(1,k,:)-aSE1(2,k,:));
        aICC1(k,:,i) = AceSE1;
        aMatHse(1,k,i) = (prctile(AceSE1,UpB)-prctile(AceSE1,LwB))/2;
        AceSE2 = squeeze(aSE2(1,k,:)-aSE2(2,k,:));
        aICC2(k,:,i) = AceSE2;
        aMatHse(2,k,i) = (prctile(AceSE2,UpB)-prctile(AceSE2,LwB))/2;
        TerceSE1 = squeeze(TerSE1(2,k,:)-TerSE1(1,k,:));
        TerICC1(k,:,i) = TerceSE1;
        TerMatHse(1,k,i) = (prctile(TerceSE1,UpB)-prctile(TerceSE1,LwB))/2;
        TerceSE2 = squeeze(TerSE2(2,k,:)-TerSE2(1,k,:));
        TerICC2(k,:,i) = TerceSE2;
        TerMatHse(2,k,i) = (prctile(TerceSE2,UpB)-prctile(TerceSE2,LwB))/2;
    end
end

%% Plot the results: EZD modeling parameters
cmap1 = parula(3);
cmap2 = spring(3);
cmap3 = cool(3);
ceTmp1 = squeeze(HrtMat1(2,:,:)-HrtMat1(1,:,:));  % inc-con
ceTmp2 = squeeze(HrtMat2(2,:,:)-HrtMat2(1,:,:));
ceTmpa1 = squeeze(HaccMat1(1,:,:)-HaccMat1(2,:,:));  % con-inc
ceTmpa2 = squeeze(HaccMat2(1,:,:)-HaccMat2(2,:,:));
vTmp1 = squeeze(vMatH1(1,:,:)-vMatH1(2,:,:));  % con-inc
vTmp2 = squeeze(vMatH2(1,:,:)-vMatH2(2,:,:));
aTmp1 = squeeze(aMatH1(1,:,:)-aMatH1(2,:,:));
aTmp2 = squeeze(aMatH2(1,:,:)-aMatH2(2,:,:));
TerTmp1  = squeeze(TerMatH1(2,:,:)-TerMatH1(1,:,:));  % inc-con
TerTmp2 = squeeze(TerMatH2(2,:,:)-TerMatH2(1,:,:));

% x,y lim
xL = nan(drawLH,2,3);  % 2:min,max; 3:v,a,Ter
yL = nan(drawLH,2,3);
for k = 1:drawLH
    xL(k,1,1) = min(vTmp1(k,:)-squeeze(vMatHse(1,k,:))');
    xL(k,2,1) = max(vTmp1(k,:)+squeeze(vMatHse(1,k,:))');
    xL(k,1,2) = min(aTmp1(k,:)-squeeze(aMatHse(1,k,:))');
    xL(k,2,2) = max(aTmp1(k,:)+squeeze(aMatHse(1,k,:))');
    xL(k,1,3) = min(TerTmp1(k,:)-squeeze(TerMatHse(1,k,:))');
    xL(k,2,3) = max(TerTmp1(k,:)+squeeze(TerMatHse(1,k,:))');
    yL(k,1,1) = min(vTmp2(k,:)-squeeze(vMatHse(2,k,:))');
    yL(k,2,1) = max(vTmp2(k,:)+squeeze(vMatHse(2,k,:))');
    yL(k,1,2) = min(aTmp2(k,:)-squeeze(aMatHse(2,k,:))');
    yL(k,2,2) = max(aTmp2(k,:)+squeeze(aMatHse(2,k,:))');
    yL(k,1,3) = min(TerTmp2(k,:)-squeeze(TerMatHse(2,k,:))');
    yL(k,2,3) = max(TerTmp2(k,:)+squeeze(TerMatHse(2,k,:))');
end
vminX = min(xL(:,1,1));
vmaxX = max(xL(:,2,1));
vminY = min(yL(:,1,1));
vmaxY = max(yL(:,2,1));
aminX = min(xL(:,1,2));
amaxX = max(xL(:,2,2));
aminY = min(yL(:,1,2));
amaxY = max(yL(:,2,2));
TerminX = min(xL(:,1,3));
TermaxX = max(xL(:,2,3));
TerminY = min(yL(:,1,3));
TermaxY = max(yL(:,2,3));

% WITHOUT confidence intervals for the modeling parameters
for k = 1:drawLH
    % 1) RT CE
    figure
    scatter(ceTmp1(k,:),ceTmp2(k,:),60,[0 0.4 0.9],'filled')
    xlim([min([ceTmp1(:);ceTmp2(:)]) max([ceTmp1(:);ceTmp2(:)])])
    ylim([min([ceTmp1(:);ceTmp2(:)]) max([ceTmp1(:);ceTmp2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [Icc] = ICC([ceTmp1(k,:)',ceTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc);
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'FontSize',40,'FontWeight','bold','verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontSize',35,'FontWeight','bold')
    %title('Reaction time','FontSize',30)

    % 2) Acc CE
    figure
    scatter(ceTmpa1(k,:),ceTmpa2(k,:),60,[0 0.2 1],'filled')
    xlim([min([ceTmpa1(:);ceTmpa2(:)]) max([ceTmpa1(:);ceTmpa2(:)])])
    ylim([min([ceTmpa1(:);ceTmpa2(:)]) max([ceTmpa1(:);ceTmpa2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [Icc] = ICC([ceTmpa1(k,:)',ceTmpa2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'FontSize', 40,'FontWeight','bold','verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontSize',35,'FontWeight','bold')
    %title('Accuracy','FontSize',30)

    % 3) v
    figure
    scatter(vTmp1(k,:),vTmp2(k,:),60,cmap1(1,:),'filled')
    xlim([min([vTmp1(:);vTmp2(:)]) max([vTmp1(:);vTmp2(:)])])
    ylim([min([vTmp1(:);vTmp2(:)]) max([vTmp1(:);vTmp2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [Icc] = ICC([vTmp1(k,:)',vTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'FontSize',40,'FontWeight','bold','verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontSize',35,'FontWeight','bold')
    %title('Drift rate','FontSize',30)

    % 4) a
    figure
    scatter(aTmp1(k,:),aTmp2(k,:),60,cmap1(2,:),'filled')
    xlim([min([aTmp1(:);aTmp2(:)]) max([aTmp1(:);aTmp2(:)])])
    ylim([min([aTmp1(:);aTmp2(:)]) max([aTmp1(:);aTmp2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [Icc] = ICC([aTmp1(k,:)',aTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'FontSize',40,'FontWeight','bold','verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontSize',35,'FontWeight','bold')
    %title('Boundary separation','FontSize',30)

    % 5) Ter
    figure
    scatter(TerTmp1(k,:),TerTmp2(k,:),60,cmap1(3,:),'filled')
    xlim([min([TerTmp1(:);TerTmp2(:)]) max([TerTmp1(:);TerTmp2(:)])])
    ylim([min([TerTmp1(:);TerTmp2(:)]) max([TerTmp1(:);TerTmp2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [Icc] = ICC([TerTmp1(k,:)',TerTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'FontSize',40,'FontWeight','bold','verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontSize',35,'FontWeight','bold')
    %title('Nondecision time','FontSize',30)
end

%% Plot ICC and its 95% confidence interval as a function of number of trials
myICC = nan(drawLH,3);  % 3: EZ-diffusion model parameters
ICC_ci = nan(drawLH,3);
for k = 1:drawLH
    myICC(k,1) = ICC([vTmp1(k,:)',vTmp2(k,:)'],'A-k',0.05);  % empirical data
    myICC(k,2) = ICC([aTmp1(k,:)',aTmp2(k,:)'],'A-k',0.05);
    myICC(k,3) = ICC([TerTmp1(k,:)',TerTmp2(k,:)'],'A-k',0.05);
    tmp1 = nan(drawLH,nBt);  % drift rate
    tmp2 = nan(drawLH,nBt);  % boundary separation
    tmp3 = nan(drawLH,nBt);  % nondecision time
    for i_b = 1:nBt
        tmp1(k,i_b) = ICC([squeeze(vICC1(k,i_b,:)),squeeze(vICC2(k,i_b,:))],'A-k',0.05);
        tmp2(k,i_b) = ICC([squeeze(aICC1(k,i_b,:)),squeeze(aICC2(k,i_b,:))],'A-k',0.05);
        tmp3(k,i_b) = ICC([squeeze(TerICC1(k,i_b,:)),squeeze(TerICC2(k,i_b,:))],'A-k',0.05);
    end
    ICC_ci(k,1) = (prctile(tmp1(k,:),UpB)-prctile(tmp1(k,:),LwB))/2;
    ICC_ci(k,2) = (prctile(tmp2(k,:),UpB)-prctile(tmp2(k,:),LwB))/2;
    ICC_ci(k,3) = (prctile(tmp3(k,:),UpB)-prctile(tmp3(k,:),LwB))/2;
end
% 95% confidence interval
x = 1:drawLH;
xconf = [x x(end:-1:1)];
%figure
for pa = 1:3  % 3 EZD parameters
    %subplot(3,1,pa)
    figure
    CIl = myICC(:,pa)-ICC_ci(:,pa);
    CIu = myICC(:,pa)+ICC_ci(:,pa);
    u = CIu;
    yconf = [CIl u(end:-1:1)];
    pl = fill(xconf,yconf(:),'b');
    %pl.FaceColor = [1 0.5 0];
    pl.FaceColor = cmap1(pa,:);
    pl.EdgeColor = 'none';
    hold on
    %p2 = plot(x,myICC(:,pa),'Color',cmap1(pa,:),'LineWidth',2.4);
    p2 = plot(x,myICC(:,pa),'Color',[1 0.5 0],'LineWidth',2.4);
    set(gca,'FontSize',18)
    xlabel('Number of trials','FontSize',19)
    xticks(1:drawLH)
    xticklabels([50 100 200 400])
    legend('95% confidence interval','ICC','Location','southeast','FontSize',15)
    if pa == 1
        title('Drift rate','FontSize',20)
    elseif pa == 2
        title('Boundary separation','FontSize',20)
    elseif pa == 3
        title('Nondecision time','FontSize',20)
    end
    ylim([-0.4 1])
end
%sgtitle('Hedge et al.''s (2018) Data','FontSize',15)