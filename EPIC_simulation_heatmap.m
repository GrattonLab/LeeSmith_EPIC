%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
%
% To run this script:
% You need Hedge et al.'s (2018) flanker task data
% (download at https://osf.io/cwzds/)
% You also need ICC.m (download at
% https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
% What this script does:
% Rank order consistency as a function of sample size and number of trials
%
% What this script outputs:
% Figure. 7: Heatmap of the ICC for the rank orders between the true mean
% and the apparent mean
%
% Created on 05/02/2023 by HJ Lee
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
nTask = 1;  % just use the flanker task data; otherwise 2
taskStrng = {'Flanker','Stroop'};
%nSession = 2;
nBlocks = 5;
nTrials = 144;  % per block
%nTrialsCE = 96;  % number of trials for calculating CE (excluding neutral)

nTdrawn = [50, 100, 200, 400, 800, 1600, 3200]/nCond;  % key manipulation
%n_subs = 100;  % will manipulate sample size below for simulation 1
bs_dist_n = 10000;
ws_dist_n = 10000;
numTest = 100;
nzsigr = 0.24;  % noise sigma; determined by evaluating SSE of ICCs
nzsiga = 0.22;

%% 1. Load and preprocess Hedge et al.'s (2018) flanker task data
GrtMat = nan(nCond,nSubjT);  % grand mean RT; need this for simulation (copulas)
GaccMat = nan(nCond,nSubjT);  % accuracy
stdMatr = nan(nCond,nSubjT);  % within-subject STD; need this also for simulation
stdMata = nan(nCond,nSubjT);
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
        % T
        acc_1_ids = find(T.acc==1);
        if length(find(T.acc==1))/length(table2array(T))*100 < 70
            disp(['Study1 ' taskStrng{t} ' task subject#' num2str(subjID1(i)) ' has below 70% accuracy'])
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
        GrtMat(:,i) = TGMrt.mean_rt;
        stdMatr(:,i) = TGSrt.std_rt;
        % Calculate mean Acc
        TGM_total = varfun(@mean,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGM_corr = varfun(@mean,T(T.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS_total = varfun(@std,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        GaccMat(1,i) = TGM_corr.GroupCount(1)/TGM_total.GroupCount(1);
        GaccMat(2,i) = TGM_corr.GroupCount(2)/TGM_total.GroupCount(2);
        stdMata(:,i) = TGS_total.std_acc;
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
    end
end
% Within-subject standard deviation - for simulation
m_ws_std_RT(1,1) = mean(stdMatr(1,:));
m_ws_std_RT(2,1) = mean(stdMatr(2,:));
m_ws_std_Acc(1,1) = mean(stdMata(1,:));
m_ws_std_Acc(2,1) = mean(stdMata(2,:));

%% 2. Simulate data
%% Copulas: Simulating dependent multivariate data
% RT
tau_r = corr(GrtMat(1,:)',GrtMat(2,:)','type','kendall');  % Kendall's rank order
nu = 1;  % degree of freedom
rho_r = copulaparam('t',tau_r,nu,'type','kendall');  % corresponding linear correlation parameter for the t copula
U_r = copularnd('t',[1 rho_r; rho_r 1], nu, bs_dist_n);
X1_r = ksdensity(GrtMat(1,:),U_r(:,1),'function','icdf','width',.01);
X2_r = ksdensity(GrtMat(2,:),U_r(:,2),'function','icdf','width',.01);
% Acc
tau_a = corr(GaccMat(1,:)',GaccMat(2,:)','type','kendall');  % Kendall's rank order
rho_a = copulaparam('t',tau_a,nu,'type','kendall');  % corresponding linear correlation parameter for the t copula
U_a = copularnd('t',[1 rho_a; rho_a 1], nu, bs_dist_n);
X1_a = ksdensity(GaccMat(1,:),U_a(:,1),'function','icdf','width',.008);
X2_a = ksdensity(GaccMat(2,:),U_a(:,2),'function','icdf','width',.01);

%% Simulation - Rank order consistency as a function of sample size and numbers of trials per subject (Figure. 6)
% Preassignment
nSubject = [50 100 200 300 400 500 1000 2000 4000];
rnkCCs = nan(length(nTdrawn),length(nSubject),numTest);  % Spearman correlation of the ranks between true mean and apparent mean
rnkCCk = nan(length(nTdrawn),length(nSubject),numTest);  % Kendall's tau
ICCr = nan(length(nTdrawn),length(nSubject),numTest);  % ICC
rnkBN = nan(length(nTdrawn),length(nSubject),numTest);  % 1(same)/0(diff) binary

% Between-subject joint distribution
bs_dist_RT_c = X1_r;  % RT con
bs_dist_RT_i = X2_r;  % inc
for l = 1:numTest
    for k = 1:length(nTdrawn)  % var1: number of trials
        for n = 1:length(nSubject)  % var2: sample size
            % preassignment
            trueM = nan(nSubject(n),1);  % true mean
            appM = nan(nSubject(n),1);
            randid = randperm(nSubject(n));
            for i = 1:nSubject(n)
                disp(['Testing ' num2str(l) ' - Sampling size: ' num2str(k) ', Sample size: ' num2str(n) ' in progress of subject#' num2str(i)])
                ws_mean_RT_c = bs_dist_RT_c(randid(i));
                ws_mean_RT_i = bs_dist_RT_i(randid(i));
                trueM(i,1) = ws_mean_RT_i-ws_mean_RT_c;  % store true mean CE
                % within-subject distribution
                ws_dist_RT_c = randn(1,ws_dist_n)*m_ws_std_RT(1)+ws_mean_RT_c;
                ws_dist_RT_i = randn(1,ws_dist_n)*m_ws_std_RT(2)+ws_mean_RT_i;
                rtC = [];
                rtI = [];
                for j = 1:nTdrawn(k)  % number of trials drawn
                    rtC = [rtC; ws_dist_RT_c(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];
                    rtI = [rtI; ws_dist_RT_i(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];
                end
                %appM(i,1) = trimmean(rtI,15)-trimmean(rtC,15);
                %appM(i,1) = mean(rmoutliers(rtI,"mean"))-mean(rmoutliers(rtC,"mean"))
                appM(i,1) = mean(rtI)-mean(rtC);
            end
            % rank order
            [~,p1] = sort(trueM,'descend');
            r1 = 1:length(trueM);  % nSubject(n)
            r1(p1) = r1;
            [~,p2] = sort(appM,'descend');
            r2 = 1:length(appM);
            r2(p2) = r2;
            % Spearman's correlation
            rnkCCs(k,n,l) = corr(r1',r2','type','Spearman');
            % Kendall's tau
            rnkCCk(k,n,l) = corr(r1',r2','type','kendall');
            % ICC
            [r,LB,UB,F,df1,df2,p] = ICC([r1',r2'],'C-k',0.05);  % 'A-k'
            ICCr(k,n,l) = r;
            % Binary map
            rnkBN(k,n,l) = sum(r1==r2)/length(r1);
        end
    end
end
%save('rankOrderConsistencySIM','rnkCCs','rnkCCk','ICCr','rnkBN')
%load rankOrderConsistencySIM

%% Plot - heatmap
% ICC
figure
%h = heatmap(round(mean(ICCr,3).*100)/100,'ColorLimits',[0 1]);
h = heatmap(round(mean(ICCr,3).*100)/100,'ColorLimits',[0 1],'CellLabelColor','none');
h.XLabel = 'Number of subjects';
h.YLabel = 'Number of trials';
Xlabels = nSubject;
Ylabels = nCond*nTdrawn;
h.FontSize = 15;
h.XDisplayLabels = Xlabels;
h.YDisplayLabels = Ylabels;
%h.Title = 'Correlation coefficients of the true and apparent rank orders';
h.Colormap = spring;