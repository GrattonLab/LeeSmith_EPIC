%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC - Public Data 2
%
%% Instructions for running the script:
% You will need the following.
%   - Hedge et al.'s (2018) flanker task data (https://osf.io/cwzds/)
%   - Mat-file: ICC.m (https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
%% Purpose of the script:
% This script ranks participants based on their congruency effect (CE) and
% incongruent trial performance to assess consistency, particularly with
% simulations of extended trial sampling.
%
%% Outputs:
% Extended Data Figure 2C and 2F
%
% Created on 05/02/2023 by HJ Lee
% Last modified on 02/04/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')

%% Parameters
nCond = 2;  % congruent(1), incongruent(2)
nStd = 3;  % standard deviation from mean for excluding outliers

tmpSubjID1 = 1:50;
tmpSubjID2 = 1:62;
subjExcld1 = [6, 17, 37, 21, 33];  % based on Hedge et al.'s Readme.txt; additionally exclude 21 33
subjExcld2 = [25, 28, 34, 38, 56, 32];  % based on Readme.txt; additionally exclude 32
subjID1 = setdiff(tmpSubjID1, subjExcld1);
subjID2 = setdiff(tmpSubjID2, subjExcld2);
nSubj1 = length(subjID1);  % study 1
nSubj2 = length(subjID2);  % study 2
nSubjT = nSubj1+nSubj2;
nTask = 1;  % just use the flanker task data
taskStrng = {'Flanker','Stroop'};
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

%% Simulation - Rank order comparison between CE and incongruent trials (Extended Data Figure 2C, 2F)
% Preassignment
nsSubj = 500;
nsTrials = 3200;  % the amount sufficient for highly reliable results
ceMr = nan(nsSubj,numTest);  % CE RT
IncMr = nan(nsSubj,numTest);  % incongruent trial mean
ceMa = nan(nsSubj,numTest);  % accuracy
IncMa = nan(nsSubj,numTest);

% Between-subject distribution
%bs_dist_RT_c = X1_r;
%bs_dist_RT_i = X2_r;
bs_dist_Acc_c = X1_a;  % accuracy
bs_dist_Acc_i = X2_a;
bs_dist_Acc_c(find(bs_dist_Acc_c>1)) = [];
bs_dist_Acc_i(find(bs_dist_Acc_i>1)) = [];
minI = min(length(bs_dist_Acc_c),length(bs_dist_Acc_i));
for l = 1:numTest
    randIDr = randperm(nsSubj);
    randIDa = randperm(nsSubj);
    for i = 1:nsSubj
        % True mean
        ws_mean_RT_c = bs_dist_RT_c(randIDr(i));
        ws_mean_RT_i = bs_dist_RT_i(randIDr(i));
        ws_mean_Acc_c = bs_dist_Acc_c(randIDa(i));
        ws_mean_Acc_i = bs_dist_Acc_i(randIDa(i));

        disp(['Subject# ' num2str(i) '[' num2str(nsSubj) '] - Testing ' num2str(l) '[' num2str(numTest) ']'])
        % within-subject distribution
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
        for j = 1:nsTrials  % number of trials drawn
            rtC = [rtC; ws_dist_RT_c(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];  % add random noise
            rtI = [rtI; ws_dist_RT_i(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];
            accC = [accC; ws_dist_Acc_c(randi(length(ws_dist_Acc_c),1,1))+normrnd(0,nzsiga)];
            accI = [accI; ws_dist_Acc_i(randi(length(ws_dist_Acc_i),1,1))+normrnd(0,nzsiga)];
        end
        %ceMr(i,l) = trimmean(rtI,15)-trimmean(rtC,15);
        %ceMr(i,l) = mean(rmoutliers(rtI,"mean"))-mean(rmoutliers(rtC,"mean"));
        ceMr(i,l) = mean(rtI)-mean(rtC);
        %IncMr(i,l) = trimmean(rtI,15);
        %IncMr(i,l) = mean(rmoutliers(rtI,"mean"));
        IncMr(i,l) = mean(rtI);
        ceMa(i,l) = 100.*((1-mean(accI))-(1-mean(accC)));  % proportion correct (%)
        IncMa(i,l) = 100.*(1-mean(accI));
    end
end
save('rankcompareCEincSIM','ceMr','IncMr','ceMa','IncMa')
%load rankcompareCEincSIM

%% Plot
palette = parula(3);
% (1) RT (Extended Data Figure 2C)
[~,p1] = sort(mean(IncMr,2),'descend');
r1 = 1:nsSubj;
r1(p1) = r1;
[~,p2] = sort(mean(ceMr,2),'descend');
r2 = 1:nsSubj;
r2(p2) = r2;
% Kendall's tau
TAUcei = corr(r1',r2','type','kendall');
% ICC
ICCcei = ICC([r1',r2'],'A-k',0.05);  % 'A-k'
figure
scatter(mean(IncMr,2),mean(ceMr,2),'filled','MarkerFaceColor',palette(1,:),...
    'MarkerEdgeColor',[.6 .6 .6])
refline
axis square
set(gca,'FontSize',16)
str = sprintf('  Tau = %1.2f\n  ICC = %1.2f', TAUcei, ICCcei);
rText = text(min(get(gca,'xlim')),max(get(gca,'ylim')), str);
set(rText, 'fontsize', 16, 'verticalalignment', 'top', 'horizontalalignment', 'left');
%xlim([0.3 1.1])
xlabel('Incongruent trial reaction time (sec)','FontSize',16)
%ylim([0 0.2])
ylabel('Congruency effect reaction time (sec)','FontSize',16)
%title('C) Simulated data: Extended sampling','FontSize',19)

% (2) Percent correct (Extended Data Figure 2F)
clear r1 r2
[~,p1] = sort(mean(IncMa,2),'descend');
r1 = 1:nsSubj;
r1(p1) = r1;
[~,p2] = sort(mean(ceMa,2),'descend');
r2 = 1:nsSubj;
r2(p2) = r2;
% Kendall's tau
TAUcei = corr(r1',r2','type','kendall');
% ICC
ICCcei = ICC([r1',r2'],'A-k',0.05);
figure
scatter(mean(IncMa,2),mean(ceMa,2),'filled','MarkerFaceColor',palette(1,:),...
    'MarkerEdgeColor',[.6 .6 .6])
refline
axis square
set(gca,'FontSize',16)
str = sprintf('  Tau = %1.2f\n  ICC = %1.2f', TAUcei, ICCcei);
rText = text(min(get(gca,'xlim')),max(get(gca,'ylim')), str);
set(rText, 'fontsize', 16, 'verticalalignment', 'top', 'horizontalalignment', 'left');
%xlim([0 50])
xlabel('Incongruent trial percent correct (%)','FontSize',16)
ylabel('Congruency effect percent correct (%)','FontSize',16)
%title('F) Simulated data: Extended sampling','FontSize',19)
