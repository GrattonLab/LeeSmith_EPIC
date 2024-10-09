%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
%
% To run this script:
% You need Hedge et al.'s (2018) flanker and Stroop task data
% (download at https://osf.io/cwzds/)
% You also need PP_CEmat.mat to extrapolate prime-probe task data (to
% create this file, run EPIC_preprocess_primeprobe.m)
% You also need ICC.m (download at
% https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
% What this script does:
% Simulates Hedge et al.'s (2018) flanker and Stroop task data,
% extrapolates prime-probe task data, and performs confirmatory factor
% analysis (CFA) to examine how sampling size affects the reliability of
% factor analysis
%
% What the script outputs: Supp. Fig. 21
% Plot the test-retest reliability of factor scores
%
% Created on 04/26/2023 by HJ Lee
% Last modified on 07/07/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')

%% Parameter settings
nCond = 2;  % experimental conditions: congruent/incongruent
nStd = 3;  % STD from mean for excluding outliers
UpB = 97.5;  % for CI upper bound
LwB = 2.5;  % lower bound

tmpSubjID1 = 1:50;
tmpSubjID2 = 1:62;
subjExcld1 = [6, 17, 37, 21, 33];  % based on Hedge et al.'s Readme.txt; also removed 21 33
subjExcld2 = [25, 28, 34, 38, 56, 32];  % based on Readme.txt; also removed 32
subjID1 = setdiff(tmpSubjID1, subjExcld1);
subjID2 = setdiff(tmpSubjID2, subjExcld2);
nSubj1 = length(subjID1);  % study 1
nSubj2 = length(subjID2);  % study 2
nSubjT = nSubj1+nSubj2;
nTaskH = 2;  % number of task of Hedge et al.'s data
nTask = 3;  % number of task for CFA
taskStrng = {'Flanker','Stroop'};
taskStrng_p = {'Flanker','Prime-Probe','Stroop'};
%nSession = 2;
nBlocks = 5;
nTrials = 144;
%nTrialsReal = 96;  % number of trials for calculating CE (excluding neutral)

nTdrawn = [50, 100, 200, 400, 800, 1600, 3200]/2;  % number of trials drawn (per condition)
drawLH = 4;  % number of trials drawn for Hedge et al.'s data; max 480 trials (cf. nTdrawn)
n_subs = 100;  % simulated subjects
bs_dist_n = 5000;  % between-subject distribution
ws_dist_n = 10000;  % within-subject
numTest = 100;

nzsig = [0.3, 0.8, 0.8];  % noise sigma values determined after evaluating SSE
crossTaskcorr = [0.1, 0.4, 0.6, 0.8, 1];  % pre-determined cross-task correlation
nF = 1;  % number of factors; it's either 1 or 2 (bc there are only 3 tasks)

%% 1. Load and preprocess Hedge et al.s (2018) flanker and Stroop task data
GrtMat = nan(nCond,nSubjT,nTaskH);  % grand mean RT; need this for simulation (copulas)
GstdMat = nan(nCond,nSubjT,nTaskH);  % standard deviation
ciMatCEr = nan(nSubjT,nTaskH);  % 95% CI of the CE; for simulation (within-subject variance)
% Study 1 data
for i = 1:nSubj1
    for t = 1:nTaskH  % 2
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
        cd('../../')  % back to home directory
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
        GrtMat(:,i,t) = TGMrt.mean_rt;
        GstdMat(:,i,t) = TGSrt.std_rt;

        % descriptive
        rtC = T.rt(and(T.n0cong==0,T.acc==1));  % excludes incorrect and outlier trials
        rtI = T.rt(and(T.n0cong==2,T.acc==1));
        rtCE = [];
        for nt = 1:1000  % bootstrap
            tp1 = datasample(rtC,length(rtC));
            tp2 = datasample(rtI,length(rtI));
            rtCE = [rtCE; mean(tp2)-mean(tp1)];
        end
        ciMatCEr(i,t) = (prctile(rtCE,UpB)-prctile(rtCE,LwB))/2;
        clear rtCE tp1 tp2 nt
    end
end
% Study 2 data
for i = 1:nSubj2
    for t = 1:nTaskH
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
        GrtMat(:,nSubj1+i,t) = TGMrt.mean_rt;
        GstdMat(:,nSubj1+i,t) = TGSrt.std_rt;

        % descriptive (cont.)
        rtC = T.rt(and(T.n0cong==0,T.acc==1));  % con (1 is neutral)
        rtI = T.rt(and(T.n0cong==2,T.acc==1));  % inc
        rtCE = [];
        for nt = 1:1000  % bootstrap
            tp1 = datasample(rtC,length(rtC));
            tp2 = datasample(rtI,length(rtI));
            rtCE = [rtCE; mean(tp2)-mean(tp1)];
        end
        ciMatCEr(nSubj1+i,t) = (prctile(rtCE,UpB)-prctile(rtCE,LwB))/2;
    end
end
% Parameters for between-subject (estimated) TRUE distribution
% mean
bs_mean_RT = nan(nCond,nTaskH);
bs_mean_RT(1,1) = mean(GrtMat(1,:,1));  % flanker
bs_mean_RT(2,1) = mean(GrtMat(2,:,1));
bs_mean_RT(1,3) = mean(GrtMat(1,:,2));  % Stroop
bs_mean_RT(2,3) = mean(GrtMat(2,:,2));
load PP_CEmat.mat  % BI dataset means
TMPpp = [mean(PP_CrtGrandMean)/1000-0.05; mean(PP_IrtGrandMean)/1000-0.05];  % convert to secs
bs_mean_RT(:,2) = TMPpp;
% std
bs_std_RT = nan(nCond,nTaskH);
bs_std_RT(1,1) = std(GrtMat(1,:,1));  % flanker
bs_std_RT(2,1) = std(GrtMat(2,:,1));
bs_std_RT(1,3) = std(GrtMat(1,:,2));  % Stroop
bs_std_RT(2,3) = std(GrtMat(2,:,2));
TMPstd = [0.052; 0.064];  % extrapolate prime-probe
bs_std_RT(:,2) = TMPstd;
% skewness and kurtosis
tmpD(1,1) = skewness(GrtMat(1,:,2));  % Stroop
tmpD(2,1) = skewness(GrtMat(2,:,2));
tmpD(1,2) = kurtosis(GrtMat(1,:,2));
tmpD(2,2) = kurtosis(GrtMat(2,:,2));

% Extrapolate prime-probe task's between-subject distribution based on Hedge's Stroop task data
GrtMatPP = nan(nCond,nSubjT);
GrtMatPP(1,:) = pearsrnd(bs_mean_RT(1,2),bs_std_RT(1,2),tmpD(1,1),tmpD(1,2),nSubjT,1);  % con
GrtMatPP(2,:) = pearsrnd(bs_mean_RT(2,2),bs_std_RT(2,2),tmpD(2,1),tmpD(2,2),nSubjT,1);  % inc

% Within-subject SE - use this in simulation if CE distributions are
% simulated - SO WILL USE THIS ONE, NOT THE ABOVE VALUES
m_ws_std_CE(1,1) = mean(ciMatCEr(:,1));  % flanker
m_ws_std_CE(3,1) = mean(ciMatCEr(:,2));  % Stroop
m_ws_std_CE(2,1) = 0.02;

%% 2. Simulate data
%% Copulas: Simulating CE distributions of three tasks with a pre-set cross-task correlations
% Xce = cell(nTask,1);
% % Step 1. Fit to empirical data
% % Flanker
% [Fi1, xi1] = ecdf(GrtMat(2,:,1)-GrtMat(1,:,1));
% Fi1_sm = ksdensity(GrtMat(2,:,1)-GrtMat(1,:,1),xi1,'function','cdf','width',.005);
% figure
% subplot(1,nTask,1)
% stairs(xi1,Fi1,'b','LineWidth',2); hold on
% plot(xi1,Fi1_sm,'r-','LineWidth',1.5)
% xlabel('X1')
% ylabel('Cumulative Probability')
% legend('Empirical','Smoothed','Location','NW'); grid on
% title('Flanker')
% % Prime-probe
% [Fi2, xi2] = ecdf(GrtMatPP(2,:)-GrtMatPP(1,:));
% Fi2_sm = ksdensity(GrtMatPP(2,:)-GrtMatPP(1,:),xi2,'function','cdf','width',.01);
% subplot(1,nTask,2)
% stairs(xi2,Fi2,'b','LineWidth',2); hold on
% plot(xi2,Fi2_sm,'r-','LineWidth',1.5)
% xlabel('X1')
% ylabel('Cumulative Probability')
% legend('Empirical','Smoothed','Location','NW'); grid on
% title('Prime-probe')
% % Stroop
% [Fi3, xi3] = ecdf(GrtMat(2,:,2)-GrtMat(1,:,2));
% Fi3_sm = ksdensity(GrtMat(2,:,2)-GrtMat(1,:,2),xi3,'function','cdf','width',.01);
% subplot(1,nTask,3)
% stairs(xi3,Fi3,'b','LineWidth',2); hold on
% plot(xi3,Fi3_sm,'r-','LineWidth',1.5)
% xlabel('X1')
% ylabel('Cumulative Probability')
% legend('Empirical','Smoothed','Location','NW'); grid on
% title('Stroop')
% % Step 2. Calculate correlation
% nDp = n_subs;
% crossTaskcorr = [0.1, 0.4, 0.6, 0.8, 1];
% nR = length(crossTaskcorr);
% for rI = 1:nR
%     rho12 = crossTaskcorr(rI);
%     rho23 = crossTaskcorr(rI);
%     rho31 = crossTaskcorr(rI);
%     Rho = [1 rho12 rho31; rho12 1 rho23; rho31 rho23 1];
%     U = copularnd('t', Rho, 2, nDp);
%     Xce{1} = ksdensity(GrtMat(2,:,1)-GrtMat(1,:,1),U(:,1),'function','icdf','width',.005);
%     Xce{3} = ksdensity(GrtMat(2,:,2)-GrtMat(1,:,1),U(:,3),'function','icdf','width',.01);
%     Xce{2} = ksdensity(GrtMatPP(2,:)-GrtMatPP(1,:),U(:,2),'function','icdf','width',.01);
%     figure
%     subplot(1,1,1)
%     plot3(Xce{1},Xce{2},Xce{3},'.')
%     grid on
%     view([-55, 15])
%     xlabel('Flanker')
%     ylabel('Prime-probe')
%     zlabel('Stroop')
%     title(['Cross-task correlation ' num2str(crossTaskcorr(rI))])
% end

%% Run this script from here for five times with different rI values
for rI = 1:length(crossTaskcorr)
    rho12 = crossTaskcorr(rI);
    rho23 = crossTaskcorr(rI);
    rho31 = crossTaskcorr(rI);
    Rho = [1 rho12 rho31; rho12 1 rho23; rho31 rho23 1];

    %% Simulate to get test-retest reliability across varying number of trials
    % Preassignment
    % Subject X tested in day 1
    simCErt_cell = cell(length(nTdrawn),1);  % mean CE RT
    lambdaM = nan(nTask,numTest,length(nTdrawn));  % factor loadings
    psiM = nan(nTask,numTest,length(nTdrawn));  % specific variance
    FM = nan(n_subs,numTest,length(nTdrawn));  % factor scores

    % (same subject) tested in day 2
    simCErt_cell2 = cell(length(nTdrawn),1);
    lambdaM2 = nan(nTask,numTest,length(nTdrawn));
    psiM2 = nan(nTask,numTest,length(nTdrawn));
    FM2 = nan(n_subs,numTest,length(nTdrawn));

    for k = 1:length(nTdrawn)
        simCErt_mat = nan(n_subs,nTask,numTest);  % day 1
        simCErt_mat2 = nan(n_subs,nTask,numTest);  % day 2
        for l = 1:numTest
            U = copularnd('t', Rho, 2, bs_dist_n);
            Xce{1} = ksdensity(GrtMat(2,:,1)-GrtMat(1,:,1),U(:,1),'function','icdf','width',.005);
            Xce{3} = ksdensity(GrtMat(2,:,2)-GrtMat(1,:,1),U(:,3),'function','icdf','width',.01);
            Xce{2} = ksdensity(GrtMatPP(2,:)-GrtMatPP(1,:),U(:,2),'function','icdf','width',.01);
            for n = 1:n_subs
                myrank = randi(bs_dist_n,1,1);
                for t = 1:nTask  % 3 tasks
                    disp(['Sampling step:' num2str(k) ', Simulation' num2str(l) ', Task' num2str(t) ' Subject#' num2str(n)])
                    bs_dist_RT = Xce{t};  % between-subject distribution
                    ws_mean_RT = bs_dist_RT(myrank);  % the task-loop should be inside the subject-loop so that the subject's rank ordering is consistent across tasks

                    % within-subject distribution
                    % Day 1
                    ws_dist_RT = randn(1,ws_dist_n)*m_ws_std_CE(t,1)+ws_mean_RT;
                    ce = [];
                    for d = 1:nTdrawn(k)  % trials drawn
                        ce = [ce; ws_dist_RT(randi(ws_dist_n,1,1))+normrnd(0,nzsig(t))];  % add random noise
                    end
                    %simCErt_mat(n,t,l) = trimmean(ce,15);  % CE
                    simCErt_mat(n,t,l) = mean(rmoutliers(ce,"mean"));  % remove outliers 3 SD away from mean
                    %simCErt_mat(n,t,l) = mean(ce);  % doesn't remove outliers
                    % Day 2
                    ws_dist_RT = randn(1,ws_dist_n)*m_ws_std_CE(t,1)+ws_mean_RT;
                    ce = [];
                    for d = 1:nTdrawn(k)  % trials drawn
                        ce = [ce; ws_dist_RT(randi(ws_dist_n,1,1))+normrnd(0,nzsig(t))];  % add random noise
                    end
                    %simCErt_mat2(n,t,l) = trimmean(ce,15);
                    simCErt_mat2(n,t,l) = mean(rmoutliers(ce,"mean"));
                    %simCErt_mat2(n,t,l) = mean(ce);
                end
            end
            % CFA
            [lambda, psi, T, stats, F] = factoran(simCErt_mat(:,:,l),nF);
            lambdaM(:,l,k) = lambda;
            psiM(:,l,k) = psi;
            FM(:,l,k) = F;
            clear lambda psi T stats F
            [lambda,psi,T,stats,F] = factoran(simCErt_mat2(:,:,l),nF);
            lambdaM2(:,l,k) = lambda;
            psiM2(:,l,k) = psi;
            FM2(:,l,k) = F;
            clear lambda psi T stats F
        end
        simCErt_cell{k,1} = simCErt_mat;
        simCErt_cell2{k,1} = simCErt_mat2;
    end
    %save(['CFAsim_testretestResult_avg_corr' num2str(rI)],'simCErt_cell','simCErt_cell2', ...
    %    'lambdaM','psiM','FM','lambdaM2','psiM2','FM2')
end

%% Plot
% Colormap
%cmap1 = parula(3);
%cmap2 = spring(3);
cmap3 = cool(3);

% Preassignment
ICCmat = nan(length(nTdrawn),length(crossTaskcorr)*(nTask+1));  % 7*20
ICCmat2 = nan(length(nTdrawn),length(crossTaskcorr)*(nTask+1));  % in different format
for rI = 1:length(crossTaskcorr)
    dataname = ['CFAsim_testretestResult_avg_corr' num2str(rI) '.mat'];
    load(dataname)
    if rI == 1
        figure
    end
    for t = 1:nTask
        for k = 1:length(nTdrawn)
            simCE = simCErt_cell{k,1};
            simCE2 = simCErt_cell2{k,1};
            sCE = mean(squeeze(simCE(:,t,:)),2);  % n_subs*numTest; 2 bc you want each subject's 100 simulation mean
            sCE2 = mean(squeeze(simCE2(:,t,:)),2);
            subplot(length(nTdrawn),length(crossTaskcorr)*(nTask+1),length(crossTaskcorr)*(nTask+1)*(k-1)+(t+(rI-1)*(nTask+1)))
            scatter(sCE,sCE2,8,'filled')
            %lsline
            xlim([min([sCE(:);sCE2(:)]) max([sCE(:);sCE2(:)])])
            ylim([min([sCE(:);sCE2(:)]) max([sCE(:);sCE2(:)])])
            [ICc] = ICC([sCE,sCE2],'A-k',0.05);
            ICCmat(k,length(crossTaskcorr)*(t-1)+rI) = ICc;
            ICCmat2(k,t+(rI-1)*(nTask+1)) = ICc;
            str = sprintf('  ICC = %1.2f',ICc); clear ICc
            rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
            set(rText,'fontweight','bold','fontsize', 7.4,'verticalalignment','top','horizontalalignment','left');
            refline
            axis square
            if and(t == 1,rI==1)
                ylabel(['N = ' num2str(nCond*nTdrawn(k))],'FontWeight','bold','FontSize',9.5)
            end
        end
    end
    subplot(length(nTdrawn),length(crossTaskcorr)*(nTask+1),1+(rI-1)*(nTask+1)); title('Flanker CE','FontSize',8)
    subplot(length(nTdrawn),length(crossTaskcorr)*(nTask+1),2+(rI-1)*(nTask+1)); title('Prime-Probe CE','FontSize',8)
    subplot(length(nTdrawn),length(crossTaskcorr)*(nTask+1),3+(rI-1)*(nTask+1)); title('Stroop CE','FontSize',8)
    % CFA parameters - factor score
    for k = 1:length(nTdrawn)
        fs = mean(squeeze(FM(:,:,k)),2);
        fs2 = mean(squeeze(FM2(:,:,k)),2);
        subplot(length(nTdrawn),length(crossTaskcorr)*(nTask+1),length(crossTaskcorr)*(nTask+1)*(k-1)+(4+(rI-1)*(nTask+1)))
        scatter(fs,fs2,8,cmap3(3,:),'filled')
        %lsline
        xlim([min([fs(:);fs2(:)]) max([fs(:);fs2(:)])])
        ylim([min([fs(:);fs2(:)]) max([fs(:);fs2(:)])])
        [ICc] = ICC([fs,fs2],'A-k',0.05);
        ICCmat(k,length(crossTaskcorr)*nTask+rI) = ICc;
        ICCmat2(k,4+(rI-1)*(nTask+1)) = ICc;
        str = sprintf('  ICC = %1.2f',ICc); clear ICc
        rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
        set(rText,'fontweight','bold','fontsize',7,'verticalalignment','top','horizontalalignment','left');
        refline
        axis square
    end
end
subplot(length(nTdrawn),length(crossTaskcorr)*(nTask+1),4+(rI-1)*(nTask+1)); title('Factor Score')

%% Plot heatmap
figure
h = heatmap(round(ICCmat.*100)/100);  % size 7*20
h.XLabel = 'Cross-task correlation';
h.YLabel = 'Number of trials';
Xlabels = repmat(crossTaskcorr,1,nTask+1);
Ylabels = nCond*nTdrawn;
h.XDisplayLabels = Xlabels;
h.YDisplayLabels = Ylabels;
%h.Title = '';
h.Colormap = jet;  % parula

figure
h = heatmap(round(ICCmat2.*100)/100);  % size 7*20
%h.XLabel = 'Cross-task correlation';
h.YLabel = 'Number of trials';
Xlabels = repmat({'Flanker','Prime-probe','Stroop','Factor score'},1,length(crossTaskcorr));
Ylabels = nCond*nTdrawn;
%h.XDisplayLabels = Xlabels;
h.YDisplayLabels = Ylabels;
%h.Title = '';
h.Colormap = jet;
h.CellLabelFormat = '%.2f';
h.FontSize = 12;