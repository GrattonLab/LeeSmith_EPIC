%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC - Public Data 2
%
%% Instructions for running the script:
% You will need the following datasets and files.
%   - Hedge et al.'s (2018) flanker and Stroop task data (https://osf.io/cwzds/)
%   - Mat-file: ICC.m (https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
%% Purpose of the script:
% This script fits the model to Hedge et al.'s (2018) flanker and Stroop 
% task data across varying noise sigma values to best replicate the ICCs 
% of the empirical data by minimizing the sum of squared errors.
% Noise is added when simulating data for confirmatory factor analysis
% (CFA) to examine how sampling size affects the reliability of 
% factor analysis.
%
%% Outputs:
% - Noise sigma value
% - Heatmap of the sum of squared error across different noise sigma values
%
% Created on 04/26/2023 by HJ Lee
% Last modified on 02/04/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')

%% Parameters
nCond = 2;  % congruent(1), incongruent(2)
nStd = 3;  % standard deviation from mean for excluding outliers
UpB = 97.5;  % upper percentile of 95% confidence interval
LwB = 2.5;  % lower percentile

tmpSubjID1 = 1:50;
tmpSubjID2 = 1:62;
subjExcld1 = [6, 17, 37, 21, 33];  % based on Hedge et al.'s Readme.txt; additionally exclude 21 33
subjExcld2 = [25, 28, 34, 38, 56, 32];  % based on Readme.txt; additionally exclude 32
subjID1 = setdiff(tmpSubjID1, subjExcld1);
subjID2 = setdiff(tmpSubjID2, subjExcld2);
nSubj1 = length(subjID1);  % study 1
nSubj2 = length(subjID2);  % study 2
nSubjT = nSubj1+nSubj2;
nTask = 2;
taskStrng = {'Flanker','Stroop'};
nBlocks = 5;
nTrials = 144;
%nTrialsCE = 96;  % number of trials for calculating CE (excluding neutral)

nTdrawn = [50, 100, 200, 400, 800, 1600, 3200]/2;  % number of trials drawn (per condition)
drawLH = 4;  % number of trials drawn for Hedge et al.'s data; max 480 trials (cf. nTdrawn)
n_subs = nSubjT;  % number of simulated subjects; match with empirical data for ease of comparison
bs_dist_n = 5000;  % between-subject distribution
ws_dist_n = 10000;  % within-subject
numTest = 100;  % repeat simulations

% For FITTING - noise sigma to simulate trial variability
% Adjust the noise sigma by calculating the sum of squared error
rI = 5;  % change this from 1-5; this is the index for "crossTaskcorr"
crossTaskcorr = [0.1, 0.4, 0.6, 0.8, 1];  % pre-determined cross-task correlation
nzsigV = 0:0.05:1;  % 0:0.01:1
nTestingr = length(nzsigV);

%% Step 1. Load and preprocess Hedge et al.'s (2018) flanker and Stroop task data
GrtMat = nan(nCond,nSubjT,nTask);  % grand mean RT; need this for simulation (copulas)
ciMatCEr = nan(nSubjT,nTask);  % 95% CI of the CE RT; for simulation (within-subject variance)
HrtMatc1 = nan(length(nTdrawn),nSubjT,nTask);  % mean congruent trial RT in secs; split-half grand mean; need this to compare w/ simulated data
HrtMati1 = nan(length(nTdrawn),nSubjT,nTask);  % incongruent
HrtMatc2 = nan(length(nTdrawn),nSubjT,nTask);
HrtMati2 = nan(length(nTdrawn),nSubjT,nTask);
% Study 1 data
for i = 1:nSubj1
    for t = 1:nTask  % 2
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
        GrtMat(:,i,t) = TGMrt.mean_rt;

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
        rtC1 = T1.rt(and(T1.n0cong==0,T1.acc==1));
        rtI1 = T1.rt(and(T1.n0cong==2,T1.acc==1));
        rtC2 = T2.rt(and(T2.n0cong==0,T2.acc==1));
        rtI2 = T2.rt(and(T2.n0cong==2,T2.acc==1));
        for k = 1:drawLH
            if length(rtC1) < nTdrawn(k)  % this happens when k == 4
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtC1))])
                HrtMatc1(k,i,t) = mean(rtC1);  % don't need trimmean because outlier trials were marked as 99 above
            else
                HrtMatc1(k,i,t) = mean(rtC1(1:nTdrawn(k)));
            end
            if length(rtI1) < nTdrawn(k)
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtI1))])
                HrtMati1(k,i,t) = mean(rtI1);
            else
                HrtMati1(k,i,t) = mean(rtI1(1:nTdrawn(k)));
            end
            if length(rtC2) < nTdrawn(k)
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtC2))])
                HrtMatc2(k,i,t) = mean(rtC2);
            else
                HrtMatc2(k,i,t) = mean(rtC2(1:nTdrawn(k)));
            end
            if length(rtI2) < nTdrawn(k)
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtI2))])
                HrtMati2(k,i,t) = mean(rtI2);
            else
                HrtMati2(k,i,t) = mean(rtI2(1:nTdrawn(k)));
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
        GrtMat(:,nSubj1+i,t) = TGMrt.mean_rt;

        % descriptive (cont.)
        rtC = T.rt(and(T.n0cong==0,T.acc==1));
        rtI = T.rt(and(T.n0cong==2,T.acc==1));
        rtCE = [];
        for nt = 1:1000  % bootstrap
            tp1 = datasample(rtC,length(rtC));
            tp2 = datasample(rtI,length(rtI));
            rtCE = [rtCE; mean(tp2)-mean(tp1)];
        end
        ciMatCEr(nSubj1+i,t) = (prctile(rtCE,UpB)-prctile(rtCE,LwB))/2;
        clear rtCE accCE tp1 tp2 tp3 tp4 nt
        rtC1 = T1.rt(and(T1.n0cong==0,T1.acc==1));
        rtI1 = T1.rt(and(T1.n0cong==2,T1.acc==1));
        rtC2 = T2.rt(and(T2.n0cong==0,T2.acc==1));
        rtI2 = T2.rt(and(T2.n0cong==2,T2.acc==1));
        for k = 1:drawLH
            if length(rtC1) < nTdrawn(k)  % this happens when k == 4
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtC1))])
                HrtMatc1(k,nSubj1+i,t) = mean(rtC1);  % don't need trimmean because outlier trials were marked as 99 above
            else
                HrtMatc1(k,nSubj1+i,t) = mean(rtC1(1:nTdrawn(k)));
            end
            if length(rtI1) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtI1))])
                HrtMati1(k,nSubj1+i,t) = mean(rtI1);
            else
                HrtMati1(k,nSubj1+i,t) = mean(rtI1(1:nTdrawn(k)));
            end
            if length(rtC2) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtC2))])
                HrtMatc2(k,nSubj1+i,t) = mean(rtC2);
            else
                HrtMatc2(k,nSubj1+i,t) = mean(rtC2(1:nTdrawn(k)));
            end
            if length(rtI2) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtI2))])
                HrtMati2(k,nSubj1+i,t) = mean(rtI2);
            else
                HrtMati2(k,nSubj1+i,t) = mean(rtI2(1:nTdrawn(k)));
            end
        end
    end
end

% Within-subject SE - for simulation
m_ws_std_CE(1,1) = mean(ciMatCEr(:,1));  % flanker
m_ws_std_CE(2,1) = mean(ciMatCEr(:,2));  % Stroop

%% Step 2. Simulate data and fit model
%% Copulas: Simulating CE distributions with a pre-set cross-task correlations
% Xce = cell(nTask,1);
% % 1. Fit to empirical data
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
rho12 = crossTaskcorr(rI);  % adjust rI: 1-5 (line 59)
rho23 = crossTaskcorr(rI);
rho31 = crossTaskcorr(rI);
Rho = [1 rho12 rho31; rho12 1 rho23; rho31 rho23 1];

%% Find the sigma size of random noise that minimizes the sum of squared error of ICC
nzMV = nan(numTest,nTask);  % noise sigma value; the ones that have minimum SSE
sseFL = nan(nTestingr,drawLH,numTest);  % sum of squared error of ICC; flanker task
sseST = nan(nTestingr,drawLH,numTest);  % Stroop
for l = 1:numTest
    ICCe = nan(nTestingr,drawLH,nTask);  % ICC of empirical data
    ICCs = nan(nTestingr,drawLH,nTask);  % ICC of simulated data

    % Between-subject joint distributions (having pre-set cross-task correlations)
    U = copularnd('t', Rho, 2, bs_dist_n);
    Xce{1} = ksdensity(GrtMat(2,:,1)-GrtMat(1,:,1),U(:,1),'function','icdf','width',.005);  % flanker
    Xce{2} = ksdensity(GrtMat(2,:,2)-GrtMat(1,:,1),U(:,3),'function','icdf','width',.01);  % Stroop
    for t = 1:nTask
        disp(['Sigma evaluation starts for task ' num2str(t) '...'])
        bs_dist_RT = Xce{t};  % between-subject distribution
        ceDiffe1 = squeeze(HrtMati1(:,:,t)-HrtMatc1(:,:,t));  % RT CE empirical data
        ceDiffe2 = squeeze(HrtMati2(:,:,t)-HrtMatc2(:,:,t));
        for s = 1:nTestingr
            nzsig = nzsigV(s);  % noise sigma
            % Preassignment
            ceDiff1 = nan(drawLH,n_subs);
            ceDiff2 = nan(drawLH,n_subs);
            for k = 1:drawLH
                for n = 1:n_subs
                    disp(['Sampling step:' num2str(k) ', Subject#' num2str(n)])
                    myrank = randi(bs_dist_n,1,1);
                    ws_mean_RT = bs_dist_RT(myrank);
                    % within-subject distribution
                    % Day 1
                    ws_dist_RT = randn(1,ws_dist_n)*m_ws_std_CE(t,1)+ws_mean_RT;
                    ce = [];
                    for d = 1:nTdrawn(k)  % trials drawn
                        ce = [ce; ws_dist_RT(randi(ws_dist_n,1,1))+normrnd(0,nzsig)];  % add random noise
                    end
                    ceDiff1(k,n) = mean(ce); clear ce
                    % Day 2
                    ws_dist_RT = randn(1,ws_dist_n)*m_ws_std_CE(t,1)+ws_mean_RT;
                    ce = [];
                    for d = 1:nTdrawn(k)  % trials drawn
                        ce = [ce; ws_dist_RT(randi(ws_dist_n,1,1))+normrnd(0,nzsig)];  % add random noise
                    end
                    ceDiff2(k,n) = mean(ce); clear ce
                end
                % ICC
                [Icce] = ICC([ceDiffe1(k,:)',ceDiffe2(k,:)'],'A-k',0.05);
                [Iccs] = ICC([ceDiff1(k,:)',ceDiff2(k,:)'],'A-k',0.05);
                ICCe(s,k,t) = Icce;
                ICCs(s,k,t) = Iccs;
                % difference in ICC - SSE
                %sumllh(s,k,t) = log(nzsig)+0.5*log(2*pi*ones(1,1))+(((ICCe(s,k,t)-ICCs(s,k,t)).^2)./(2*nzsig.^2));
                if t == 1
                    sseFL(s,k,l) = (ICCs(s,k,t)-ICCe(s,k,t))^2;
                elseif t == 2
                    sseST(s,k,l) = (ICCs(s,k,t)-ICCe(s,k,t))^2;
                end
            end
        end
        if t == 1
            sumSSE = sum(sseFL(:,:,l),2);
            [m,id] = min(sumSSE);
        elseif t == 2
            sumSSE = sum(sseST(:,:,l),2);
            [m,id] = min(sumSSE);
        end
        nzMV(l,t) = nzsigV(id);
    end
end
disp(['[Cross-task correlation ' num2str(crossTaskcorr(rI)) '] Sigma value that maximizes log-likelihood for task ' num2str(1) 'across 100 simulations is ' num2str(mean(nzMV(:,1)))])
disp(['[Cross-task correlation ' num2str(crossTaskcorr(rI)) '] Sigma value that maximizes log-likelihood for task ' num2str(2) 'across 100 simulations is ' num2str(mean(nzMV(:,2)))])

%% Plot the output
figure
a = [mean(squeeze(sum(sseFL,2)),2),mean(squeeze(sum(sseST,2)),2)];
h = heatmap(a);
h.Xlabel = 'Task';
h.Ylabel = 'Noise sigma';
Xlabels = {'Flanker','Stroop'};
Ylabels = nzsigV;
h.XDisplayLabels = Xlabels;
h.YDisplayLabels = Ylabels;
h.Title = 'Noise sigma Sum of Squared Error (average across 100 simulations)';
h.Colormap = parula;
saveas(h,'nzsigSSE_heatmap_CFA.fig')