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
% Fits the model to Hedge et al.'s (2018) flanker task data across varying 
% sizes of noise sigma to best replicate the ICCs of the empirical data 
% (via minimizing the sum of squared error)
% The noise will be added when simulating data for EZ-diffusion modeling 
% and analyses on rank order consistency
%
% What this script outputs:
% The noise sigma values for RT and accuracy & the heatmap of the sum of
% squared error across different noise sigma values
%
% Created on 05/16/2023 by HJ Lee
% Last modified on 07/05/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
drawLH = 4;  % number of trials drawn for Hedge et al's data; max 480 trials (cf. nTdrawn)
n_subs = nSubjT;  % number of simulated subjects; match with empirical data for ease of comparison
bs_dist_n = 5000;
ws_dist_n = 10000;
numTest = 100;  % repeat simulations

% For FITTING - noise sigma to simulate trial variability
% Adjust the noise sigma by calculating the sum of squared error
nzsigVr = 0.02:0.04:0.82;  % rt; run once again by narrowing down to 0.1:0.04:0.42
nzsigVa = nzsigVr;  % acc; have this to be equal length to that of RT as I'll use only one for-loop here
nTestingr = length(nzsigVr);
nTestinga = length(nzsigVa);
if nTestingr ~= nTestinga
    error('Set nTestingr and nTestinga to have equal lengths')
end

%% Step 1. Load and preprocess Hedge et al.'s (2018) flanker task data
GrtMat = nan(nCond,nSubjT);  % grand mean RT; need this for simulation (copulas)
GaccMat = nan(nCond,nSubjT);  % accuracy
stdMatr = nan(nCond,nSubjT);  % within-subject STD; need this also for simulation
stdMata = nan(nCond,nSubjT);
HrtMat1 = nan(nCond,drawLH,nSubjT);  % mean RT in secs; split-half grand mean; need this to compare w/ simulated data
HrtMat2 = nan(nCond,drawLH,nSubjT);
HrtvMat1 = nan(nCond,drawLH,nSubjT);  % RT variance
HrtvMat2 = nan(nCond,drawLH,nSubjT);
HaccMat1 = nan(nCond,drawLH,nSubjT);  % mean accuracy
HaccMat2 = nan(nCond,drawLH,nSubjT);
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
                HrtMat1(1,k,nSubj1+i) = mean(rtC1);
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
% CE across varying numbers of trials - for noise sigma evaluation
ceDiffe1 = squeeze(HrtMat1(2,:,:)-HrtMat1(1,:,:));  % RT; session 1
ceDiffe2 = squeeze(HrtMat2(2,:,:)-HrtMat2(1,:,:));  % session 2
ceDiffeA1 = squeeze(HaccMat1(1,:,:)-HaccMat1(2,:,:));  % accuracy
ceDiffeA2 = squeeze(HaccMat2(1,:,:)-HaccMat2(2,:,:));

% Within-subject standard deviation - for simulation
m_ws_std_RT(1,1) = mean(stdMatr(1,:));  % con
m_ws_std_RT(2,1) = mean(stdMatr(2,:));  % inc
m_ws_std_Acc(1,1) = mean(stdMata(1,:));
m_ws_std_Acc(2,1) = mean(stdMata(2,:));

%% Step 2. Simulate data and fit model
%% Copulas: Generating dependent multivariate data to simulate the correlation between congruent and incongruent trials
% RT
% 2. Calculate correlation between congruent and incongruent trials
tau_r = corr(GrtMat(1,:)',GrtMat(2,:)','type','kendall');  % Kendall's rank order
nu = 1;  % degree of freedom
rho_r = copulaparam('t',tau_r,nu,'type','kendall');  % corresponding linear correlation parameter for the t copula
% 3. Generate random values from the t copula
U_r = copularnd('t',[1 rho_r; rho_r 1], nu, bs_dist_n);
X1_r = ksdensity(GrtMat(1,:),U_r(:,1),'function','icdf','width',.01);
X2_r = ksdensity(GrtMat(2,:),U_r(:,2),'function','icdf','width',.01);
% Acc
% 2. Calculate correlation
tau_a = corr(GaccMat(1,:)',GaccMat(2,:)','type','kendall');  % Kendall's rank order
rho_a = copulaparam('t',tau_a,nu,'type','kendall');  % corresponding linear correlation parameter for the t copula
% 3. Generate random values from the t copula
U_a = copularnd('t',[1 rho_a; rho_a 1], nu, bs_dist_n);
X1_a = ksdensity(GaccMat(1,:),U_a(:,1),'function','icdf','width',.008);
X2_a = ksdensity(GaccMat(2,:),U_a(:,2),'function','icdf','width',.01);

%% Now run simulation
% Between-subject joint distributions
bs_dist_RT_c = X1_r;  % RT con
bs_dist_RT_i = X2_r;  % inc
bs_dist_Acc_c = X1_a;  % accuracy
bs_dist_Acc_i = X2_a;
bs_dist_Acc_c(find(bs_dist_Acc_c>1)) = [];
bs_dist_Acc_i(find(bs_dist_Acc_i>1)) = [];
minI = min(length(bs_dist_Acc_c),length(bs_dist_Acc_i));

% Preassignment
ICCre = nan(nTestingr,drawLH,numTest);  % rt; ICC of empirical data
ICCae = nan(nTestinga,drawLH,numTest);  % acc
ICCrs = nan(nTestingr,drawLH,numTest);  % ICC of simulated data
ICCas = nan(nTestinga,drawLH,numTest);
sser = nan(nTestingr,drawLH,numTest);  % rt; sum of squared error
ssea = nan(nTestinga,drawLH,numTest);  % acc
nzMV = nan(2,numTest);  % 2: nzsigr & nzsiga; store the noise sigma values that have minimum SSE
for ts = 1:numTest
    for s = 1:nTestingr
        nzsigr = nzsigVr(s);  % noise sigma
        nzsiga = nzsigVa(s);
        % Preassignment
        simCErt_mat1 = nan(nCond,drawLH,n_subs);  % RT; subject X tested in day 1
        simCErt_mat2 = nan(nCond,drawLH,n_subs);  % (same subject) tested in day 2
        simCEacc_mat1 = nan(nCond,drawLH,n_subs);  % acc
        simCEacc_mat2 = nan(nCond,drawLH,n_subs);
        for k = 1:drawLH
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
                for l = 1:nTdrawn(k)  % number of trials
                    disp(['TESTING ' num2str(ts) ' [sigma step ' num2str(s) '] Sampling step:' num2str(k) ', Subject#' num2str(n) ', Simulation trial#' num2str(l)])
                    rtC = [rtC; ws_dist_RT_c(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];
                    rtI = [rtI; ws_dist_RT_i(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];
                    accC = [accC; ws_dist_Acc_c(randi(length(ws_dist_Acc_c),1,1))+normrnd(0,nzsiga)];
                    accI = [accI; ws_dist_Acc_i(randi(length(ws_dist_Acc_i),1,1))+normrnd(0,nzsiga)];
                end
                %rtC = rmoutliers(rtC,"mean");
                %rtI = rmoutliers(rtI,"mean");
                simCErt_mat1(1,k,n) = mean(rtC);
                simCErt_mat1(2,k,n) = mean(rtI);
                simCEacc_mat1(1,k,n) = mean(accC);
                simCEacc_mat1(2,k,n) = mean(accI);
                clear rtC rtI accC accI

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
                for l = 1:nTdrawn(k) % number of trials
                    disp(['TESTING ' num2str(ts) ' [sigma step ' num2str(s) '] Sampling step:' num2str(k) ', Subject#' num2str(n) ', Simulation trial#' num2str(l)])
                    rtC = [rtC; ws_dist_RT_c(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];
                    rtI = [rtI; ws_dist_RT_i(randi(ws_dist_n,1,1))+normrnd(0,nzsigr)];
                    accC = [accC; ws_dist_Acc_c(randi(length(ws_dist_Acc_c),1,1))+normrnd(0,nzsiga)];
                    accI = [accI; ws_dist_Acc_i(randi(length(ws_dist_Acc_i),1,1))+normrnd(0,nzsiga)];
                end
                %rtC = rmoutliers(rtC,"mean");
                %rtI = rmoutliers(rtI,"mean");
                simCErt_mat2(1,k,n) = mean(rtC);
                simCErt_mat2(2,k,n) = mean(rtI);
                simCEacc_mat2(1,k,n) = mean(accC);
                simCEacc_mat2(2,k,n) = mean(accI);
                clear rtC rtI accC accI
            end
            ceDiff1 = squeeze(simCErt_mat1(2,k,:)-simCErt_mat1(1,k,:));  % test 1 RT inc-con
            ceDiff2 = squeeze(simCErt_mat2(2,k,:)-simCErt_mat2(1,k,:));  % test 2
            ceDiffA1 = squeeze(simCEacc_mat1(1,k,:)-simCEacc_mat1(2,k,:));  % accuracy con-inc
            ceDiffA2 = squeeze(simCEacc_mat2(1,k,:)-simCEacc_mat2(2,k,:));
            % ICC
            % RT
            [Icce] = ICC([ceDiffe1(k,:)',ceDiffe2(k,:)'],'A-k',0.05);
            [Iccs] = ICC([ceDiff1,ceDiff2],'A-k',0.05);
            ICCre(s,k,ts) = Icce;
            ICCrs(s,k,ts) = Iccs; clear Icce Iccs
            % Acc
            [Icce] = ICC([ceDiffeA1(k,:)',ceDiffeA2(k,:)'],'A-k',0.05);
            [Iccs] = ICC([ceDiffA1,ceDiffA2],'A-k',0.05);
            ICCae(s,k,ts) = Icce;
            ICCas(s,k,ts) = Iccs; clear Icce Iccs
            % a) difference in ICC
            %sumllhr(s,k,ts) = log(nzsigr)+0.5*log(2*pi*ones(1,1))+(((ICCrs(s,k,ts)-ICCre(s,k,ts)).^2)./(2*nzsigr.^2));
            %sumllha(s,k,ts) = log(nzsiga)+0.5*log(2*pi*ones(1,1))+(((ICCas(s,k,ts)-ICCae(s,k,ts)).^2)./(2*nzsiga.^2));
            sser(s,k,ts) = (ICCrs(s,k,ts)-ICCre(s,k,ts))^2;
            ssea(s,k,ts) = (ICCas(s,k,ts)-ICCae(s,k,ts))^2;
        end
    end
    % Find the minimum of the SSE
    sumSSEr = sum(sser(:,:,ts),2);  % sum of different nTdrawn
    sumSSEa = sum(ssea(:,:,ts),2);
    [Mr,IDr] = min(sumSSEr,[],'all');
    [Ma,IDa] = min(sumSSEa,[],'all');
    nzMV(:,ts) = [nzsigVr(IDr);nzsigVa(IDa)];
end
save('nzsigSSEresults_EZD','nzMV','sser','ssea','ICCre','ICCae','ICCrs','ICCas')
nzsigVr_avg = mean(nzMV(1,:));  % average across repeated testing (numTest)
nzsigVa_avg = mean(nzMV(2,:));
disp(['Sigma value that maximizes log-likelihood (across simulations) is ' num2str(nzsigVr_avg) ' for RT and ' num2str(nzsigVa_avg) ' for Acc'])
% 0.26 for RT and 0.23 for accuracy

%% Plot the ouput
figure
a = [mean(squeeze(sum(sser,2)),2),mean(squeeze(sum(ssea,2)),2)];
h = heatmap(a);
h.XLabel = 'Measure';
h.YLabel = 'Noise sigma';
Xlabels = {'RT','Acc'};
Ylabels = nzsigVr;
h.XDisplayLabels = Xlabels;
h.YDisplayLabels = Ylabels;
h.Title = 'Noise Sigma Sum of Squared Error (average across 100 simulations)';
h.Colormap = parula;
saveas(h,'nzsigSSE_heatmap_EZD.fig')