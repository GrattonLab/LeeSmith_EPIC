%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
% 3. PrimeProbe task
%
% To run this script:
% You need the EPIC dataset and "session_numbering.xlsx"
% (download both at https://osf.io/jk9nb)
%
% What this script does:
% Preprocessing and statistical analysis
% IV 1: session (1-18)
% IV 2: current trial congruency (n0cong) - congruent(1) vs. incongruent(2)
% DV 1: RT
% DV 2: Accuracy
% DV 3: Inverse Efficiency Score (IES)
% 
% What this script outputs:
% This script calculates the grand mean of all sessions and its 95% CI
% and outputs a mat-file to use it in EPIC_violinPlotGrandmeanPlot.m 
% The script calculates RT mean, RT variance, accuracy mean across sessions
% and outputs an excel file.
%
% Created on 10/02/2022 by HJ Lee
% Last modified on 06/29/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')

%% Parameter settings
% Indexing participant IDs
tmpSubjID = 1:12;
subjExcld = [1,2,11,9];  % excluded participants; #9 was later excluded due to noisy data
subjID = setdiff(tmpSubjID, subjExcld);
nSubj = length(subjID);

% Task information
nTask = 4;
taskIndx = 3;  % 1:Flanker, 2:GoNogo, 3:PrimeProbe, 4:STROOP
taskStrng = 'PrimeProbe';
nBlocks = 4;
numTrials = 96;  % per block
nTrials = nBlocks*numTrials;  % total number of trials
nSession = 18;  % total number of sessions for each task
sessID = 1:nSession;
nCond = 2;  % congruent(1) vs. incongruent(2)
nStd = 3;  % outlier criterion; 3 standard deviation from mean

% 95% confidence interval
UpB = 97.5;
LwB = 2.5;
nSim = 1000;

% Table header
header_output = {'sess','n0cong','RT','var','acc'};

%% Preassignment
% Matrices for plotting mean RT/acc across conditions, weeks, and participants
rtMat = nan(nCond,nSession,nSubj);
accMat = nan(nCond,nSession,nSubj);
seHMat = nan(nSession,nSubj);  % standard error 18 sessions

RTdist = cell(nCond,nSession,nSubj);
accdist = cell(nCond,nSession,nSubj);
cRTdist = cell(nCond,nSession,nSubj);  % RT distribution of correctly responded RTs
inRTdist = cell(nCond,nSession,nSubj);  % RT distribution after removing incorrect trials and outliers
outRTdist = cell(nCond,nSession,nSubj);  % outlier RTs
num_removedT = cell(nSubj,1);  % Track numbers of excluded trials

% for each participant
for i = 1:nSubj
    %% Load session indexing information
    dir = ['rawData/EPIC' num2str(subjID(i))];
    % session information
    opts = detectImportOptions('session_numbering.xlsx');
    opts.Sheet = ['subj' num2str(subjID(i))];
    tmpMat = readmatrix('session_numbering.xlsx',opts);
    sessionID = tmpMat(:,taskIndx);

    %% Organizing session length
    sessionID = sessionID(1:18);  % throws away excessive session for EPIC 10

    % for each session
    nRemovedT = nan(1,nSession);
    for j = 1:nSession
        % 0) Load data
        fileID = [dir '/' taskStrng '/run' num2str(sessionID(j))];
        load(fileID)
        allDataT = allData;

        % Create session-based list
        sess_list = array2table(zeros(nCond,length(header_output)),'VariableNames',header_output);
        sess_list.sess(:) = j;
        sess_list.n0cong = [1;2];

        % Current trial congruency
        matn0cong = nan(length(allDataT),1);
        for k = 1:length(allDataT)
            if strcmp(allDataT(k,1),'CONGRUENT')
                matn0cong(k) = 1;
            elseif strcmp(allDataT(k,1),'INCONGRUENT')
                matn0cong(k) = 2;
            else
                error('n0cong string error: unable to identify current trial congruency')
            end
        end

        % Create raw data table
        trialNum = (1:length(allDataT))';
        T = table(trialNum,matn0cong,double(cell2mat(allDataT(:,4))),double(cell2mat(allDataT(:,5))),...
            'VariableNames',["trialNum","n0cong","acc","rt"]);

        % Store data - for calculating grand mean
        RTdist{1,j,i} = T.rt(T.n0cong==1);  % congruent trials
        RTdist{2,j,i} = T.rt(T.n0cong==2);  % incongruent trials
        accdist{1,j,i} = T.acc(T.n0cong==1);
        accdist{2,j,i} = T.acc(T.n0cong==2);
        cRTdist{1,j,i} = T.rt(and(T.n0cong==1,T.acc==1));
        cRTdist{2,j,i} = T.rt(and(T.n0cong==2,T.acc==1));

        %% Calculate Congruency Effect
        % 1) RT mean
        % (1) Select correct trials
        acc_1_ids = find(T.acc==1);
        T_rm = T(acc_1_ids,:);
        nRemovedT(j) = size(T,1)-size(T_rm,1);

        % (2) Remove outlier trials
        % (2)-1 Set outlier criteria
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
        % (2)-2 Get outlier trial indices
        out_ids = [];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==1 & T_rm.rt<inbound(1))];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==2 & T_rm.rt<inbound(2))];

        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==1 & T_rm.rt>outbound(1))];
        out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==2 & T_rm.rt>outbound(2))];
        nRemovedT(j) = nRemovedT(j)+size(out_ids,1);
        % (2)-3 Mark outlier trials' accuracy as 99
        T_rm.acc(ismember(T_rm.trialNum,out_ids)) = 99;
        T.acc(ismember(T.trialNum,out_ids)) = 99;
        outRTdist{1,j,i} = T_rm.rt(and(T_rm.n0cong==1,T_rm.acc==99));
        outRTdist{2,j,i} = T_rm.rt(and(T_rm.n0cong==2,T_rm.acc==99));
        inRTdist{1,j,i} = T_rm.rt(and(T_rm.n0cong==1,T_rm.acc==1));
        inRTdist{2,j,i} = T_rm.rt(and(T_rm.n0cong==2,T_rm.acc==1));
        % (3) Calculate mean RTs (of correct and inbound trials)
        TGMrt = varfun(@mean,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        sess_list.RT = TGMrt.mean_rt;
        % 2) RT variance
        TGVrt = varfun(@var,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        sess_list.var = TGVrt.var_rt;
        % 3) Accuracy mean
        TGM_total = varfun(@mean,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGM_corr = varfun(@mean,T(T.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        sess_list.acc = TGM_corr.GroupCount./TGM_total.GroupCount;

        % 95% CI
        simCEsess = [];
        for si = 1:nSim
            tmpCdata = 1000*inRTdist{1,j,i};  % convert to msecs
            tmpIdata = 1000*inRTdist{2,j,i};  % convert to msecs
            simConsess = datasample(tmpCdata,length(tmpCdata));
            simIncsess = datasample(tmpIdata,length(tmpIdata));
            simCEsess = [simCEsess, mean(simIncsess)-mean(simConsess)];
        end

        %% Store data
        sess_cell{j,1} = table2array(sess_list);
        rtMat(1,j,i) = sess_list.RT(sess_list.n0cong==1);  % con
        rtMat(2,j,i) = sess_list.RT(sess_list.n0cong==2);  % inc
        accMat(1,j,i) = sess_list.acc(sess_list.n0cong==1);
        accMat(2,j,i) = sess_list.acc(sess_list.n0cong==2);
        seHMat(j,i) = (prctile(simCEsess,UpB)-prctile(simCEsess,LwB))/2;
    end
    num_removedT{i} = nRemovedT;

    %% Integrate session data
    Wdata = array2table(cell2mat(sess_cell),'VariableNames',header_output);

    %% Export to Excel
    savename = ['outputs/' taskStrng '_subj' num2str(subjID(i)) '_sess.xlsx'];
    writetable(Wdata, savename)
end
% 18 session means - for violin plots
CErtMat = squeeze(rtMat(2,:,:)-rtMat(1,:,:));  % inc-con
CEaccMat = squeeze(accMat(1,:,:)-accMat(2,:,:));  % con-inc
iesMat = (1000*rtMat)./accMat;
CEiesMat = squeeze(iesMat(2,:,:)-iesMat(1,:,:));

plot_rtMat = 1000*rtMat;  % convert to msecs
plot_CErtMat = 1000*CErtMat;  % convert to msecs

%% Output number of removed trials
tmpNR = sum([num_removedT{:}]);
disp([num2str(100*tmpNR/(nTrials*nSubj*nSession)) '% of participants removed for ' taskStrng ' task'])

%% Grand Mean and 95% Confidence Intervals
% RT
CErtGM = nan(nSubj,1);  % grand mean
CErtCH = nan(nSubj,1);  % confidence interval/half length
CrtGM = nan(nSubj,1);  % congruent
IrtGM = nan(nSubj,1);  % incongruent
CrtCH = nan(nSubj,1);
IrtCH = nan(nSubj,1);
fxdC = cell(nSubj,1);  % RTs after excluding outliers
fxdI = cell(nSubj,1);
for i = 1:nSubj
    % Remove outliers for RT data
    tmpC = [];
    tmpI = [];
    for s = 1:nSession
        tmpC = [tmpC; cRTdist{1,s,i}];
        tmpI = [tmpI; cRTdist{2,s,i}];
    end
    if isnan(tmpC)
        disp('NaN exists in tmpC')
    elseif isnan(tmpI)
        disp('NaN exists in tmpI')
    end
    tmpCgm = mean(tmpC,'omitnan');
    tmpIgm = mean(tmpI,'omitnan');
    tmpCsd = std(tmpC,'omitnan');
    tmpIsd = std(tmpI,'omitnan');
    lbound(1) = tmpCgm-nStd*tmpCsd;
    ubound(1) = tmpCgm+nStd*tmpCsd;
    lbound(2) = tmpIgm-nStd*tmpIsd;
    ubound(2) = tmpIgm+nStd*tmpIsd;
    fxdC{i} = tmpC(find(and(tmpC>=lbound(1),tmpC<=ubound(1))));
    fxdI{i} = tmpI(find(and(tmpI>=lbound(2),tmpI<=ubound(2))));
    CErtGM(i) = 1000*(mean(fxdI{i},'omitnan')-mean(fxdC{i},'omitnan'));  % convert to msec
    CrtGM(i) = 1000*mean(fxdC{i},'omitnan');
    IrtGM(i) = 1000*mean(fxdI{i},'omitnan');

    % 95% confidence interval by bootstrapping
    simCE = [];
    simC = [];
    simI = [];
    for si = 1:nSim
        simCon = datasample(fxdC{i},length(fxdC{i}));  % use tmpC if you don't want to remove outliers
        simInc = datasample(fxdI{i},length(fxdI{i}));
        simCE = [simCE, mean(simInc)-mean(simCon)];
        simC = [simC, mean(simCon)];
        simI = [simI, mean(simInc)];
    end
    CErtCH(i) = 1000*(prctile(simCE,UpB)-prctile(simCE,LwB))/2;
    CrtCH(i) = 1000*(prctile(simC,UpB)-prctile(simC,LwB))/2;
    IrtCH(i) = 1000*(prctile(simI,UpB)-prctile(simI,LwB))/2;
end

% Acc
CEaccGM = nan(nSubj,1);  % grand mean
CEaccCH = nan(nSubj,1);  % confidence interval/half length
CaccGM = nan(nSubj,1);  % congruent
IaccGM = nan(nSubj,1);  % incongruent
CaccCH = nan(nSubj,1);
IaccCH = nan(nSubj,1);
for i = 1:nSubj
    cAcc = [];
    iAcc = [];
    for s = 1:nSession
        cAcc = [cAcc; accdist{1,s,i}];
        iAcc = [iAcc; accdist{2,s,i}];
    end
    if isnan(cAcc)
        disp('NaN exists in cAcc')
    elseif isnan(iAcc)
        disp('NaN exists in iAcc')
    end
    CEaccGM(i) = mean(cAcc)-mean(iAcc);
    CaccGM(i) = mean(cAcc);
    IaccGM(i) = mean(iAcc);

    % 95% confidence interval by bootstrapping
    simCE = [];
    simC = [];
    simI = [];
    for si = 1:nSim
        simCon = datasample(cAcc,length(cAcc));
        simInc = datasample(iAcc,length(iAcc));
        simCE = [simCE, (mean(simCon)-mean(simInc))];
        simC = [simC, mean(simCon)];
        simI = [simI, mean(simInc)];
    end
    CEaccCH(i) = (prctile(simCE,UpB)-prctile(simCE,LwB))/2;
    CaccCH(i) = (prctile(simC,UpB)-prctile(simC,LwB))/2;
    IaccCH(i) = (prctile(simI,UpB)-prctile(simI,LwB))/2;
end

% IES
CiesGM = CrtGM./CaccGM;  % grand mean
IiesGM = IrtGM./IaccGM;
CEiesGM = IiesGM-CiesGM;
CEiesCH = nan(nSubj,1);  % confidence interval/half length
CiesCH = nan(nSubj,1);
IiesCH = nan(nSubj,1);
for i = 1:nSubj
    simCE = [];
    simC = [];
    simI = [];
    for si = 1:nSim
        simConRT = datasample(fxdC{i},length(fxdC{i}));
        simIncRT = datasample(fxdI{i},length(fxdI{i}));
        simConACC = datasample(cAcc,length(cAcc));
        simIncACC = datasample(iAcc,length(iAcc));
        simConIES = 1000*mean(simConRT)/mean(simConACC);
        simIncIES = 1000*mean(simIncRT)/mean(simIncACC);
        simCE = [simCE,simIncIES-simConIES];
        simC = [simC,simConIES];
        simI = [simI,simIncIES];
    end
    CEiesCH(i) = (prctile(simCE,UpB)-prctile(simCE,LwB))/2;
    CiesCH(i) = (prctile(simC,UpB)-prctile(simC,LwB))/2;
    IiesCH(i) = (prctile(simI,UpB)-prctile(simI,LwB))/2;
end

%% Save variables
% CE
PP_CErtMat = plot_CErtMat;  % RT 18 session means
PP_CEseMat = seHMat;
PP_CEaccMat = CEaccMat;  % Accuracy
PP_CEiesMat = CEiesMat;  % IES
PP_CErtGrandMean = CErtGM;  % Grand mean of all sessions
PP_CEaccGrandMean = CEaccGM;
PP_CEiesGrandMean = CEiesGM;
PP_CErtConHalf = CErtCH;  % 95% CI
PP_CEaccConHalf = CEaccCH;
PP_CEiesConHalf = CEiesCH;
% Congruent/Incongruent
PP_CrtMat = squeeze(plot_rtMat(1,:,:));
PP_IrtMat = squeeze(plot_rtMat(2,:,:));
PP_CrtGrandMean = CrtGM;
PP_IrtGrandMean = IrtGM;
PP_CrtConHalf = CrtCH;
PP_IrtConHalf = IrtCH;
PP_CaccMat = squeeze(accMat(1,:,:));
PP_IaccMat = squeeze(accMat(2,:,:));
PP_CaccGrandMean = CaccGM;
PP_IaccGrandMean = IaccGM;
PP_CaccConHalf = CaccCH;
PP_IaccConHalf = IaccCH;
PP_CiesMat = PP_CrtMat./PP_CaccMat;
PP_IiesMat = PP_IrtMat./PP_IaccMat;
PP_CiesGrandMean = CiesGM;
PP_IiesGrandMean = IiesGM;
PP_CiesConHalf = CiesCH;
PP_IiesConHalf = IiesCH;
save('PP_CEmat.mat','PP_CErtMat','PP_CEseMat','PP_CEaccMat','PP_CEiesMat',...
    'PP_CErtGrandMean','PP_CEaccGrandMean','PP_CEiesGrandMean',...
    'PP_CErtConHalf','PP_CEaccConHalf','PP_CEiesConHalf',...
    'PP_CrtMat','PP_IrtMat',...
    'PP_CrtGrandMean','PP_IrtGrandMean','PP_CrtConHalf','PP_IrtConHalf',...
    'PP_CaccMat','PP_IaccMat','PP_CaccGrandMean','PP_IaccGrandMean',...
    'PP_CaccConHalf','PP_IaccConHalf',...
    'PP_CiesMat','PP_IiesMat',...
    'PP_CiesGrandMean','PP_IiesGrandMean','PP_CiesConHalf','PP_IiesConHalf')