%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC - Flanker task analysis
%
%% Instructions for running the script:
% You will need the EPIC dataset and the "session_numbering.xlsx" file,
% both available for download at https://osf.io/jk9nb.
%
%% Purpose of the script:
% This script performs preprocessing on the flanker task data.
%   - Independent Variables:
%       Session: Session number (1-18).
%       Current trial congruency: 1 for congruent and 2 for incongruent.
%   - Dependent Variables:
%       Reaction time (RT)
%       Accuracy (acc)
%       Inverse Efficiency Score (IES)
%
%% Outputs of the script:
%   - Excel file: Ouptuts the mean and variance in an Excel file
%   - MAT file: Outputs a .mat file, which can be used in the
%   EPIC_violinPlotGrandmeanPlot.m script.
%
% Created on 10/02/2022 by HJ Lee
% Last modified on 02/03/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')

%% Experimental and statistical parameters
% Participant information
tmpSubjID = 1:12;
subjExcld = [1,2,11,9];  % EPIC 9 excluded
subjID = setdiff(tmpSubjID, subjExcld);
nSubj = length(subjID);

% Task information
nTask = 4;
taskIndx = 1;  % 1:Flanker, 2:GoNogo, 3:PrimeProbe, 4:STROOP
taskStrng = 'Flanker';
nBlocks = 4;
numTrials = 100;  % per block
nTrials = nBlocks*numTrials;  % total number of trials per session
nSession = 18;  % number of sessions for each task
nCond = 2;  % congruent(1) vs. incongruent(2)
nStd = 3;  % outlier criterion: 3 standard deviation from mean

% 95% confidence interval
UpB = 97.5;  % Upper percentile for 95% CI
LwB = 2.5;  % Lower percentile for 95% CI
nSim = 1000;

% Table header
header_output = {'sess','n0cong','RT','var','acc'};

%% Preassigment
% Matrices for mean RT/acc across conditions, sessions, and participants
rtMat = nan(nCond,nSession,nSubj);
accMat = nan(nCond,nSession,nSubj);
seHMat = nan(nSession,nSubj);  % standard error 18 sessions (95% CI half)
accdist = cell(nCond,nSession,nSubj);  % acc
cRTdist = cell(nCond,nSession,nSubj);  % correctly responded RTs
inRTdist = cell(nCond,nSession,nSubj);  % RT after removing incorrect trials and outliers
num_removedT = cell(nSubj,1);  % Track numbers of excluded trials
% for each participant
for i = 1:nSubj
    dir = ['rawData/EPIC' num2str(subjID(i))];
    %% Load session index information
    opts = detectImportOptions('session_numbering.xlsx');
    opts.Sheet = ['subj' num2str(subjID(i))];
    tmpMat = readmatrix('session_numbering.xlsx',opts);
    sessionID = tmpMat(:,taskIndx);
    sessionID = sessionID(1:nSession);  % throw away excessive session for EPIC 10

    % for each session
    nRemovedT = nan(1,nSession);
    for j = 1:nSession
        % 0) Load data
        fileID = [dir '/' taskStrng '/run' num2str(sessionID(j))];
        load(fileID)
        allDataT = allData;
        if length(allData) ~= nTrials
            error(['Check allData length for participant ' num2str(subjID(i))])
        end

        % Create a session-level list
        sess_list = array2table(zeros(nCond,length(header_output)),'VariableNames',header_output);
        sess_list.sess(:) = j;
        sess_list.n0cong = [1;2];

        % Current trial congruency (n0cong)
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
        T = table(trialNum,matn0cong,double(cell2mat(allDataT(:,4))),double(cell2mat(allDataT(:,5))),'VariableNames',["trialNum","n0cong","acc","rt"]);

        % Store data - for calculating the grand mean
        accdist{1,j,i} = T.acc(T.n0cong==1);
        accdist{2,j,i} = T.acc(T.n0cong==2);
        cRTdist{1,j,i} = T.rt(and(T.n0cong==1,T.acc==1));
        cRTdist{2,j,i} = T.rt(and(T.n0cong==2,T.acc==1));

        %% Calculate Congruency Effect
        % 1) RT mean
        % (1) Select correct trials
        acc_1_ids = find(T.acc==1);
        T_rm = T(acc_1_ids,:);
        nRemovedT(j) = size(T,1)-size(T_rm,1);  % # of removed incorrect trials

        % (2) Remove outlier trials
        % (2)-1 Set outlier criterion
        TGM = varfun(@mean,T_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');  % mean
        TGS = varfun(@std,T_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');  % standard deviation
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
        nRemovedT(j) = nRemovedT(j)+size(out_ids,1);  % # of outliers
        % (2)-3 Mark outlier trials' accuracy as 99
        T_rm.acc(ismember(T_rm.trialNum,out_ids)) = 99;
        %T.acc(ismember(T.trialNum,out_ids)) = 99;
        inRTdist{1,j,i} = T_rm.rt(and(T_rm.n0cong==1,T_rm.acc==1));
        inRTdist{2,j,i} = T_rm.rt(and(T_rm.n0cong==2,T_rm.acc==1));
        % (3) Calculate mean RTs (of correct and inbound trials)
        TGMrt = varfun(@mean,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        sess_list.RT = TGMrt.mean_rt;
        % 2) RT variance
        TGVrt = varfun(@var,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        sess_list.var = TGVrt.var_rt;
        % 3) Accuracy mean
        TGM_total = varfun(@mean,T,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGM_corr = varfun(@mean,T(T.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
        sess_list.acc = TGM_corr.GroupCount./TGM_total.GroupCount;

        % 95% CI of the congruency effect
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

%% Display the number of removed trials
tmpNR = sum([num_removedT{:}]);
disp([num2str(100*tmpNR/(nTrials*nSubj*nSession)) '% of participants removed for ' taskStrng ' task'])

%% Grand Mean and 95% Confidence Intervals
% RT
rtGM = nan(nSubj,1);  % overall grand mean
rtCH = nan(nSubj,1);  % overal confidence interval/half length
CErtGM = nan(nSubj,1);  % grand mean
CErtCH = nan(nSubj,1);  % confidence interval/half length
CrtGM = nan(nSubj,1);  % congruent
IrtGM = nan(nSubj,1);  % incongruent
CrtCH = nan(nSubj,1);
IrtCH = nan(nSubj,1);
fxdC = cell(nSubj,1);  % RTs after excluding outliers
fxdI = cell(nSubj,1);
fxdA = cell(nSubj,1);  % all trials RT regardless of trial conditions
for i = 1:nSubj
    % Remove outliers for RT data
    tmpC = [];
    tmpI = [];
    tmpA = [];
    for s = 1:nSession
        tmpC = [tmpC; cRTdist{1,s,i}];
        tmpI = [tmpI; cRTdist{2,s,i}];
        tmpA = [tmpA; cRTdist{1,s,i}; cRTdist{2,s,i}];
    end
    if any(isnan(tmpC))
        disp('NaN exists in tmpC')
    elseif any(isnan(tmpI))
        disp('NaN exists in tmpI')
    end
    tmpCgm = mean(tmpC,'omitnan');
    tmpIgm = mean(tmpI,'omitnan');
    tmpAgm = mean(tmpA,'omitnan');
    tmpCsd = std(tmpC,'omitnan');
    tmpIsd = std(tmpI,'omitnan');
    tmpAsd = std(tmpA,'omitnan');
    lbound(1) = tmpCgm-nStd*tmpCsd;
    ubound(1) = tmpCgm+nStd*tmpCsd;
    lbound(2) = tmpIgm-nStd*tmpIsd;
    ubound(2) = tmpIgm+nStd*tmpIsd;
    lbound(3) = tmpAgm-nStd*tmpAsd;
    ubound(3) = tmpAgm+nStd*tmpAsd;
    fxdC{i} = tmpC(tmpC>=lbound(1) & tmpC<=ubound(1));
    fxdI{i} = tmpI(tmpI>=lbound(2) & tmpI<=ubound(2));
    fxdA{i} = tmpA(tmpA>=lbound(3) & tmpA<=ubound(3));
    rtGM(i) = 1000*mean(fxdA{i},'omitnan');
    CErtGM(i) = 1000*(mean(fxdI{i},'omitnan')-mean(fxdC{i},'omitnan'));  % convert to msec
    CrtGM(i) = 1000*mean(fxdC{i},'omitnan');
    IrtGM(i) = 1000*mean(fxdI{i},'omitnan');

    % 95% confidence interval by bootstrapping
    simCE = [];
    simC = [];
    simI = [];
    simA = [];
    for si = 1:nSim
        simCon = datasample(fxdC{i},length(fxdC{i}));  % use tmpC if you don't want to remove outliers
        simInc = datasample(fxdI{i},length(fxdI{i}));
        simAll = datasample(fxdA{i},length(fxdA{i}));
        simCE = [simCE, mean(simInc)-mean(simCon)];
        simC = [simC, mean(simCon)];
        simI = [simI, mean(simInc)];
        simA = [simA, mean(simAll)];
    end
    rtCH(i) = 1000*(prctile(simA,UpB)-prctile(simA,LwB))/2;
    CErtCH(i) = 1000*(prctile(simCE,UpB)-prctile(simCE,LwB))/2;
    CrtCH(i) = 1000*(prctile(simC,UpB)-prctile(simC,LwB))/2;
    IrtCH(i) = 1000*(prctile(simI,UpB)-prctile(simI,LwB))/2;
end

% Acc
accGM = nan(nSubj,1);  % overall grand mean
accCH = nan(nSubj,1);  % overall confidence interval/half length
CEaccGM = nan(nSubj,1);  % grand mean
CEaccCH = nan(nSubj,1);  % confidence interval/half length
CaccGM = nan(nSubj,1);  % congruent
IaccGM = nan(nSubj,1);  % incongruent
CaccCH = nan(nSubj,1);
IaccCH = nan(nSubj,1);
fxdaccC = cell(nSubj,1);
fxdaccI = cell(nSubj,1);
fxdaccA = cell(nSubj,1);
for i = 1:nSubj
    cAcc = [];
    iAcc = [];
    aAcc = [];
    for s = 1:nSession
        cAcc = [cAcc; accdist{1,s,i}];
        iAcc = [iAcc; accdist{2,s,i}];
        aAcc = [aAcc; accdist{1,s,i}; accdist{2,s,i}];
    end
    if isnan(cAcc)
        disp('NaN exists in cAcc')
    elseif isnan(iAcc)
        disp('NaN exists in iAcc')
    end
    fxdaccC{i} = cAcc;
    fxdaccI{i} = iAcc;
    fxdaccA{i} = aAcc;
    accGM(i) = mean(aAcc);
    CEaccGM(i) = mean(cAcc)-mean(iAcc);
    CaccGM(i) = mean(cAcc);
    IaccGM(i) = mean(iAcc);

    % 95% confidence interval by bootstrapping
    simCE = [];
    simC = [];
    simI = [];
    simA = [];
    for si = 1:nSim
        simCon = datasample(cAcc,length(cAcc));
        simInc = datasample(iAcc,length(iAcc));
        simAll = datasample(aAcc,length(aAcc));
        simCE = [simCE, (mean(simCon)-mean(simInc))];
        simC = [simC, mean(simCon)];
        simI = [simI, mean(simInc)];
        simA = [simA, mean(simAll)];
    end
    accCH(i) = (prctile(simA,UpB)-prctile(simA,LwB))/2;
    CEaccCH(i) = (prctile(simCE,UpB)-prctile(simCE,LwB))/2;
    CaccCH(i) = (prctile(simC,UpB)-prctile(simC,LwB))/2;
    IaccCH(i) = (prctile(simI,UpB)-prctile(simI,LwB))/2;
end

% IES
iesGM = rtGM./accGM;  % overall grand mean
iesCH = nan(nSubj,1);  % overall confidence interval/half length
CiesGM = CrtGM./CaccGM;  % grand mean
IiesGM = IrtGM./IaccGM;
CEiesGM = IiesGM-CiesGM;
CEiesCH = nan(nSubj,1);  % confidence interval/half length
CiesCH = nan(nSubj,1);
IiesCH = nan(nSubj,1);
% 95% confidence interval by bootstrapping
for i = 1:nSubj
    simCE = [];
    simC = [];
    simI = [];
    simA = [];
    for si = 1:nSim
        simConRT = datasample(fxdC{i},length(fxdC{i}));  % tmpC
        simIncRT = datasample(fxdI{i},length(fxdI{i}));
        simRT = datasample(fxdA{i},length(fxdA{i}));
        simConACC = datasample(fxdaccC{i},length(fxdaccC{i}));
        simIncACC = datasample(fxdaccI{i},length(fxdaccI{i}));
        simACC = datasample(fxdaccA{i},length(fxdaccA{i}));
        simConIES = 1000*mean(simConRT)/mean(simConACC);
        simIncIES = 1000*mean(simIncRT)/mean(simIncACC);
        simIES = 1000*mean(simRT)/mean(simACC);
        simCE = [simCE,simIncIES-simConIES];
        simC = [simC,simConIES];
        simI = [simI,simIncIES];
        simA = [simA,simIES];
    end
    iesCH(i) = (prctile(simA,UpB)-prctile(simA,LwB))/2;
    CEiesCH(i) = (prctile(simCE,UpB)-prctile(simCE,LwB))/2;
    CiesCH(i) = (prctile(simC,UpB)-prctile(simC,LwB))/2;
    IiesCH(i) = (prctile(simI,UpB)-prctile(simI,LwB))/2;
end

%% Save variables to draw violin plots and error bar plots across tasks
% Overall
FL_rtGrandMean = rtGM;
FL_rtConHalf = rtCH;
FL_accGrandMean = accGM;
FL_accConHalf = accCH;
FL_iesGrandMean = iesGM;
FL_iesConHalf = iesCH;
% CE
FL_CErtMat = plot_CErtMat;  % RT 18 session means
FL_CEseMat = seHMat;
FL_CEaccMat = CEaccMat;  % Accuracy
FL_CEiesMat = CEiesMat;  % IES
FL_CErtGrandMean = CErtGM;  % Grand mean of all sessions
FL_CEaccGrandMean = CEaccGM;
FL_CEiesGrandMean = CEiesGM;
FL_CErtConHalf = CErtCH;  % 95% CI
FL_CEaccConHalf = CEaccCH;
FL_CEiesConHalf = CEiesCH;
% Congruent/Incongruent
FL_CrtMat = squeeze(plot_rtMat(1,:,:));
FL_IrtMat = squeeze(plot_rtMat(2,:,:));
FL_CrtGrandMean = CrtGM;
FL_IrtGrandMean = IrtGM;
FL_CrtConHalf = CrtCH;
FL_IrtConHalf = IrtCH;
FL_CaccMat = squeeze(accMat(1,:,:));
FL_IaccMat = squeeze(accMat(2,:,:));
FL_CaccGrandMean = CaccGM;
FL_IaccGrandMean = IaccGM;
FL_CaccConHalf = CaccCH;
FL_IaccConHalf = IaccCH;
FL_CiesMat = FL_CrtMat./FL_CaccMat;
FL_IiesMat = FL_IrtMat./FL_IaccMat;
FL_CiesGrandMean = CiesGM;
FL_IiesGrandMean = IiesGM;
FL_CiesConHalf = CiesCH;
FL_IiesConHalf = IiesCH;
save('FL_CEmat.mat', 'FL_rtGrandMean','FL_accGrandMean','FL_iesGrandMean',...
    'FL_rtConHalf','FL_accConHalf','FL_iesConHalf',...
    'FL_CErtMat','FL_CEseMat','FL_CEaccMat','FL_CEiesMat',...
    'FL_CErtGrandMean','FL_CEaccGrandMean','FL_CEiesGrandMean',...
    'FL_CErtConHalf','FL_CEaccConHalf','FL_CEiesConHalf',...
    'FL_CrtMat','FL_IrtMat',...
    'FL_CrtGrandMean','FL_IrtGrandMean','FL_CrtConHalf','FL_IrtConHalf',...
    'FL_CaccMat','FL_IaccMat','FL_CaccGrandMean','FL_IaccGrandMean',...
    'FL_CaccConHalf','FL_IaccConHalf',...
    'FL_CiesMat','FL_IiesMat',...
    'FL_CiesGrandMean','FL_IiesGrandMean','FL_CiesConHalf','FL_IiesConHalf')