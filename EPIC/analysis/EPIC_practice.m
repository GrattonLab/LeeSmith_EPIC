%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC - Flanker, Prime-probe, and Stroop Tasks
%
%% Instructions for running the script:
% You will need the EPIC dataset and the "session_numbering.xlsx" file,
% both available for download at https://osf.io/jk9nb.
%
%% Purpose of the script:
% This script examines within-subject variability in temporal effects,
% specifically the improvement in congruency effect over time.
%
%% Outputs:
% Line graphs of the congruency effect over time (Supp.Fig.7)
%
% Created on 10/02/2022 by HJ Lee
% Last modified on 02/03/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experimental and statistical parameters
% Participant information
tmpSubjID = 1:12;
subjExcld = [1,2,11,9];  % EPIC 9 excluded
subjID = setdiff(tmpSubjID, subjExcld);
nSubj = length(subjID);

% Task information
nTask = 3;
taskIndx = [1,3,4];  % 1:Flanker, 2:GoNogo, 3:PrimeProbe, 4:STROOP
taskStrng = {'Flanker','PrimeProbe','STROOP'};
nTrials = [100,96,108]*4;
nSession = 18;  % number of sessions for each task
nSessionperSet = 3;  % concatenate 3 consecutive sessions to have reliable results ~1000 trials
nSet = nSession/nSessionperSet;
setID = 1:nSet;
nCond = 2;  % congruent(1) vs. incongruent(2)
nStd = 3;  % outlier criterion: 3 standard deviation from mean

% Table header
header_output1 = {'set','n0cong','RT','var','acc'};

% Colormap for plots
cmap = turbo(nSubj);

%% Preassignment
rtMat = nan(nCond,nSet,nSubj);  % mean RT across conditions, sets, and participants
accMat = nan(nCond,nSet,nSubj);
for t = 1:nTask
    for i = 1:nSubj
        dir = ['rawData/EPIC' num2str(subjID(i))];
        %% Load session index information
        opts = detectImportOptions('session_numbering.xlsx');
        opts.Sheet = ['subj' num2str(subjID(i))];
        tmpMat = readmatrix('session_numbering.xlsx',opts);
        sessionID = tmpMat(:,taskIndx(t));
        sessionID = sessionID(1:nSession);
        for j = 1:nSet
            % 0) Load data
            if ~isnan(sessionID(3*(setID(j)-1)+1))  % skip data loading if sessionID is nan
                fileID1 = [dir '/' taskStrng{t} '/run' num2str(sessionID(3*(setID(j)-1)+1))];
                load(fileID1)
                allData1 = allData;
            else
                allData1 = [];
            end
            if ~isnan(sessionID(3*(setID(j)-1)+2))
                fileID2 = [dir '/' taskStrng{t} '/run' num2str(sessionID(3*(setID(j)-1)+2))];
                load(fileID2)
                allData2 = allData;
            else
                allData2 = [];
            end
            if ~isnan(sessionID(3*(setID(j)-1)+3))
                fileID3 = [dir '/' taskStrng{t} '/run' num2str(sessionID(3*(setID(j)-1)+3))];
                load(fileID3)
                allData3 = allData;
            else
                allData3 = [];
            end
            allDataT = [allData1;allData2;allData3];
            if length(allDataT) ~= (nTrials(t)*nSessionperSet)
                error(['Check allDataT length for participant ' num2str(subjID(i))])
            end
            % Create a set-level list
            set_list = array2table(zeros(nCond,length(header_output1)),'VariableNames',header_output1);
            set_list.set(:) = j;
            set_list.n0cong = [1;2];
            trialNum = (1:length(allDataT))';

            % Current trial congruency
            matn0cong = nan(length(allDataT),1);
            if t == 3  % Stroop
                for m = 1:length(allDataT)
                    if allDataT{m,3} == 0
                        matn0cong(m) = 3;  % netural
                    elseif allDataT{m,3} == 1
                        matn0cong(m) = 2;  % incongruent
                    elseif allDataT{m,3} == 2
                        matn0cong(m) = 1;  % congruent
                    else
                        error('n0cong string error: unable to identify current trial congruency')
                    end
                end
                T = table(trialNum,matn0cong,double(cell2mat(allDataT(:,6))),double(cell2mat(allDataT(:,5))),'VariableNames',["trialNum","n0cong","acc","rt"]);
            else  % Flanker & Prime-Probe
                for m = 1:length(allDataT)
                    if strcmp(allDataT(m,1),'CONGRUENT')
                        matn0cong(m) = 1;
                    elseif strcmp(allDataT(m,1),'INCONGRUENT')
                        matn0cong(m) = 2;
                    else
                        error('n0cong string error: unable to identify current trial congruency')
                    end
                end
                T = table(trialNum,matn0cong,double(cell2mat(allDataT(:,4))),double(cell2mat(allDataT(:,5))),'VariableNames',["trialNum","n0cong","acc","rt"]);
            end

            %% Calculate Congruency Effect
            % 1) RT mean
            % (1) Select correct trials
            acc_1_ids = find(T.acc==1);
            T_rm = T(acc_1_ids,:);
            % (2) Remove outlier trials
            % (2)-1 Set outlier criterion
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
            % (2)-2 Get outlier trial indices
            out_ids = [];
            out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==1 & T_rm.rt<inbound(1))];
            out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==2 & T_rm.rt<inbound(2))];
            out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==1 & T_rm.rt>outbound(1))];
            out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==2 & T_rm.rt>outbound(2))];
            % (2)-3 Mark outlier trials' accuracy as 99
            T_rm.acc(ismember(T_rm.trialNum,out_ids)) = 99;
            %T.acc(ismember(T.trialNum,out_ids)) = 99;
            % (3) Calculate mean RTs (of correct and inbound trials)
            TGMrt = varfun(@mean,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
            set_list.RT = TGMrt.mean_rt(1:2);
            % 2) RT variance
            TGVrt = varfun(@var,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
            set_list.var = TGVrt.var_rt(1:2);
            % 3) Accuracy mean
            TGM_total = varfun(@mean,T,'GroupingVariables',{'n0cong'},'OutputFormat','table');
            TGM_corr = varfun(@mean,T(T.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
            tmp = TGM_corr.GroupCount./TGM_total.GroupCount;
            set_list.acc = tmp(1:2);
            %% Store data
            rtMat(1,j,i) = set_list.RT(set_list.n0cong==1);  % cong
            rtMat(2,j,i) = set_list.RT(set_list.n0cong==2);  % incong
            accMat(1,j,i) = set_list.acc(set_list.n0cong==1);
            accMat(2,j,i) = set_list.acc(set_list.n0cong==2);
        end

        %% Integrate set data
        CErtMat = squeeze(rtMat(2,:,:)-rtMat(1,:,:));  % incong-cong
        CEaccMat = squeeze(accMat(1,:,:)-accMat(2,:,:));  % cong-incong
    end
    %% Plot Congruency Effect across 6 Sets
    plot_CErtMat = 1000*CErtMat;  % convert to msecs
    % RT
    figure
    for l = 1:nSubj
        plot(setID,plot_CErtMat(:,l),'Color',cmap(l,:),'LineWidth',2.8); hold on
    end
    set(gca,'FontSize',18)
    xlabel('Three consecutive sessions','FontSize',18)
    xticks(1:nSet)
    ylim([0 150])  % [-20 200]
    %legend('EPIC 03','EPIC 04','EPIC 05','EPIC 06',...
    %    'EPIC 07','EPIC 08','EPIC 10','EPIC 12','FontSize',14,'Location','best')

    % Accuracy
    figure
    for l = 1:nSubj
        plot(1:nSet,CEaccMat(:,l),'Color',cmap(l,:),'LineWidth',2.8); hold on
    end
    set(gca,'FontSize',18)
    xlabel('Three consecutive sessions','FontSize',18)
    xticks(1:nSet)
    ylim([-0.05 0.12])
end