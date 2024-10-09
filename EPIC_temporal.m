%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
% Impairment effects
%
% To run this script:
% You need the EPIC dataset and "session_numbering.xlsx"
% (download both at https://osf.io/jk9nb)
%
% What this script does:
% This script examines any changes in reaction time (RT) and accuracy that
% display performance impairment over time
%
% What this script outputs:
% Dot plots of overall RT/accuracy across sessions (Supp.Fig.8)
% Bar graphs of first vs. second halves of each session (Supp.Fig.9)
%
% Created on 12/18/2023 by HJ Lee
% Last modified on 01/24/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')

%% Parameter settings
% Participant information
tmpSubjID = 1:12;
subjExcld = [1,2,11,9];  % excluded participants
subjID = setdiff(tmpSubjID, subjExcld);
nSubj = length(subjID);

% Task information
nTask = 3;
taskIndx=[1,3,4];  % 1:Flanker, 2:GoNogo, 3:PrimeProbe, 4:STROOP
taskStrng = {'Flanker','PrimeProbe','STROOP'};
taskStrng_print = {'Flanker','Prime-Probe','Stroop'};  % for printing in figures
nBlocks = 4;
numTrials = [100; 96; 108];  % per block
nTrials = nBlocks*numTrials;  % total number of trials
nSession = 18;  % total number of sessions for each task
nCond = 2;  % congruent(1) vs. incongruent(2)
nStd = 3;  % outlier criterion; 3 standard deviation from mean
nExclT = 50;  % number of trials in the beginning of each session to exclude

% 95% confidence interval
UpB = 97.5;  % percentile; upper boundary of the confidence interval
LwB = 2.5;  % lower boundary
nSim = 1000;

% Table header
header_output1 = {'sess','RT','var','acc'};  % overall RT/accuracy
header_output2 = {'sess','n0cong','RT','var','acc'};  % CE RT/accuracy
header_output3 = {'hlf','RT','acc'};
header_output4 = {'sess','hlf','RT','acc'};  % overall RT/accuracy halves

%% Preassignment
RTmat = nan(nSession,nTask,nSubj);  % overall RT
ACCmat = nan(nSession,nTask,nSubj);  % overall accuracy
cRTmat = nan(nSession,nTask,nSubj);  % congruent trials RT
iRTmat = nan(nSession,nTask,nSubj);  % incongruent
cACCmat = nan(nSession,nTask,nSubj);  % congruent trials accuracy
iACCmat = nan(nSession,nTask,nSubj);  % incongruent
ceRTmat = nan(nSession,nTask,nSubj);  % CE RT
ceACCmat = nan(nSession,nTask,nSubj);  % CE accuracy
hlfRTmat1 = nan(nSession,nTask,nSubj);  % first half
hlfRTmat2 = nan(nSession,nTask,nSubj);  % second half
hlfACCmat1 = nan(nSession,nTask,nSubj);
hlfACCmat2 = nan(nSession,nTask,nSubj);

for i = 1:nSubj
    dir = ['Data_backups/IndivVar00' num2str(subjID(i))];  % rawData/EPIC
    % session information
    opts = detectImportOptions('session_numbering.xlsx');
    opts.Sheet = ['subj' num2str(subjID(i))];
    sessionMat = readmatrix('session_numbering.xlsx',opts);
    for t = 1:nTask
        %% Load session indexing information
        sessionID = sessionMat(:,taskIndx(t));
        sessionID = sessionID(1:nSession);
        for j = 1:nSession
            %% Load data
            fileID = [dir '/' taskStrng{t} '/run' num2str(sessionID(j))];
            load(fileID)
            %% Current trial congruency
            matn0cong = nan(length(allData),1);
            trialNum = (1:length(allData))';
            if t == 3  % Stroop
                for m = 1:length(allData)
                    if allData{m,3} == 0
                        matn0cong(m) = 3;  % netural
                    elseif allData{m,3} == 1
                        matn0cong(m) = 2;  % incongruent
                    elseif allData{m,3} == 2
                        matn0cong(m) = 1;  % congruent
                    else
                        error('n0cong string error: unable to identify current trial congruency')
                    end
                end
                T = table(trialNum,matn0cong,double(cell2mat(allData(:,6))),double(cell2mat(allData(:,5))),'VariableNames',["trialNum","n0cong","acc","rt"]);
            else  % Flanker & Prime-Probe
                for m = 1:length(allData)
                    if strcmp(allData(m,1),'CONGRUENT')
                        matn0cong(m) = 1;
                    elseif strcmp(allData(m,1),'INCONGRUENT')
                        matn0cong(m) = 2;
                    else
                        error('n0cong string error: unable to identify current trial congruency')
                    end
                end
                T = table(trialNum,matn0cong,double(cell2mat(allData(:,4))),double(cell2mat(allData(:,5))),'VariableNames',["trialNum","n0cong","acc","rt"]);
            end

            % Create session-based lists
            sess_list1 = array2table(zeros(1,length(header_output1)),'VariableNames',header_output1);
            sess_list1.sess = j;
            sess_list2 = array2table(zeros(nCond,length(header_output2)),'VariableNames',header_output2);
            sess_list2.sess(:) = j;
            sess_list2.n0cong = [1;2];
            sess_list3 = array2table(zeros(2,length(header_output3)),'VariableNames',header_output3);
            sess_list3.hlf = [1;2];
            sess_list4 = array2table(zeros(2,length(header_output4)),'VariableNames',header_output4);
            sess_list4.sess(:) = j;
            sess_list4.hlf = [1;2];

            %% Preprocess and store data
            % Overall RT/accuracy
            % (1) Select correct trials
            acc_1_ids = find(T.acc==1);
            T_rm = T(acc_1_ids,:);
            T_2 = T;
            % (2) Remove outlier trials
            Tgm = mean(T_rm.rt,'omitnan');
            Tgs = std(T_rm.rt,'omitnan');
            inbound = Tgm-Tgs*nStd;
            if inbound < 0
                inbound = 0;
            end
            outbound = Tgm+Tgs*nStd;
            out_ids = [];
            out_ids = [out_ids; T_rm.trialNum(T_rm.rt<inbound)];
            out_ids = [out_ids; T_rm.trialNum(T_rm.rt>outbound)];
            T_rm.acc(ismember(T_rm.trialNum,out_ids)) = 99;
            T_2.acc(ismember(T_2.trialNum,out_ids)) = 99;
            % (3) Caclulate RTs (of correct and inbound trials) and accuracies
            RTmat(j,t,i) = mean(T_rm.rt(T_rm.acc==1))*1000;  % convert to msecs
            ACCmat(j,t,i) = length(find(T_2.acc==1)) / length(find(T_2.acc~=99));
            sess_list1.RT = RTmat(j,t,i);
            sess_list1.var = var(T_rm.rt(T_rm.acc==1))*1000;  % convert to msecs
            sess_list1.acc = ACCmat(j,t,i);

            % Each halves
            tmpT_rm = T_rm(T_rm.acc==1,:);
            hlfRTmat1(j,t,i) = mean(tmpT_rm.rt(1+nExclT:floor(size(tmpT_rm,1)/2+nExclT/2)))*1000;
            hlfRTmat2(j,t,i) = mean(tmpT_rm.rt(floor(size(tmpT_rm,1)/2+nExclT/2)+1:end))*1000;
            tmpT_2_1 = table2array(T_2(1+nExclT:floor(size(T_2,1)/2+nExclT/2),3));
            tmpT_2_2 = table2array(T_2(floor(size(T_2,1)/2+nExclT/2)+1:end,3));
            hlfACCmat1(j,t,i) = length(find(tmpT_2_1==1)) / length(find(tmpT_2_1~=99));
            hlfACCmat2(j,t,i) = length(find(tmpT_2_2==1)) / length(find(tmpT_2_2~=99));
            sess_list4.RT(1) = hlfRTmat1(j,t,i);
            sess_list4.RT(2) = hlfRTmat2(j,t,i);
            sess_list4.acc(1) = hlfACCmat1(j,t,i);
            sess_list4.acc(2) = hlfACCmat2(j,t,i);
            clear T_rm inbound outbound out_ids

            % CE RT/accuracy
            % (1) Select correct trials
            T_rm = T(acc_1_ids,:);
            % (2) Remove outlier trials
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
            out_ids = [];
            out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==1 & T_rm.rt<inbound(1))];
            out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==2 & T_rm.rt<inbound(2))];

            out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==1 & T_rm.rt>outbound(1))];
            out_ids = [out_ids; T_rm.trialNum(T_rm.n0cong==2 & T_rm.rt>outbound(2))];
            T_rm.acc(ismember(T_rm.trialNum,out_ids)) = 99;
            T.acc(ismember(T.trialNum,out_ids)) = 99;
            % (3) Caclulate RTs (of correct and inbound trials) and accuracies
            TGMrt = varfun(@mean,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
            sess_list2.RT = TGMrt.mean_rt(1:2)*1000;  % convert to msecs
            TGVrt = varfun(@var,T_rm(T_rm.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
            sess_list2.var = TGVrt.var_rt(1:2)*1000;  % convert to msecs
            TGM_total = varfun(@mean,T(T.acc~=99,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
            TGM_corr = varfun(@mean,T(T.acc==1,:),'GroupingVariables',{'n0cong'},'OutputFormat','table');
            sess_list2.acc = TGM_corr.GroupCount(1:2) ./ TGM_total.GroupCount(1:2);
            cRTmat(j,t,i) = sess_list2.RT(sess_list2.n0cong==1);
            iRTmat(j,t,i) = sess_list2.RT(sess_list2.n0cong==2);
            cACCmat(j,t,i) = sess_list2.acc(sess_list2.n0cong==1);
            iACCmat(j,t,i) = sess_list2.acc(sess_list2.n0cong==2);
            ceRTmat(j,t,i) = sess_list2.RT(sess_list2.n0cong==2)-sess_list2.RT(sess_list2.n0cong==1);
            ceACCmat(j,t,i) = sess_list2.acc(sess_list2.n0cong==1)-sess_list2.acc(sess_list2.n0cong==2);
        end
        sess_list3.RT(1) = mean(hlfRTmat1(:,t,i));
        sess_list3.RT(2) = mean(hlfRTmat2(:,t,i));
        sess_list3.acc(1) = mean(hlfACCmat1(:,t,i));
        sess_list3.acc(2) = mean(hlfACCmat2(:,t,i));
    end
end

%% Connected dot plots
% Performance impairment across sessions (Supp.Fig.8)
i_cmap = turbo(nSubj);
for t = 1:nTask
    % RT
    figure
    for i = 1:nSubj
        plot(1:nSession,RTmat(:,t,i),'Marker','o','Color',i_cmap(i,:),'LineWidth',2); hold on
    end
    set(gca,'FontSize',18)
    grid on
    xlabel('Session','FontSize',25)
    xticks(1:nSession)
    ylabel('Reaction time (ms)','FontSize',25)
    if t == 1
        ylim([300 650])
    elseif t == 2
        ylim([500 850])
    elseif t == 3
        ylim([400 850])
    end
    legend('Participant 03','Participant 04','Participant 05','Participant 06',...
        'Participant 07','Participant 08','Participant 10','Participant 12','FontSize',11)
    % Accuracy
    figure
    for i = 1:nSubj
        plot(1:nSession,ACCmat(:,t,i),'Marker','o','Color',i_cmap(i,:),'LineWidth',2); hold on
    end
    set(gca,'FontSize',18)
    grid on
    xlabel('Session','FontSize',25)
    xticks(1:nSession)
    ylabel('Accuracy','FontSize',25)
    ylim([0.88 1])
    %legend('Participant 03','Participant 04','Participant 05','Participant 06',...
    %    'Participant 07','Participant 08','Participant 10','Participant 12','FontSize',15)
end

%% Comparing the first and second halves of each session
% Performance impairment within session (Supp.Fig.9)
% RT
for t = 1:nTask
    figure
    i_mSess1 = squeeze(mean(hlfRTmat1(:,t,:)));
    i_mSess2 = squeeze(mean(hlfRTmat2(:,t,:)));
    i_mSess = [i_mSess1,i_mSess2];
    g_mSess1 = mean(i_mSess1);
    g_mSess2 = mean(i_mSess2);
    x = categorical({'First half','Second half'});
    y = [g_mSess1,g_mSess2];
    b = bar(x,y,'white'); hold on
    b.FaceColor = [0.5 0.5 0.5];
    set(gca,'FontSize',20)
    ylabel('Reaction time (ms)','FontSize',25)
    for i = 1:nSubj
        plot(1:2,i_mSess(i,:),'Marker','o','Color',i_cmap(i,:),'LineWidth',2)
    end
    %scatter(1,i_mSess1,80,i_cmap,'filled','jitter','on','jitterAmount',0.3)
    %scatter(2,i_mSess2,80,i_cmap,'filled','jitter','on','jitterAmount',0.3)
    %legend('','Participant 03','Participant 04','Participant 05','Participant 06',...
    %    'Participant 07','Participant 08','Participant 10','Participant 12')
    if t == 1
        %legend('','03','04','05','06','07','08','10','12')
        ylim([300 510])
    elseif t == 2
        ylim([500 710])
    elseif t == 3
        ylim([500 710])
    end
end
% Accuracy
for t = 1:nTask
    figure
    i_mSess1 = squeeze(mean(hlfACCmat1(:,t,:)));
    i_mSess2 = squeeze(mean(hlfACCmat2(:,t,:)));
    i_mSess = [i_mSess1,i_mSess2];
    g_mSess1 = mean(i_mSess1);
    g_mSess2 = mean(i_mSess2);
    x = categorical({'First half','Second half'});
    y = [g_mSess1,g_mSess2];
    b = bar(x,y,'white'); hold on
    b.FaceColor = [0.5 0.5 0.5];
    set(gca,'FontSize',20)
    ylabel('Accuracy','FontSize',25)
    for i = 1:nSubj
        plot(1:2,i_mSess(i,:),'Marker','o','Color',i_cmap(i,:),'LineWidth',2)
    end
    ylim([0.9 1])
end