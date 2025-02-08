%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC - Flanker, Prime-probe, and Stroop Tasks
%
%% Instructions for running the script:
% You will need the EPIC dataset and the "session_numbering.xlsx" file,
% both available for download at https://osf.io/jk9nb.
%
%% Purpose of the script:
% This script examines within-subject variability in temporal effects,
% specifically the impairment in reaction time and accuracy over time.
%
%% Outputs:
%   - Dot plots of overall reaction time/accuracy across sessions (Supp.Fig.8)
%   - Bar graphs comparing the first and second halves of each session (Supp.Fig.9)
%
% Created on 12/18/2023 by HJ Lee
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
taskIndx=[1,3,4];  % 1:Flanker, 2:GoNogo, 3:PrimeProbe, 4:STROOP
taskStrng = {'Flanker','PrimeProbe','STROOP'};
nSession = 18;  % number of sessions for each task
nStd = 3;  % outlier criterion: 3 standard deviation from mean
nExclT = 50;  % number of trials in the beginning of each session to exclude

%% Preassignment
RTmat = nan(nSession,nTask,nSubj);  % overall RT
ACCmat = nan(nSession,nTask,nSubj);  % overall accuracy
hlfRTmat1 = nan(nSession,nTask,nSubj);  % first half
hlfRTmat2 = nan(nSession,nTask,nSubj);  % second half
hlfACCmat1 = nan(nSession,nTask,nSubj);
hlfACCmat2 = nan(nSession,nTask,nSubj);
for i = 1:nSubj
    dir = ['rawData/EPIC' num2str(subjID(i))];
    %% Load session index information
    opts = detectImportOptions('session_numbering.xlsx');
    opts.Sheet = ['subj' num2str(subjID(i))];
    sessionMat = readmatrix('session_numbering.xlsx',opts);
    for t = 1:nTask
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

            %% Overall RT/accuracy
            % (1) Select correct trials
            acc_1_ids = find(T.acc==1);
            T_rm = T(acc_1_ids,:);
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
            % (3) Calculate RTs (of correct and inbound trials) and accuracies
            RTmat(j,t,i) = mean(T_rm.rt(T_rm.acc==1))*1000;  % convert to msecs
            ACCmat(j,t,i) = mean(T.acc);

            %% Each halves
            tmpT_rm = T_rm(T_rm.acc==1,:);
            hlfRTmat1(j,t,i) = mean(tmpT_rm.rt(1+nExclT:floor(size(tmpT_rm,1)/2+nExclT/2)))*1000;  % exclude the first 50 trials
            hlfRTmat2(j,t,i) = mean(tmpT_rm.rt(floor(size(tmpT_rm,1)/2+nExclT/2)+1:end))*1000;
            tmpT_2_1 = table2array(T(1+nExclT:floor(size(T,1)/2+nExclT/2),3));
            tmpT_2_2 = table2array(T(floor(size(T,1)/2+nExclT/2)+1:end,3));
            hlfACCmat1(j,t,i) = mean(tmpT_2_1);
            hlfACCmat2(j,t,i) = mean(tmpT_2_2);
        end
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
    %legend('EPIC 03','EPIC 04','EPIC 05','EPIC 06',...
    %    'EPIC 07','EPIC 08','EPIC 10','EPIC 12','FontSize',15)
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
end

%% Comparing the first and second halves of each session
% Performance impairment within session (Supp.Fig.9)
% Reaction time
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
    if t == 1
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