%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC - Public Data 2
%
%% Instructions for running the script:
% You will need the following datasets and file.
%   - EPIC flanker task data (https://osf.io/jk9nb)
%   - Excel file: "session_numbering.xlsx" (https://osf.io/jk9nb)
%   - Hedge et al.'s (2018) flanker task data (https://osf.io/cwzds/)
%
%% Purpose of the script:
% This script runs three simulations with a fixed between-subject mean and 
% standard deviation, while varying within-subject standard deviation to 
% demonstrate how increasing within-subject variability affects 
% between-subject variability.
% The within-subject distribution parameters are based on EPIC data, and
% the between-subject distribution parameters are derived from Hedge et 
% al.'s data.
% Simulation 1. Fixed between-subject mean and standard deviation
% Simulation 2. Varying between-subject standard deviation (5ms 10 20 30...)
% Simulation 3. Different numbers of simulated subjects (100 200 300 400 500...)
%
%% Outputs:
% Figure 6.
%
% Created on 12/21/2023 by HJ Lee
% Last modified on 02/04/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')  % reset the seed of the random number generator

%% Parameter settings
% For simulations
n_subs = 100;  % number of simulated subjects
n_subs_list = [50 100 200 300 400 500 1000];  % to loop through different subject number (Sim.3)
bs_std_list = [5 10 20 30 40 50 60];  % to loop through different true between-subject standard deviation (Sim.2)
bs_dist_n = 5000;  % number of data points of between-subject distribution
ws_dist_n = 100;  % number of data points of within-subject distribution
numTest = 1000;  % repetitions to get mean and standard error of simulation results
UpB = 97.5;  % upper percentile of 95% confidence interval
LwB = 2.5;  % lower percentile
nStd = 3;  % outlier exclusion criterion

%% Within-subject distribution parameters - EPIC flanker task data
tmpSubjID = 1:12;
subjExcld = [1,2,11,9];  % EPIC 9 excluded
subjID = setdiff(tmpSubjID, subjExcld);
nSubj = length(subjID);
taskIndx = 1;  % flanker
taskStrng = 'flanker';
nSession = 18;
nCond = 2;  % congruent (1), incongruent (2)
nSampling = [50 100 200 400 800 1600 3200 6400]./2;  % number of trial sampling; to calculate within-subject variability
fxdC = cell(nSubj,1);  % correct trial RTs after excluding outliers
fxdI = cell(nSubj,1);
nBt = 1000;  % repetitions to get 95% confidence interval of the mean CE
wsdCE = nan(3,length(nSampling),nSubj);  % within-subject congruency effect distribution parameters. 3: 95% CI, skewness, kurtosis
for i = 1:nSubj
    dir = ['rawData/EPIC' num2str(subjID(i))];
    opts = detectImportOptions('session_numbering.xlsx');
    opts.Sheet = ['subj' num2str(subjID(i))];
    tmpMat = readmatrix('session_numbering.xlsx',opts);
    sessionID = tmpMat(:,taskIndx);
    sessionID = sessionID(1:nSession);  % throw away excessive sessions

    tmpC = [];  % accumulate correct congruent trials of all sessions
    tmpI = [];  % incongruent trials
    for j = 1:nSession
        fileID = [dir '/' taskStrng '/run' num2str(sessionID(j))];
        load(fileID)
        matn0cong = nan(length(allData),1);
        for k = 1:length(allData)
            if strcmp(allData(k,1),'CONGRUENT')
                matn0cong(k) = 1;
            elseif strcmp(allData(k,1),'INCONGRUENT')
                matn0cong(k) = 2;
            else
                error('n0cong string error: unable to identify current trial congruency')
            end
        end
        trialNum = (1:length(allData))';
        T = table(trialNum,matn0cong,double(cell2mat(allData(:,4))),double(cell2mat(allData(:,5))),'VariableNames',["trialNum","n0cong","acc","rt"]);
        tmpC = [tmpC; T.rt(and(T.n0cong==1,T.acc==1))];
        tmpI = [tmpI; T.rt(and(T.n0cong==2,T.acc==1))];
    end
    % Remove outlier trials
    tmpCgm = mean(tmpC,'omitnan');
    tmpIgm = mean(tmpI,'omitnan');
    tmpCsd = std(tmpC,'omitnan');
    tmpIsd = std(tmpI,'omitnan');
    lbound(1) = tmpCgm-nStd*tmpCsd;
    ubound(1) = tmpCgm+nStd*tmpCsd;
    lbound(2) = tmpIgm-nStd*tmpIsd;
    ubound(2) = tmpIgm+nStd*tmpIsd;
    fxdC{i,1} = tmpC(find(and(tmpC>=lbound(1),tmpC<=ubound(1)))).*1000;  % convert to msec
    fxdI{i,1} = tmpI(find(and(tmpI>=lbound(2),tmpI<=ubound(2)))).*1000;
    % Bootstrapping
    for i_n = 1:length(nSampling)
        simCE = [];
        for si = 1:nBt
            simCon = datasample(fxdC{i,1},nSampling(i_n));  % half
            simInc = datasample(fxdI{i,1},nSampling(i_n));
            simCE = [simCE, mean(simInc)-mean(simCon)];
        end
        wsdCE(1,i_n,i) = prctile(simCE,UpB)-prctile(simCE,LwB);
        wsdCE(2,i_n,i) = skewness(simCE);
        wsdCE(3,i_n,i) = kurtosis(simCE);
    end
end
ws_std_listCE = mean(wsdCE(1,:,:),3);
ws_skewnessCE = mean(mean(wsdCE(2,:,:),3));
ws_kurtosisCE = mean(mean(wsdCE(3,:,:),3));

%% Between-subject distribution parameters - Hedge et al.'s (2018) flanker task data
tmpSubjID1 = 1:50;
tmpSubjID2 = 1:62;
subjExcld1 = [6, 17, 37, 21, 33];  % based on Hedge et al.'s Readme.txt; additionally exclude 21 33
subjExcld2 = [25, 28, 34, 38, 56, 32];  % based on Readme.txt; additionally exclude 32
subjID1 = setdiff(tmpSubjID1, subjExcld1);
subjID2 = setdiff(tmpSubjID2, subjExcld2);
nSubj1 = length(subjID1);  % study 1
nSubj2 = length(subjID2);  % study 2
nSubjT = nSubj1+nSubj2;
nBlocks = 5;
nTrials = 144;  % per block
GrtMat = nan(nCond,nSubjT);  % grand mean RT
% Study 1 data
for i = 1:nSubj1
    dir = ['HedgeData/Study1-' taskStrng];
    cd(dir)
    % session 1
    fileID  = ['Study1_P' num2str(subjID1(i)) taskStrng num2str(1) '.csv'];
    rawData = readmatrix(fileID);
    T1 = table((1:nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
    T1(T1.n0cong==1,:) = [];  % exclude neutral condition
    % session 2
    fileID  = ['Study1_P' num2str(subjID1(i)) taskStrng num2str(2) '.csv'];
    rawData = readmatrix(fileID);
    T2 = table((1+nBlocks*nTrials:2*nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
    T2(T2.n0cong==1,:) = [];  % exclude neutral condition
    T = [T1;T2];  % combine session 1&2
    cd('../../')  % home directory
    % Exclude outliers
    acc_1_ids = find(T.acc==1);
    if length(find(T.acc==1))/length(table2array(T))*100 < 70
        disp(['Study1 ' taskStrng ' task subject#' num2str(subjID1(i)) ' has below 70% accuracy'])
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
end
% Study 2 data
for i = 1:nSubj2
    dir = ['HedgeData/Study2-' taskStrng];
    cd(dir)
    % session 1
    fileID  = ['Study2_P' num2str(subjID2(i)) taskStrng num2str(1) '.csv'];
    rawData = readmatrix(fileID);
    T1 = table((1:nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
    T1(T1.n0cong==1,:) = [];  % exclude neutral condition
    % session 2
    fileID  = ['Study2_P' num2str(subjID2(i)) taskStrng num2str(2) '.csv'];
    rawData = readmatrix(fileID);
    T2 = table((1+nBlocks*nTrials:2*nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
    T2(T2.n0cong==1,:) = [];  % exclude neutral condition
    T = [T1;T2];
    cd('../../')  % home directory
    % Exclude outliers
    acc_1_ids = find(T.acc==1);
    if length(find(T.acc==1))/length(table2array(T))*100 < 70
        disp(['Study2 ' taskStrng ' task subject#' num2str(subjID2(i)) ' has below 70% accuracy'])
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
end
%GrtMat = GrtMat.*1000;  % covert to msec; covert after simulating data
bs_meanCE = 1000*mean(GrtMat(2,:)-GrtMat(1,:));
bs_stdCE = 1000*std(GrtMat(2,:)-GrtMat(1,:)); % ~18ms
bs_skewCE = skewness(GrtMat(2,:)-GrtMat(1,:));
bs_kurtCE = kurtosis(GrtMat(2,:)-GrtMat(1,:));

%% Simulation 1. Fixed between-subject mean and standard deviation
% Preassignment
ws_exper_mean = nan(n_subs,numTest,length(ws_std_listCE));  % mean CE of a subject (across varying ws std)
bs_exper_std = nan(numTest,length(ws_std_listCE));  % (observed) CE between-subject std
for l = 1:numTest
    % Generate a distribution of random values with the between-subject estimates above
    bs_dist = pearsrnd(bs_meanCE,bs_stdCE,bs_skewCE,bs_kurtCE,bs_dist_n,1);

    % loop through different within-subject variances
    for ws = 1:length(ws_std_listCE)
        ws_std = ws_std_listCE(ws);

        for n = 1:n_subs
            disp(['Testing' num2str(l) ': Subject#' num2str(n) ' step' num2str(ws)])
            % simulate a session as a random draw from this distribution
            ws_mean = bs_dist(randi(bs_dist_n,1,1));
            ws_dist = pearsrnd(ws_mean,ws_std,ws_skewnessCE,ws_kurtosisCE,ws_dist_n,1);
            ws_exper_mean(n,l,ws) = ws_dist(randi(ws_dist_n,1,1));
        end
    end
end
bs_exper_std(:,:) = squeeze(std(ws_exper_mean,'omitnan'));
mean_bs_exper_std2 = mean(bs_exper_std,1);  % mean std of 100 testings
CI_bs_exper_std2 = (prctile(bs_exper_std,UpB)-prctile(bs_exper_std,LwB))/2;  % 95% CI

% Plot the results (Figure 6A)
cmap = hsv(length(bs_std_list));  % colormap
figure
a = errorbar(ws_std_listCE,mean_bs_exper_std2,CI_bs_exper_std2,'LineWidth',2.4); hold on
a.Marker = 'o';
%a.MarkerSize = 4;
a.Color = cmap(3,:);  % 3rd is 20 ms
a.CapSize = 15;
yline(bs_stdCE,'r--','LineWidth',2.4);
set(gca,'FontSize',16)
xlabel('Within-subject standard deviation (ms)','FontSize',20)
ylabel({'Apparent between-subject'; 'standard deviation (ms)'},'FontSize',20)
ylim([0 92])
%legend('','True between-subject standard deviation','FontSize',17)

%% Simulation 2. Varying between-subject standard deviation
% Preassignment
bs_std_cell = nan(length(bs_std_list),length(ws_std_listCE));
conIntvlh_cell = nan(length(bs_std_list),length(ws_std_listCE));
% loop through different true bs variances
for b = 1:length(bs_std_list)
    bs_std = bs_std_list(b);

    % loop through different ws variances
    ws_exper_sim = nan(n_subs,numTest,length(ws_std_listCE));
    for l = 1:numTest
        % generate a distribution of random values with the between-subject estimates above
        bs_dist = pearsrnd(bs_meanCE,bs_std,bs_skewCE,bs_kurtCE,bs_dist_n,1);
        % loop through different within-subject variances
        for ws = 1:length(ws_std_listCE)
            ws_std = ws_std_listCE(ws);
            for n = 1:n_subs  % sample a sub at random from the bs distribution with replacement
                disp(['BS step' num2str(b) ' - Testing' num2str(l) ': Subject#' num2str(n) ' WS step' num2str(ws)])
                % set that to the 'true' ws_mean for this sub
                ws_mean = bs_dist(randi(bs_dist_n,1,1));
                ws_dist = pearsrnd(ws_mean,ws_std,ws_skewnessCE,ws_kurtosisCE,ws_dist_n,1);
                ws_exper_sim(n,l,ws) = ws_dist(randi(ws_dist_n,1,1));
            end
        end
    end

    % now calculate the apparent bs variance from random draw
    sim_bs_std = squeeze(std(ws_exper_sim,'omitnan'));
    bs_std_cell(b,:) = mean(sim_bs_std);
    conIntvlh_cell(b,:) = (prctile(sim_bs_std,UpB)-prctile(sim_bs_std,LwB))/2;
end

% Plot results (Figure 6B)
figure
for b = 1:length(bs_std_list)
    a = errorbar(ws_std_listCE,bs_std_cell(b,:),conIntvlh_cell(b,:),'LineWidth',2.5); hold on  % random draw 원하면 bs_std_cell 대신 bs_exper_std_RD
    a.Marker = 'o';
    %a.MarkerSize = 4;
    a.Color = cmap(b,:);
    a.CapSize = 15;
end
set(gca,'FontSize',16)
xlabel('Within-subject standard deviation (ms)','FontSize',20)
ylabel({'Apparent between-subject'; 'standard deviation (ms)'},'FontSize',20);
ylim([0 92])
%lgd = legend('5 ms','10 ms','20 ms','30 ms','40 ms','50 ms','60 ms','Location','southeast','FontSize',12);
%title(lgd,{'True between-subject'; 'standard deviation'},'FontSize',12)

%% Simulation 3. Different numbers of simulated subjects
% Preassignment
bs_std_cell2_2 = nan(length(n_subs_list),length(ws_std_listCE));
conIntvlh_cell2_2 = nan(length(n_subs_list),length(ws_std_listCE));
% loop through different true bs variances
for s = 1:length(n_subs_list)
    n_subs = n_subs_list(s);

    % loop through different ws variances
    ws_exper_sim = nan(n_subs,numTest,length(ws_std_listCE));

    for l = 1:numTest
        bs_dist = pearsrnd(bs_meanCE,bs_stdCE,bs_skewCE,bs_kurtCE,bs_dist_n,1);
        for ws = 1:length(ws_std_listCE)
            ws_std = ws_std_listCE(ws);
            for n = 1:n_subs  % sample a sub at random from the bs distribution with replacement
                disp(['Sample size step' num2str(s) ' - Testing' num2str(l) ': Subject#' num2str(n) ' WS step' num2str(ws)])
                % set that to the 'true' ws_mean for this sub
                ws_mean = bs_dist(randi(bs_dist_n,1,1));
                ws_dist = pearsrnd(ws_mean,ws_std,ws_skewnessCE,ws_kurtosisCE,ws_dist_n,1);
                ws_exper_sim(n,l,ws) = ws_dist(randi(ws_dist_n,1,1));
            end
        end
    end

    % now calculate the apparent bs variance
    % BS calculated as the mean of 100 simulations
    sim_bs_std = squeeze(std(ws_exper_sim,'omitnan'));
    bs_std_cell2_2(s,:) = mean(sim_bs_std);
    conIntvlh_cell2_2(s,:) = (prctile(sim_bs_std,UpB)-prctile(sim_bs_std,LwB))/2;
end

% Plot the results (Figure 6C)
figure
cmap = turbo(length(n_subs_list));
for s = 1:length(n_subs_list)
    a = errorbar(ws_std_listCE,bs_std_cell2_2(s,:),conIntvlh_cell2_2(s,:),'LineWidth',2.4); hold on
    a.Marker = 'o';
    %a.MarkerSize = 4;
    a.Color = cmap(s,:);
    a.CapSize = 15;
end
set(gca,'FontSize',16)
xlabel('Within-subject standard deviation (ms)','FontSize',20)
ylabel({'Apparent between-subject'; 'standard deviation (ms)'},'FontSize',20)
ylim([0 92])
lgd = legend('50','100','200','300','400','500','1,000','Location','northwest','FontSize',16);
%title(lgd,'Sample size','FontSize',16)

%% Save data
save('sim_varResults_hlf.mat',...
    'mean_bs_exper_std2','CI_bs_exper_std2',...
    'bs_std_cell','conIntvlh_cell',...
    'bs_std_cell2_2','conIntvlh_cell2_2')