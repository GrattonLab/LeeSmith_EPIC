%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC - Public Data 1
%
%% Instructions for running the script:
% You will need Robinson & Steyvers' (2023) flanker task data, available
% for download at https://osf.io/6hjwv/
%
%% Purpose of the script:
% This script simulates the selected participants from Robinson and 
% Steyvers (2023) (185 participants with >= 2500 trials) using the Pearson 
% system to replicate and extend Figure 4. The goal is to examine the
% impact of repeated measures on within- and between-subject variability by
% comparing small and large trial samples.
%
%% Outputs:
% Figure 5. Comparison between model 1 (large within-subject variability) 
% and model 2 (small within-subject variability)
%
% Created on 01/12/2023 by HJ Lee
% Last modified on 02/05/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')

%% Load Robinson and Steyvers' (2023) flanker task data
load FlankerData_learn

%% Parameters
tmpnSubj = length(d_flanker_tab);  % number of total participants: 495
nTrialMat = nan(tmpnSubj,1);
for i = 1:tmpnSubj
    nTrialMat(i) = size(d_flanker_tab{i,1},1);
end
UpB = 97.5;  % upper percentile of 95% confidence interval
LwB = 2.5;  % lower percentile
numTest = 1000;
nCond = 2;  % congruent(1), incongruent(2)

%% 1. Preprocess data
% Matrix definition and preassignment
nSess = nan(tmpnSubj,1);  % number of sessions each participant had (this differs across participants)
for i = 1:tmpnSubj
    nSess(i) = length(find(d_flanker_tab{i,1}.trial_num==1));
end
[Max,IMax] = max(nSess);
nStep = Max;
conMSrt = cell(nStep,tmpnSubj);
incMSrt = cell(nStep,tmpnSubj);
conMSacc = cell(nStep,tmpnSubj);
incMSacc = cell(nStep,tmpnSubj);

% Divide data across sessions and store them in matrices
for i = 1:tmpnSubj
    tmpT = table(d_flanker_tab{i,1}.trial_num,...
        double(d_flanker_tab{i,1}.compatible),d_flanker_tab{i,1}.accuracy,...
        d_flanker_tab{i,1}.response_time,'VariableNames',["trial_num","n0cong","acc","rt"]);
    indx = find(tmpT.trial_num==1);
    indx = [indx; length(table2array(tmpT))+1];

    for b = 1:length(indx)-1
        tmpTd = tmpT(indx(b):indx(b+1)-1,:);
        conMSrt{b,i} = tmpTd.rt(and(tmpTd.n0cong==1,tmpTd.acc==1));
        incMSrt{b,i} = tmpTd.rt(and(tmpTd.n0cong==0,tmpTd.acc==1));
        conMSacc{b,i} = tmpTd.acc(tmpTd.n0cong==1);
        incMSacc{b,i} = tmpTd.acc(tmpTd.n0cong==0);
    end
end
stepTL = nSess;

%% Exclusion criteria (a): Remove participants with below 70% accuracy
cAccmat = nan(tmpnSubj,1);
iAccmat = nan(tmpnSubj,1);
for i = 1:tmpnSubj
    cM = d_flanker_tab{i,1}.accuracy(double(d_flanker_tab{i,1}.compatible)==1);
    iM = d_flanker_tab{i,1}.accuracy(double(d_flanker_tab{i,1}.compatible)==0);
    cAccmat(i) = mean(cM,'omitnan');
    iAccmat(i) = mean(iM,'omitnan');
end
thshlda = .7;
subjIDa = find(and(cAccmat>thshlda,iAccmat>thshlda));  % length: 487 (8 excluded)

%% Exclusion criteria (b): Remove participants with 0% accuracy session
% RT matrix will be empty if all responses are incorrect
excldID = [];
for i = 1:tmpnSubj
    for j = 1:nSess(i)
        if isempty(conMSrt{j,i})
            excldID = [excldID; i];
        end
        if isempty(incMSrt{j,i})
            excldID = [excldID; i];
        end
    end
end
subjIDb = unique(excldID);
subjID = setdiff(subjIDa,subjIDb);
subjExcld = setdiff(1:tmpnSubj,subjID);
stepTL(subjExcld,:) = [];
conMSrt(:,subjExcld) = [];
incMSrt(:,subjExcld) = [];
conMSacc(:,subjExcld) = [];
incMSacc(:,subjExcld) = [];
nSubj = length(subjID);  % 448

%% Concatenate data
tCrt = cell(nSubj,1);
tIrt = cell(nSubj,1);
nTrialMatExl = nan(nSubj,1);
for i = 1:nSubj
    tCrt{i} = cell2mat(conMSrt(1:stepTL(i),i));
    tIrt{i} = cell2mat(incMSrt(1:stepTL(i),i));
    %tCacc{i} = cell2mat(conMSacc(1:stepTL(i),i));
    %tIacc{i} = cell2mat(incMSacc(1:stepTL(i),i));
    nTrialMatExl(i) = length(tCrt{i})+length(tIrt{i});
end

%% Exclusion criteria (c): Remove participants with less than or equal to 2500 CORRECT trials for RT data
lrg = 2500;
sGrpID = find(nTrialMatExl>lrg);
l_sGrp = length(sGrpID);  % 185
RTdistC = tCrt(sGrpID);
RTdistI = tIrt(sGrpID);

% Within-subject and between-subject variability
tCEmat = nan(l_sGrp,1);
tmeanmat = nan(l_sGrp,nCond);  % nCond: congruent, incongruent mean separately
twsSE = nan(l_sGrp,1);  % within-subject variability - 95% CI
twsSD = nan(l_sGrp,nCond);
for i = 1:l_sGrp
    tCEmat(i,1) = mean(RTdistI{i,1})-mean(RTdistC{i,1});
    tmeanmat(i,1) = mean(RTdistC{i,1});
    tmeanmat(i,2) = mean(RTdistI{i,1});
    twsSD(i,1) = std(RTdistC{i,1});
    twsSD(i,2) = std(RTdistI{i,1});
    simCEtmp = [];
    for k = 1:numTest
        simContmp = datasample(RTdistC{i,1},length(RTdistC{i,1}));
        simInctmp = datasample(RTdistI{i,1},length(RTdistI{i,1}));
        simCEtmp = [simCEtmp, mean(simInctmp)-mean(simContmp)];
    end
    twsSE(i,1) = (prctile(simCEtmp,UpB)-prctile(simCEtmp,LwB))/2;
end
tmean = [];
tbsSD = std(tCEmat);  % between-subject variability: bs SD

%% 2. Simulate data
% Generate data using the Pearson system (https://www.mathworks.com/help/stats/pearsrnd.html)
% Parameter settings
nTu = [40 80 160 320 640 1000 2000 4000 8000 16000]./nCond;  % number of trials sampled from distribution
nSteps = length(nTu);

% Preassignment
simCEmat = nan(l_sGrp,numTest,nSteps);  % simulated CE
conIntvl = nan(l_sGrp,nSteps);  % 95% confidence interval of the mean CE
sCcell = cell(l_sGrp,1);
sIcell = cell(l_sGrp,1);
% Simulate
for i = 1:l_sGrp
    conM = nan(numTest,nSteps);
    incM = nan(numTest,nSteps);
    % Generate data - Distribution parameters
    sC= {mean(RTdistC{i}),std(RTdistC{i}),skewness(RTdistC{i}),kurtosis(RTdistC{i})}; sCcell{i,1} = sC;
    sI= {mean(RTdistI{i}),std(RTdistI{i}),skewness(RTdistI{i}),kurtosis(RTdistI{i})}; sIcell{i,1} = sI;
    for k = 1:numTest
        for s = 1:nSteps
            disp(['Subject#' num2str(i) ': Simulation' num2str(k) ' step' num2str(s)])
            % Draw samples from distribution simulated with Pearson system
            rC = pearsrnd(sC{:},nTu(s),1);
            rI = pearsrnd(sI{:},nTu(s),1);
            conM(k,s) = mean(rC);
            incM(k,s) = mean(rI);
        end
        % Simulated mean
        simCEmat(i,k,:) = incM(k,:)-conM(k,:);
    end
    % 95% confidence interval
    conIntvl(i,:) = prctile(incM-conM,UpB)-prctile(incM-conM,LwB);
end
save('simulation_RobinsonRT_final','simCEmat','conIntvl')
%load simulation_RobinsonRT_final

%% Between-subject standard deviation
simbsSD = squeeze(std(simCEmat,'omitnan'));
UBsd = prctile(simbsSD,UpB);
LBsd = prctile(simbsSD,LwB);
conIntvlsd = UBsd-LBsd;
conIntvlsdh = conIntvlsd/2;

bsSDtmp = squeeze(std(simCEmat));
bsSDmean = mean(bsSDtmp);

%% Plot data (Figure 5)
nCase = 2;  % large / small within-subject variability
cRg = cell(2,1);
cRg{1} = 1:5;
cRg{2} = 6:10;
plotnSubj = l_sGrp;
for mycase = 1:nCase
    figure
    % Within-subject variability
    subplot(2,1,1)  % width of 95% CI
    x = nCond*nTu(cRg{mycase});
    for i = 1:plotnSubj
        h = plot(x,conIntvl(i,cRg{mycase}),'LineWidth',1,'Color',[0.5 0.5 0.5]); hold on
    end
    h2 = plot(x,median(conIntvl(:,cRg{mycase})),'Color','r','LineWidth',2);
    set(gca,'FontSize',14)  % default is 10
    set(get(h2,'Parent'),'XScale','log')
    xlabel('Number of trials','FontSize',15)
    xticks(x)
    xlim([x(1) x(end)])
    ylim([0 500])
    ylabel({'Width of 95%'; 'confidence interval (ms)'},'FontSize',15)
    grid on
    if mycase == 1
        title('A) Within-Subject Variability','FontSize',19)
    elseif mycase == 2
        title('C) Within-Subject Variability','FontSize',19)
        legend('1 simulated subject','Group median','FontSize',14)
    end

    % Between-subject variability
    subplot(2,1,2)
    xconf = [x x(end:-1:1)];
    CIl = bsSDmean(cRg{mycase})-conIntvlsdh(cRg{mycase});  % replace bsSDmean with bsSD if you want 1 random draw
    CIu = bsSDmean(cRg{mycase})+conIntvlsdh(cRg{mycase});
    u = CIu;
    yconf = [CIl u(end:-1:1)];
    p = fill(xconf,yconf,'b');
    p.FaceColor = [0.8 0.8 1];
    p.EdgeColor = 'none';
    hold on
    p2 = plot(x,bsSDmean(cRg{mycase}),'Color','b','LineWidth',2);
    set(gca,'FontSize',14)
    set(get(p2,'Parent'),'XScale','log')
    xlabel('Number of trials','FontSize',15)
    xticks(x)
    xlim([x(1) x(end)])
    ylim([20 70])
    ylabel({'Between-subject'; 'standard deviation (ms)'},'FontSize',15)
    grid on
    if mycase == 1
        title('B) Between-Subject Variability','FontSize',19)
    elseif mycase == 2
        title('D) Between-Subject Variability','FontSize',19)
        legend('95% confidence interval','FontSize',14)
    end
%     if mycase == 1
%         sgtitle('Model 1: Small Trial Sampling','FontSize',19)
%     elseif mycase == 2
%         sgtitle('Model 2: Large Trial Sampling','FontSize',19)
%     end
end