%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC - Public Data 1
%
%% Instructions for running the script:
% You will need Robinson & Steyvers' (2023) flanker task data, available
% for download at https://osf.io/6hjwv/
%
%% Purpose of the script:
% This script generates stability curves using Method 2 (calculating the 
% width of the 95% confidence interval of the mean congruency effect [CE]) 
% for reaction time (RT) and accuracy data. It also plots the 
% between-subject standard deviation as a function of the number of trials.
%
% The analysis excludes the first four sessions.
%
% Participant exclusion criteria:
% (a) Accuracy below 70% in either experimental conditions (congruent or incongruent)
% (b) 0% accuracy in any session
% (c) Fewer than 2500 correctly responded trials
%
%% Outputs:
% Figure 4.
%
% Created on 12/17/2022 by HJ Lee
% Last modified on 02/05/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')

%% Load Robinson & Steyvers' (2023) flanker task data
load FlankerData_learn

%% Parameters
tmpnSubj = length(d_flanker_tab);  % number of total participants: 495
taskStrng = 'Flanker';
UpB = 97.5;  % upper percentile of 95% confidence interval
LwB = 2.5;  % lower percentile
numTest = 5000;  % number of repetitions; for consistency with the EPIC data
nCond = 2;  % congruent(1), incongruent(2)

%% Plot the distribution. X-axis: number of trials, y-axis: mean CE - BEFORE excluding participants
nTrialMat = nan(tmpnSubj,1);  % number of trials; this differ across participants
nTrialMatc = nan(tmpnSubj,1);  % number of correct trials
mCEmat = nan(tmpnSubj,1);  % individual mean CE
iRTdist = cell(nCond,tmpnSubj);  % trial-level RTs
for i = 1:tmpnSubj
    nTrialMat(i) = size(d_flanker_tab{i,1},1);
    cM = d_flanker_tab{i,1}.response_time(and(double(d_flanker_tab{i,1}.compatible)==1,...
        d_flanker_tab{i,1}.accuracy==1));  % correct congruent
    iM = d_flanker_tab{i,1}.response_time(and(double(d_flanker_tab{i,1}.compatible)==0,...
        d_flanker_tab{i,1}.accuracy==1));  % correct incongruent
    mCEmat(i) = mean(iM,'omitnan')-mean(cM,'omitnan');
    iRTdist{1,i} = cM;
    iRTdist{2,i} = iM;
    nTrialMatc(i) = length(cM)+length(iM);
end
%figure
%scatter(nTrialMatc,mCEmat)
%xlabel('Number of correctly responded trials')
%ylabel('Mean congruency effect')

%% 1. Preprocess data
% Matrix definition and preassignment
nSess = nan(tmpnSubj,1);  % number of sessions each participant had (this differs across participants)
for i = 1:tmpnSubj
    nSess(i) = length(find(d_flanker_tab{i,1}.trial_num==1));
end
[Max,Imax] = max(nSess);
nStep = Max;  % to store them in a unified size of cell; this is the # of steps because one session data (~64 trials) will be added to the test set sample on each step
conMSrt = cell(nStep,tmpnSubj);
incMSrt = cell(nStep,tmpnSubj);
conMSacc = cell(nStep,tmpnSubj);
incMSacc = cell(nStep,tmpnSubj);

% Divide data across sessions and store them in matrices
ageBin = nan(tmpnSubj,1);  % Each participant's age
for i = 1:tmpnSubj
    tmpT = table(d_flanker_tab{i,1}.agebin,d_flanker_tab{i,1}.trial_num,...
        double(d_flanker_tab{i,1}.compatible),d_flanker_tab{i,1}.accuracy,...
        d_flanker_tab{i,1}.response_time,'VariableNames',["ageBin","trial_num","n0cong","acc","rt"]);
    indx = find(tmpT.trial_num==1);
    indx = [indx; length(table2array(tmpT))+1];
    ageBin(i) = unique(tmpT.ageBin);  % if more than one value, this should shoot error

    for b = 1:length(indx)-1
        tmpTd = tmpT(indx(b):indx(b+1)-1,:);
        conMSrt{b,i} = tmpTd.rt(and(tmpTd.n0cong==1,tmpTd.acc==1));  % congruent RT
        incMSrt{b,i} = tmpTd.rt(and(tmpTd.n0cong==0,tmpTd.acc==1));  % incongruent
        conMSacc{b,i} = tmpTd.acc(tmpTd.n0cong==1);  % accuracy
        incMSacc{b,i} = tmpTd.acc(tmpTd.n0cong==0);
    end
end
stepTL = nSess;

% Check how many trials on average are tested in each session
mean(nTrialMat./nSess)  % ~49

% Update the following variables after excluding participants: stepTL, conMSrt, incMSrt, conMSacc, incMSacc, ageBin
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
ageBin(subjExcld,:) = [];
nSubj = length(subjID);  % 448

%% Exclusion criteria (c): Remove participants with less than or equal to 2500 CORRECT trials for RT data
lrg = 2500;
nTrialMatexcld = nan(nSubj,1);
for i = 1:nSubj
    nTrialMatexcld(i) = length(cell2mat(conMSrt(1:stepTL(i),i)))+length(cell2mat(incMSrt(1:stepTL(i),i)));
end
sGrp = find(nTrialMatexcld>lrg);
l_sGrp = length(sGrp);  % 185

%% Plot the learning curves to see if any practice effects exist in Robinson's data
% ceLC = nan(nStep,l_sGrp);
% for i = 1:l_sGrp
%     for j = 1:stepTL(sGrp(i))
%         ceLC(j,i) = mean(cell2mat(incMSrt(j,sGrp(i))),'omitnan')-mean(cell2mat(conMSrt(j,sGrp(i))),'omitnan');
%     end
% end
% cmapL = turbo(l_sGrp);
% figure
% for i = 1:l_sGrp
%     plot(ceLC(1:stepTL(sGrp(i)),i),'Color',cmapL(i,:),'LineWidth',1); hold on
%     xlabel('Session')
%     ylabel('RT CE (ms)')
% end

%% Remove the initial 4 sessions
rmvL = 1:4;
nStep = nStep-length(rmvL);
stepTL = stepTL-length(rmvL);
conMSrt(rmvL,:) = [];
incMSrt(rmvL,:) = [];
conMSacc(rmvL,:) = [];
incMSacc(rmvL,:) = [];

%% 2. Method 2: Stability Curve with the Width of 95% Confidence Interval
% Matrix definition and preassignment
CEmat = nan(nStep,nSubj);  % estimated CE from bootstrapping
realCEmat = nan(nStep,nSubj);  % actual CE
conIntvl = nan(nStep,nSubj);  % 95% confidence interval of the mean CE
conNTi = nan(nStep,nSubj);  % number of trials added on each step, congruent
incNTi = nan(nStep,nSubj);  % incongruent

%% RT or Accuracy?
rtAcc = 1;  % rt:1, acc:2

%% Bootstrapping
for i = 1:nSubj
    if rtAcc == 1
        % RT
        conMS = conMSrt(1:stepTL(i),i);
        incMS = incMSrt(1:stepTL(i),i);
    elseif rtAcc == 2
        % Acc
        conMS = conMSacc(1:stepTL(i),i);
        incMS = incMSacc(1:stepTL(i),i);
    end
    conM = nan(numTest,stepTL(i));
    incM = nan(numTest,stepTL(i));
    conNT = nan(numTest,stepTL(i));  % number of accumulated trials on each step; will use this for x-axis values when plotting
    incNT = nan(numTest,stepTL(i));
    for k = 1:numTest
        %% Make test set: Add a block on every step
        TestCon = [];
        TestInc = [];

        testorder = datasample((1:stepTL(i)),stepTL(i));
        for s = 1:stepTL(i)
            %% Add block on every step
            TestCon = [TestCon;cell2mat(conMS(testorder(s),1))];
            TestInc = [TestInc;cell2mat(incMS(testorder(s),1))];
            conM(k,s) = mean(TestCon);
            incM(k,s) = mean(TestInc);

            %% Track the number of trials added
            conNT(k,s) = length(TestCon);
            incNT(k,s) = length(TestInc);
        end
    end
    % Vectors to store actual CE - Put this outside the numTest-loop
    tCon = [];
    tInc = [];
    for s = 1:stepTL(i)
        %% Actual CE
        tCon = [tCon;cell2mat(conMS(s,1))];
        tInc = [tInc;cell2mat(incMS(s,1))];
        if rtAcc == 1
            realCEmat(s,i) = mean(tInc)-mean(tCon);
        elseif rtAcc == 2
            realCEmat(s,i) = mean(tCon)-mean(tInc);
        end
    end
    %% Magnitude of the CE
    if rtAcc == 1
        % RT
        % Mean
        CEmat(1:stepTL(i),i) = mean(incM-conM);
        % 95% CI
        conIntvl(1:stepTL(i),i) = prctile(incM-conM,UpB)-prctile(incM-conM,LwB);
    elseif rtAcc == 2
        % Acc
        % Mean
        CEmat(1:stepTL(i),i) = mean(conM-incM);
        % 95% CI
        conIntvl(1:stepTL(i),i) = prctile(conM-incM,UpB)-prctile(conM-incM,LwB);
    end

    %% Mean number of trials added one ach step
    conNTi(1:stepTL(i),i) = mean(conNT);
    incNTi(1:stepTL(i),i) = mean(incNT);
end
% Save data
if rtAcc == 1
    save('RobinsonRTstabilityCurve2Variables','CEmat','realCEmat','conIntvl','conNTi','incNTi')
elseif rtAcc == 2
    save('RobinsonACCstabilityCurve2Variables','CEmat','realCEmat','conIntvl','conNTi','incNTi')
end

%% Load data for plotting
% if rtAcc == 1
%     load RobinsonRTstabilityCurve2Variables
% elseif rtAcc == 2
%     load RobinsonACCstabilityCurve2Variables
% end

%% Between-Subject Variance
% Use only the selected group (Participants with more than 2500)
sgX = conNTi(:,sGrp)+incNTi(:,sGrp);
[mX,iX] = max(nTrialMatexcld(sGrp));  % iX=55;
sgXx = sgX(:,iX);  % index of the longest (and the longest step length is 160)
sgRealCEmat = realCEmat(:,sGrp)';
bsSD = std(sgRealCEmat,'omitnan');  % between-subject standard deviation

%% Plot within-subject and between-subject variability (Figure 4)
tmpSG = nan(l_sGrp,1);  % to get the number 50??
for i = 1:l_sGrp
    tmpSG(i) = length(conIntvl(1:stepTL(sGrp(i)),sGrp(i)));
end
idx = min(tmpSG);  % 50: the max number without nan. 50steps~2500 trials, 수정 전에 이렇고 수정 후엔 39, initial two blocks drop한 다음엔 35
WI = nan(l_sGrp,idx);  % to get the group median of within-subject variability
figure
% Within-subject variability
subplot(2,1,1)  % width of 95% CI
for i = 1:l_sGrp
    id = sGrp(i);
    x = conNTi(1:stepTL(id),id)+incNTi(1:stepTL(id),id);
    h = plot(x,conIntvl(1:stepTL(id),id),'Color',[0.5 0.5 0.5],'LineWidth',1); hold on
    WI(i,:) = conIntvl(1:idx,id)';
end
h2 = plot(x(1:idx),median(WI),'Color','r','LineWidth',2);  % median
set(gca,'FontSize',14)
set(get(h2,'Parent'),'XScale','log')
xlabel('Number of trials','FontSize',15)
%xlim([0 3500])  % As it is
xlim([0 x(idx)])  % Cut off where there is fewer participants; less than 20 have more than 3500 trials
xticks([50 100 200 400 800 1600])
if rtAcc == 1
    ylabel({'Width of 95%'; 'confidence interval (ms)'},'FontSize',15)
elseif rtAcc == 2
    ylabel({'Width of 95%'; 'confidence interval'},'FontSize',15)
end
grid on
if rtAcc == 1
    title('A) Within-Subject Variability','FontSize',19)
elseif rtAcc == 2
    title('C) Within-Subject Variability','FontSize',19)
end

% Between-subject variability
subplot(2,1,2)
p = plot(sgXx,bsSD,'Color','b','LineWidth',2);
set(gca,'FontSize',14)
set(get(p,'Parent'),'XScale','log')
xlabel('Number of trials','FontSize',15)
%xlim([0 3500])  % As it is
xlim([0 x(idx)])  % Cut off where there is fewer participants
xticks([50 100 200 400 800 1600])
if rtAcc == 1
    ylabel({'Between-subject'; 'standard deviation (ms)'},'FontSize',15)
    ylim([20 70])  % [20 80]
elseif rtAcc == 2
    ylabel({'Between-subject'; 'standard deviation'},'FontSize',15)
    %ylim([0 0.055])
end
grid on
if rtAcc == 1
    title('B) Between-Subject Variability','FontSize',19)
elseif rtAcc == 2
    title('D) Between-Subject Variability','FontSize',19)
end
% if rtAcc == 1
%     sgtitle('Reaction Time','FontSize',25)
% elseif rtAcc == 2
%     sgtitle('Accuracy','FontSize',25)
% end