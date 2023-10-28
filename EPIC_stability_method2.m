%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
% flanker, prime-probe, and Stroop tasks
%
% To run this script:
% You need the EPIC dataset and "session_numbering.xlsx"
% (download both at https://osf.io/jk9nb)
%
% What this script does:
% Draws stability curves with method 2 (width of 95% confidence interval of
% the mean congruency effect [CE]) for RT, accuracy, and IES data
% Drops the INITIAL two blocks (for all measures: RT, Acc, & IES CE) and
% regresses out practice effects with a linear model (for RT CE only) 
% STEPS:
% 1. Break data up into 144 trial units (4blocks*2(in half)*18sessions)
% 2. Add a unit (randomly choose one with replacement)
% 3. Calculate the mean congruency effect (CE)
% 4. Repeat steps 2-3 for 144 steps
% 5. Repeat steps 2-4 5000 times and get the 95% CI of the 5000 scores
% 6. Plot the width of the CI
% 7. Find the knee point or abruptly changing point%
% n trials added on evey step(=144)
% n = 50 trials for Flanker (con25 + inc25)
% n = 48 trials for Prime-probe (con24 + inc24)
% n = 36 trials for Stroop task (con18 + inc18)
%
% What this script outputs:
% A figure showing before and after linear regression (practice effects)
% Supp. Fig. 7. CE RT
% Supp. Fig. 8-2 CE RT without regressing practice effect
% Supp. Fig. 9-2 CE accuracy
% Supp. Fig. 10-2 CE IES
%
% Created on 10/28/2022 by HJ Lee
% Last modifed on 06/30/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')
%% Parameter settings
tmpSubjID = 1:12;
subjExcld = [1,2,11,9];  % excluded participants
subjID = setdiff(tmpSubjID, subjExcld);
nSubj = length(subjID);

nTask = 4-1;  % 1: Flanker, 2: Prime-Probe, 3: Stroop
taskIndx=[1,3,4];  % Three of the 4 tasks (1=FL, 3=PP, 4=ST) session_numbering.xlsx
taskStrng = {'Flanker','PrimeProbe','STROOP'};
nSession = 18;  % People do it for 9 weeks and each task is used during two sessions each week
nBlocks = 4;
nBlocksT = nSession*nBlocks;  % 72
bD = 2;  % block division; half the block
nTrials = [100; 96; 108]./bD;  % half a block; this amount is added on each step
nSteps = nBlocksT*bD;  % number of steps
numTest = 5000;  % number of repetition
nCond = 2;  % congruent(1) vs. incongruent(2)
UpB = 97.5;
LwB = 2.5;

%% Load raw data and break them up by half the blocks
% Preallocation
% (1) RT - Store only the correctly responded RTs
conMSrt = cell(nSteps,nTask,nSubj);  % 144
incMSrt = cell(nSteps,nTask,nSubj);
% (2) Acc
conMSacc = cell(nSteps,nTask,nSubj);
incMSacc = cell(nSteps,nTask,nSubj);
for i = 1:nSubj
    for t = 1:nTask
        %% Load session indexing information
        dir = ['rawData/EPIC' num2str(subjID(i))];
        % session information
        opts = detectImportOptions('session_numbering.xlsx');
        opts.Sheet = ['subj' num2str(subjID(i))];
        tmpMat = readmatrix('session_numbering.xlsx',opts);
        sessionID = tmpMat(:,taskIndx(t));
        sessionID = sessionID(1:18);
        for j = 1:nSession
            %% Load data
            fileID = [dir '/' taskStrng{t} '/run' num2str(sessionID(j))];
            load(fileID)
            %% Current trial congruency
            matn0cong = nan(length(allData),1);
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
                T = table(matn0cong,double(cell2mat(allData(:,6))),double(cell2mat(allData(:,5))),'VariableNames',["n0cong","acc","rt"]);
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
                T = table(matn0cong,double(cell2mat(allData(:,4))),double(cell2mat(allData(:,5))),'VariableNames',["n0cong","acc","rt"]);
            end
            %% Break up data
            T1 = T(1:nTrials(t),:);
            T2 = T(nTrials(t)+1:2*nTrials(t),:);
            T3 = T(2*nTrials(t)+1:3*nTrials(t),:);
            T4 = T(3*nTrials(t)+1:4*nTrials(t),:);
            T5 = T(4*nTrials(t)+1:5*nTrials(t),:);
            T6 = T(5*nTrials(t)+1:6*nTrials(t),:);
            T7 = T(6*nTrials(t)+1:7*nTrials(t),:);
            T8 = T(7*nTrials(t)+1:8*nTrials(t),:);
            %% Data assignment
            % RT
            conMSrt{(j-1)*nBlocks*bD+1,t,i} = 1000*T1.rt(and(T1.n0cong==1,T1.acc==1));  % convert to msecs
            incMSrt{(j-1)*nBlocks*bD+1,t,i} = 1000*T1.rt(and(T1.n0cong==2,T1.acc==1));
            conMSrt{(j-1)*nBlocks*bD+2,t,i} = 1000*T2.rt(and(T2.n0cong==1,T2.acc==1));
            incMSrt{(j-1)*nBlocks*bD+2,t,i} = 1000*T2.rt(and(T2.n0cong==2,T2.acc==1));
            conMSrt{(j-1)*nBlocks*bD+3,t,i} = 1000*T3.rt(and(T3.n0cong==1,T3.acc==1));
            incMSrt{(j-1)*nBlocks*bD+3,t,i} = 1000*T3.rt(and(T3.n0cong==2,T3.acc==1));
            conMSrt{(j-1)*nBlocks*bD+4,t,i} = 1000*T4.rt(and(T4.n0cong==1,T4.acc==1));
            incMSrt{(j-1)*nBlocks*bD+4,t,i} = 1000*T4.rt(and(T4.n0cong==2,T4.acc==1));
            conMSrt{(j-1)*nBlocks*bD+5,t,i} = 1000*T5.rt(and(T5.n0cong==1,T5.acc==1));
            incMSrt{(j-1)*nBlocks*bD+5,t,i} = 1000*T5.rt(and(T5.n0cong==2,T5.acc==1));
            conMSrt{(j-1)*nBlocks*bD+6,t,i} = 1000*T6.rt(and(T6.n0cong==1,T6.acc==1));
            incMSrt{(j-1)*nBlocks*bD+6,t,i} = 1000*T6.rt(and(T6.n0cong==2,T6.acc==1));
            conMSrt{(j-1)*nBlocks*bD+7,t,i} = 1000*T7.rt(and(T7.n0cong==1,T7.acc==1));
            incMSrt{(j-1)*nBlocks*bD+7,t,i} = 1000*T7.rt(and(T7.n0cong==2,T7.acc==1));
            conMSrt{(j-1)*nBlocks*bD+8,t,i} = 1000*T8.rt(and(T8.n0cong==1,T8.acc==1));
            incMSrt{(j-1)*nBlocks*bD+8,t,i} = 1000*T8.rt(and(T8.n0cong==2,T8.acc==1));
            % Acc
            conMSacc{(j-1)*nBlocks*bD+1,t,i} = T1.acc(T1.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+1,t,i} = T1.acc(T1.n0cong==2);
            conMSacc{(j-1)*nBlocks*bD+2,t,i} = T2.acc(T2.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+2,t,i} = T2.acc(T2.n0cong==2);
            conMSacc{(j-1)*nBlocks*bD+3,t,i} = T3.acc(T3.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+3,t,i} = T3.acc(T3.n0cong==2);
            conMSacc{(j-1)*nBlocks*bD+4,t,i} = T4.acc(T4.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+4,t,i} = T4.acc(T4.n0cong==2);
            conMSacc{(j-1)*nBlocks*bD+5,t,i} = T5.acc(T5.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+5,t,i} = T5.acc(T5.n0cong==2);
            conMSacc{(j-1)*nBlocks*bD+6,t,i} = T6.acc(T6.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+6,t,i} = T6.acc(T6.n0cong==2);
            conMSacc{(j-1)*nBlocks*bD+7,t,i} = T7.acc(T7.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+7,t,i} = T7.acc(T7.n0cong==2);
            conMSacc{(j-1)*nBlocks*bD+8,t,i} = T8.acc(T8.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+8,t,i} = T8.acc(T8.n0cong==2);
        end
    end
end

%% Inverse Efficiency Scores (IES)
conMSies = nan(nTask,nSteps,nSubj);
incMSies = nan(nTask,nSteps,nSubj);
for i = 1:nSubj
    for s = 1:nSteps
        conMSies(1,s,i) = mean(conMSrt{s,1,i})/mean(conMSacc{s,1,i});  % FL
        conMSies(2,s,i) = mean(conMSrt{s,2,i})/mean(conMSacc{s,2,i});  % PP
        conMSies(3,s,i) = mean(conMSrt{s,3,i})/mean(conMSacc{s,3,i});  % ST
        incMSies(1,s,i) = mean(incMSrt{s,1,i})/mean(incMSacc{s,1,i});
        incMSies(2,s,i) = mean(incMSrt{s,2,i})/mean(incMSacc{s,2,i});
        incMSies(3,s,i) = mean(incMSrt{s,3,i})/mean(incMSacc{s,3,i});
    end
end
bCEmatIES = incMSies-conMSies;

%% Mean CE before linear regression - RT
bCEmatRT = nan(nTask,nSteps,nSubj);  % before (nSteps==144)
bCEmatACC = nan(nTask,nSteps,nSubj);  % this will not be regressed; CE acc does not have practice effects
for i = 1:nSubj
    for s = 1:nSteps
        bCEmatRT(1,s,i) = mean(incMSrt{s,1,i})-mean(conMSrt{s,1,i});
        bCEmatRT(2,s,i) = mean(incMSrt{s,2,i})-mean(conMSrt{s,2,i});
        bCEmatRT(3,s,i) = mean(incMSrt{s,3,i})-mean(conMSrt{s,3,i});
        bCEmatACC(1,s,i) = mean(conMSacc{s,1,i})-mean(incMSacc{s,1,i});
        bCEmatACC(2,s,i) = mean(conMSacc{s,2,i})-mean(incMSacc{s,2,i});
        bCEmatACC(3,s,i) = mean(conMSacc{s,3,i})-mean(incMSacc{s,3,i});
    end
end

%% Regress the practice effect - RT
aCEmatRT = nan(nTask,nSteps,nSubj);  % after
for t = 1:nTask
    tempM = squeeze(bCEmatRT(t,:,:))';
    [r,c] = size(tempM);
    cons = ones(c,1);
    linreg = [ones(c,1) linspace(0,1,c)'];
    tempMcell = num2cell(tempM(:,logical(cons))',1);
    linregcell = repmat({linreg(logical(cons),:)},[1 r]);
    tempbetas = cellfun(@mldivide,linregcell,tempMcell,'uniformoutput',0);
    tempbetas = cell2mat(tempbetas);
    tempbetas = tempbetas';
    tmpintvals = tempbetas*linreg';
    %tempT = tempT-tmpintvals;
    regM = tempM-tmpintvals+repmat(mean(tempM,2),[1,c]);
    aCEmatRT(t,:,:) = regM';
end

%% Plot the mean CE (after linear regression)
cmap = turbo(nSubj);
for t = 1:nTask
    figure
    subplot(1,2,1)
    for i = 1:nSubj
        plot(1:nSteps,bCEmatRT(t,:,i),'Color',cmap(i,:),'LineWidth',1.2); hold on
    end
    xlabel('Trial units')
    ylabel('RT CE (ms)')
    title('Before linear regression')
    subplot(1,2,2)
    for i = 1:nSubj
        plot(1:nSteps,aCEmatRT(t,:,i),'Color',cmap(i,:),'LineWidth',1.2); hold on
    end
    xlabel('Trial units')
    ylabel('RT CE (ms)')
    title('After')
    sgtitle([taskStrng{t} ' Task'])
end

%% Remove the initial 4 steps due to large variance in the beginning of experiment
rmvL = 1:4;  % nTrials*4 = the amount removing 200(FL), 192(PP), 216(ST)
conMSrt(rmvL,:,:) = [];
incMSrt(rmvL,:,:) = [];
conMSacc(rmvL,:,:) = [];
incMSacc(rmvL,:,:) = [];
aCEmatRT(:,rmvL,:) = [];  % rt; practice effect regressed
bCEmatRT(:,rmvL,:) = [];  % rt; practice effect NOT regressed
bCEmatACC(:,rmvL,:) = [];  % acc; practice effect NOT regressed
bCEmatIES(:,rmvL,:) = [];  % ies; practice effect NOT regressed
nSteps = nSteps-length(rmvL);  % update nSteps

%% Matrix definition and preassignment
conIntvl = nan(nTask,nSteps,nSubj);
conIntvl_test = nan(nTask,nSteps,nSubj);  % for testing
conNTi = nan(nTask,nSteps,nSubj);  % mean number of trials added; for testing and plotting
incNTi = nan(nTask,nSteps,nSubj);

%% RT or Accuracy? Run this script FOUR times from this line; change rtAcc and practice
rtAcc = 1;  % rt:1, acc:2, IES:3
practice = 1;  % regressed:1 (only when rtAcc=1 and you want to regress out practice effect), no:0

if rtAcc == 1
    if practice == 1
        trialUnit = aCEmatRT;  % CE after LINEAR REGRESSION
    elseif practice == 0
        trialUnit = bCEmatRT;  % if you don't want to regress out the practice effect
    end
elseif rtAcc == 2
    trialUnit = bCEmatACC;
elseif rtAcc == 3
    trialUnit = bCEmatIES;
end
%% Bootstrapping
for i = 1:nSubj
    CEtmpMat = nan(numTest,nSteps,nTask);
    CEtmpMat_test = nan(numTest,nSteps,nTask);  % for testing
    conNT = nan(numTest,nSteps,nTask);  % for testing and plotting
    incNT = nan(numTest,nSteps,nTask);  % for testing and plotting
    for k = 1:numTest
        % Make test set
        TestFL = [];
        TestPP = [];
        TestST = [];
        TestConFL = [];  % for testing
        TestIncFL = [];
        TestConPP = [];
        TestIncPP = [];
        TestConST = [];
        TestIncST = [];
        
        testorder = datasample((1:nSteps),nSteps);  % randomize(w/replacement) order of adding a trial unit on each step
        for s = 1:nSteps
            TestFL = [TestFL; trialUnit(1,testorder(s),i)];
            TestPP = [TestPP; trialUnit(2,testorder(s),i)];
            TestST = [TestST; trialUnit(3,testorder(s),i)];
            CEtmpMat(k,s,1) = mean(TestFL);
            CEtmpMat(k,s,2) = mean(TestPP);
            CEtmpMat(k,s,3) = mean(TestST);
            TestConFL = [TestConFL; conMSrt{testorder(s),1,i}];  % for testing
            TestIncFL = [TestIncFL; incMSrt{testorder(s),1,i}];
            TestConPP = [TestConPP; conMSrt{testorder(s),2,i}];
            TestIncPP = [TestIncPP; incMSrt{testorder(s),2,i}];
            TestConST = [TestConST; conMSrt{testorder(s),3,i}];
            TestIncST = [TestIncST; incMSrt{testorder(s),3,i}];
            CEtmpMat_test(k,s,1) = mean(TestIncFL)-mean(TestConFL);  % for testing
            CEtmpMat_test(k,s,2) = mean(TestIncPP)-mean(TestConPP);
            CEtmpMat_test(k,s,3) = mean(TestIncST)-mean(TestConST);
            conNT(k,s,1) = length(TestConFL);  % for testing and plotting
            incNT(k,s,1) = length(TestIncFL);
            conNT(k,s,2) = length(TestConPP);
            incNT(k,s,2) = length(TestIncPP);
            conNT(k,s,3) = length(TestConST);
            incNT(k,s,3) = length(TestIncST);
        end
    end
    % 95% Confidence Intervals of the simulated scores
    conIntvl(1,:,i) = prctile(CEtmpMat(:,:,1),UpB)-prctile(CEtmpMat(:,:,1),LwB);
    conIntvl(2,:,i) = prctile(CEtmpMat(:,:,2),UpB)-prctile(CEtmpMat(:,:,2),LwB);
    conIntvl(3,:,i) = prctile(CEtmpMat(:,:,3),UpB)-prctile(CEtmpMat(:,:,3),LwB);
    conIntvl_test(1,:,i) = prctile(CEtmpMat_test(:,:,1),UpB)-prctile(CEtmpMat_test(:,:,1),LwB);  % for testing
    conIntvl_test(2,:,i) = prctile(CEtmpMat_test(:,:,2),UpB)-prctile(CEtmpMat_test(:,:,2),LwB);
    conIntvl_test(3,:,i) = prctile(CEtmpMat_test(:,:,3),UpB)-prctile(CEtmpMat_test(:,:,3),LwB);
    conNTi(1,:,i) = mean(conNT(:,:,1));  % for testing and plotting
    incNTi(1,:,i) = mean(incNT(:,:,1));
    conNTi(2,:,i) = mean(conNT(:,:,2));
    incNTi(2,:,i) = mean(incNT(:,:,2));
    conNTi(3,:,i) = mean(conNT(:,:,3));
    incNTi(3,:,i) = mean(incNT(:,:,3));
end

% if and(rtAcc==1,practice==1)
%     conIntvl_1_1 = conIntvl;
%     save('RTstabilityCurveVariables2_initial_practice','conIntvl_1_1','conNTi','incNTi')
% elseif and(rtAcc==1,practice==0)
%     conIntvl_1_2 = conIntvl;
%     save('RTstabilityCurveVariables2_initial','conIntvl_1_2','conNTi','incNTi')
% elseif rtAcc==2
%     conIntvl_2 = conIntvl;
%     save('ACCstabilityCurveVariables2_initial','conIntvl_2','conNTi','incNTi')
% elseif rtAcc==3
%     conIntvl_3 = conIntvl;
%     save('IESstabilityCurveVariables2_initial','conIntvl_3','conNTi','incNTi')
% end

% if and(rtAcc==1,practice==1)
%     load RTstabilityCurveVariables2_initial_practice
%     conIntvl = conIntvl_1_1;
% elseif and(rtAcc==1,practice==0)
%     load RTstabilityCurveVariables2_initial
%     conIntvl = conIntvl_1_2;
% elseif rtAcc==2
%     load ACCstabilityCurveVariables2_initial
%     conIntvl = conIntvl_2;
% elseif rtAcc==3
%     load IESstabilityCurveVariables2_initial
%     conIntvl = conIntvl_3;
% end

%% Plot the results
% Stability curve
taskStrngfx = {'Flanker','Prime-Probe','Stroop'};
stp = 140;  % 144-4; nSteps
nihTB = 40;  % number of trials of NIH toolbox flanker task
Ebg = 96;  % number of trials of Eisenberg et al.'s (2019) Stroop task
for t = 1:nTask
    figure
    for i = 1:nSubj
        x = conNTi(t,1:stp,i)+incNTi(t,1:stp,i);  % # of trials of the test set sample
        h = plot(x,conIntvl(t,1:stp,i),'Color',cmap(i,:),'LineWidth',1.2); hold on
        set(get(h,'Parent'),'XScale','log')
    end
    xlabel('Number of trials','FontSize',14)
    xlim([0 8000])  % 18session*400trials = 7200
    if rtAcc == 1
        ylabel('Width of 95% CI (ms)','FontSize',14)
        if t == 1
            ylim([0 110])
        elseif t == 2
            ylim([0 150])
        elseif t == 3
            ylim([0 250])
        end
    elseif rtAcc == 2
        ylabel('Width of 95% CI','FontSize',14)
        ylim([0 0.35]);
    elseif rtAcc == 3
        ylabel('Width of 95% CI')
        if t == 1
            ylim([0 200])
        elseif t == 2
            ylim([0 300])
        elseif t == 3
            ylim([0 400])
        end
    end
    xticks([50 100 200 400 800 1600 3200 6400])
    title([taskStrngfx{t} ' Task'],'FontSize',15)
    if x(1) < nihTB
        xline(nihTB,'k','LineWidth',0.7,'LineStyle','--')
    else
        xline(x(1),'k','LineWidth',0.7,'LineStyle','--')
    end
    xline(Ebg,'k','LineWidth',0.7,'LineStyle','-.')
    grid on

    if t == 1
        legend('Participant 03','Participant 04','Participant 05','Participant 06',...
            'Participant 07','Participant 08','Participant 10','Participant 12',....
            'NIH toolbox','Eisenberg et al.','FontSize',12)
    end
end