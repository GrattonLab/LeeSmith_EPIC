%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
% flanker, prime-probe, and Stroop tasks
%
% To run this script:
% You need the EPIC dataset and "session_numbering.xlsx"
% (download both at https://osf.io/jk9nb)
%
% What this script does:
% Draws stability curves with method 1 (absolute difference between the
% reference set and the test set sample) for RT, accuracy, and IES data
% Drops the INITIAL two blocks (for all measures: RT, Acc, & IES CE) and
% regresses out practice effects with a linear model (for RT CE only)
% STEPS:
% 1. Breaks data up into 144 trial units (4blocks*2(in half)*18sessions)
% 2. Assign half the units to a reference set(=72units) and the rest to a
% test set
% 3. Add a trial unit from the test set by randomly choosing one with
% replacement
% 4. Calculate the absolute difference between the reference set CE and the
% test set CE
% 5. Repeat 3-4 for 72 steps (The test set will be built up step by step
% with an addition of a unit)
% 6. Repeat 2-5 for 5,000 times and get the mean and 95% CI of the 5,000
% scores of the absolute difference
% 7. Plot the data
%
% What this script outputs:
% A figure showing before and after linear regression (practice effects)
% Fig. 3 CE RT
% Supp. Fig. 8-1 CE RT without regressing practice effect
% Supp. Fig. 9-1 CE accuracy
% Supp. Fig. 10-1 CE IES
%
% Created on 10/19/2022 by HJ Lee
% Last modified on 06/30/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')
%% Parameter settings
% Subject information
tmpSubjID = 1:12;
subjExcld = [1,2,11,9];  % excluded participants
subjID = setdiff(tmpSubjID, subjExcld);
nSubj = length(subjID);

% Task information
nTask = 4-1;  % 1: Flanker, 2: Prime-Probe, 3: Stroop
taskIndx=[1,3,4];  % Three of the 4 tasks (1=FL, 3=PP, 4=ST) session_numbering.xlsx
taskStrng = {'Flanker','PrimeProbe','STROOP'};
nSession = 18;  % People do it for 9 weeks and each task is used during two sessions each week
nBlocks = 4;
nBlocksT = nSession*nBlocks;  % 72
bD = 2;  % break block into half
nTrials = [100; 96; 108]./bD;  % half a block; this amount is added on each step
nSteps = (nBlocksT*bD)/2;  % number of steps; add half the block to the test set on every step
numTest = 5000;  % number of repetition
nCond = 2;  % congruent(1) vs. incongruent(2)
UpB = 97.5;  % percentile; upper boundary of the confidence interval
LwB = 2.5;  % lower boundary

%% Load raw data and break them up by half the blocks
% Preallocation
% (1) RT - Store only the correctly responded RTs
conMSrt = cell(nSteps*2,nTask,nSubj);  % 72*2
incMSrt = cell(nSteps*2,nTask,nSubj);
% (2) Acc
conMSacc = cell(nSteps*2,nTask,nSubj);
incMSacc = cell(nSteps*2,nTask,nSubj);
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
conMSies = nan(nTask,nSteps*2,nSubj);
incMSies = nan(nTask,nSteps*2,nSubj);
for i = 1:nSubj
    for s = 1:nSteps*2
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
bCEmatRT = nan(nTask,nSteps*2,nSubj);  % before (nSteps==72)
bCEmatACC = nan(nTask,nSteps*2,nSubj);  % this will not be regressed; CE acc does not have practice effects
bCmatRTmean = nan(nTask,nSteps*2,nSubj);
bImatRTmean = nan(nTask,nSteps*2,nSubj);
bCmatRTstd = nan(nTask,nSteps*2,nSubj);
bImatRTstd = nan(nTask,nSteps*2,nSubj);
for i = 1:nSubj
    for s = 1:nSteps*2
        bCEmatRT(1,s,i) = mean(incMSrt{s,1,i})-mean(conMSrt{s,1,i});
        bCEmatRT(2,s,i) = mean(incMSrt{s,2,i})-mean(conMSrt{s,2,i});
        bCEmatRT(3,s,i) = mean(incMSrt{s,3,i})-mean(conMSrt{s,3,i});
        bCEmatACC(1,s,i) = mean(conMSacc{s,1,i})-mean(incMSacc{s,1,i});
        bCEmatACC(2,s,i) = mean(conMSacc{s,2,i})-mean(incMSacc{s,2,i});
        bCEmatACC(3,s,i) = mean(conMSacc{s,3,i})-mean(incMSacc{s,3,i});

        % To store parameter values
        bCmatRTmean(1,s,i) = mean(conMSrt{s,1,i});
        bCmatRTmean(2,s,i) = mean(conMSrt{s,2,i});
        bCmatRTmean(3,s,i) = mean(conMSrt{s,3,i});
        bCmatRTstd(1,s,i) = std(conMSrt{s,1,i});
        bCmatRTstd(2,s,i) = std(conMSrt{s,2,i});
        bCmatRTstd(3,s,i) = std(conMSrt{s,3,i});
        bImatRTmean(1,s,i) = mean(incMSrt{s,1,i});
        bImatRTmean(2,s,i) = mean(incMSrt{s,2,i});
        bImatRTmean(3,s,i) = mean(incMSrt{s,3,i});
        bImatRTstd(1,s,i) = std(incMSrt{s,1,i});
        bImatRTstd(2,s,i) = std(incMSrt{s,2,i});
        bImatRTstd(3,s,i) = std(incMSrt{s,3,i});
    end
end

%% Regress the practice effect - RT
aCEmatRT = nan(nTask,nSteps*2,nSubj);  % after
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
        plot(1:nSteps*2,bCEmatRT(t,:,i),'Color',cmap(i,:),'LineWidth',1.2); hold on
    end
    xlabel('Trial units')
    ylabel('RT CE (ms)')
    title('Before linear regression')
    subplot(1,2,2)
    for i = 1:nSubj
        plot(1:nSteps*2,aCEmatRT(t,:,i),'Color',cmap(i,:),'LineWidth',1.2); hold on
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
nSteps = nSteps-(length(rmvL))/2;  % update nSteps

%% Matrix definition and preassignment
aProgMean = nan(nTask,nSteps,nSubj);
aProgMean_test = nan(nTask,nSteps,nSubj);  % for testing
aProgSD = nan(nTask,nSteps,nSubj);
aProgUB = nan(nTask,nSteps,nSubj);
aProgLB = nan(nTask,nSteps,nSubj);
conNTi = nan(nTask,nSteps,nSubj);  % number of trials - congruent condition; for testing and plotting
incNTi = nan(nTask,nSteps,nSubj);  % incongruent

%% RT or Accuracy? Run this script FOUR times from this line; change rtAcc and practice
rtAcc = 1;  % rt:1, acc:2, IES:3
practice = 1;  % regressed:1 (only when rtAcc=1 and you want to regress out practice effect), no:0

if rtAcc == 1
    if practice == 1
        trialUnit = aCEmatRT;  % LINEAR REGRESSION
    elseif practice == 0
        trialUnit = bCEmatRT;  % if you don't want to regress out practice effect
    end
elseif rtAcc == 2
    trialUnit = bCEmatACC;
elseif rtAcc == 3
    trialUnit = bCEmatIES;
end
%% Bootstrapping
AsignBase=[zeros(1,nSteps),ones(1,nSteps)];  % Os for the reference and 1s for the test set
for i = 1:nSubj
    absdiffvec = nan(numTest,nSteps,nTask);
    absdiffvec_test = nan(numTest,nSteps,nTask);  % for testing
    conNT = nan(numTest,nSteps,nTask);  % for testing and plotting
    incNT = nan(numTest,nSteps,nTask);  % for testing and plotting
    for k = 1:numTest
        %% Make reference set & test set: Block assignment in random orders
        %Randomize the block that is added for each step (Test Set)%
        % Flanker
        blockreasign = randperm(nSteps*2); % scramble the block order
        Asigner = AsignBase(blockreasign);
        RefFL = trialUnit(1,Asigner==0,i);  % reference set
        TestFL_T = trialUnit(1,Asigner==1,i);  % test set
        conRefFL = conMSrt(Asigner==0,1,i);  % for testing
        incRefFL = incMSrt(Asigner==0,1,i);
        conTestFL_T = conMSrt(Asigner==1,1,i);
        incTestFL_T = incMSrt(Asigner==1,1,i);

        % Prime-probe
        blockreasign = randperm(nSteps*2);
        Asigner = AsignBase(blockreasign);
        RefPP = trialUnit(2,Asigner==0,i);
        TestPP_T = trialUnit(2,Asigner==1,i);
        conRefPP = conMSrt(Asigner==0,2,i);  % for testing
        incRefPP = incMSrt(Asigner==0,2,i);
        conTestPP_T = conMSrt(Asigner==1,2,i);
        incTestPP_T = incMSrt(Asigner==1,2,i);

        % Stroop
        blockreasign = randperm(nSteps*2);
        Asigner = AsignBase(blockreasign);
        RefST = trialUnit(3,Asigner==0,i);
        TestST_T = trialUnit(3,Asigner==1,i);
        conRefST = conMSrt(Asigner==0,3,i);  % for testing
        incRefST = incMSrt(Asigner==0,3,i);
        conTestST_T = conMSrt(Asigner==1,3,i);
        incTestST_T = incMSrt(Asigner==1,3,i);

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

        testorder = datasample(1:nSteps,nSteps);  % bootstrapping with replacement
        for s = 1:nSteps
            %% Add block to test set on every step
            TestFL = [TestFL; TestFL_T(testorder(s))];
            TestPP = [TestPP; TestPP_T(testorder(s))];
            TestST = [TestST; TestST_T(testorder(s))];
            TestConFL = [TestConFL; conTestFL_T{testorder(s)}];  % for testing
            TestIncFL = [TestIncFL; incTestFL_T{testorder(s)}];
            TestConPP = [TestConPP; conTestPP_T{testorder(s)}];
            TestIncPP = [TestIncPP; incTestPP_T{testorder(s)}];
            TestConST = [TestConST; conTestST_T{testorder(s)}];
            TestIncST = [TestIncST; incTestST_T{testorder(s)}];
            conNT(k,s,1) = length(TestConFL);  % for testing and plotting
            incNT(k,s,1) = length(TestIncFL);
            conNT(k,s,2) = length(TestConPP);
            incNT(k,s,2) = length(TestIncPP);
            conNT(k,s,3) = length(TestConST);
            incNT(k,s,3) = length(TestIncST);

            %% Get the absolute difference
            absdiffvec(k,s,1) = abs(mean(RefFL)-mean(TestFL));
            absdiffvec(k,s,2) = abs(mean(RefPP)-mean(TestPP));
            absdiffvec(k,s,3) = abs(mean(RefST)-mean(TestST));
            absdiffvec_test(k,s,1) = abs((mean(cell2mat(incRefFL))-mean(cell2mat(conRefFL)))-(mean(TestIncFL)-mean(TestConFL)));  % for testing
            absdiffvec_test(k,s,2) = abs((mean(cell2mat(incRefPP))-mean(cell2mat(conRefPP)))-(mean(TestIncPP)-mean(TestConPP)));
            absdiffvec_test(k,s,3) = abs((mean(cell2mat(incRefST))-mean(cell2mat(conRefST)))-(mean(TestIncST)-mean(TestConST)));
        end
    end
    aProgMean(1,:,i) = mean(absdiffvec(:,:,1),1);
    aProgMean(2,:,i) = mean(absdiffvec(:,:,2),1);
    aProgMean(3,:,i) = mean(absdiffvec(:,:,3),1);
    aProgSD(1,:,i) = std(absdiffvec(:,:,1));
    aProgSD(2,:,i) = std(absdiffvec(:,:,2));
    aProgSD(3,:,i) = std(absdiffvec(:,:,3));
    aProgUB(1,:,i) = prctile(absdiffvec(:,:,1),UpB);
    aProgUB(2,:,i) = prctile(absdiffvec(:,:,2),UpB);
    aProgUB(3,:,i) = prctile(absdiffvec(:,:,3),UpB);
    aProgLB(1,:,i) = prctile(absdiffvec(:,:,1),LwB);
    aProgLB(2,:,i) = prctile(absdiffvec(:,:,2),LwB);
    aProgLB(3,:,i) = prctile(absdiffvec(:,:,3),LwB);
    aProgMean_test(1,:,i) = mean(absdiffvec_test(:,:,1),1);  % for testing
    aProgMean_test(2,:,i) = mean(absdiffvec_test(:,:,2),1);
    aProgMean_test(3,:,i) = mean(absdiffvec_test(:,:,3),1);
    conNTi(1,:,i) = mean(conNT(:,:,1));  % for testing and plotting
    incNTi(1,:,i) = mean(incNT(:,:,1));
    conNTi(2,:,i) = mean(conNT(:,:,2));
    incNTi(2,:,i) = mean(incNT(:,:,2));
    conNTi(3,:,i) = mean(conNT(:,:,3));
    incNTi(3,:,i) = mean(incNT(:,:,3));
end

% if and(rtAcc==1,practice==1)
%     aProgMean_1_1 = aProgMean;
%     save('RTstabilityCurveVariables1_initial_practice','aProgMean_1_1','aProgSD','aProgUB','aProgLB','conNTi','incNTi')
% elseif and(rtAcc==1,practice==0)
%     aProgMean_1_2 = aProgMean;
%     save('RTstabilityCurveVariables1_initial','aProgMean_1_2','aProgSD','aProgUB','aProgLB','conNTi','incNTi')
% elseif rtAcc == 2
%     aProgMean_2 = aProgMean;
%     save('ACCstabilityCurveVariables1_initial','aProgMean_2','aProgSD','aProgUB','aProgLB','conNTi','incNTi')
% elseif rtAcc == 3
%     aProgMean_3 = aProgMean;
%     save('IESstabilityCurveVariables1_initial_practice','aProgMean_3','aProgSD','aProgUB','aProgLB','conNTi','incNTi')
% end

% if and(rtAcc==1,practice==1)
%     load RTstabilityCurveVariables1_initial_practice
%     aProgMean = aProgMean_1_1;
% elseif and(rtAcc==1,practice==0)
%     load RTstabilityCurveVariables1_initial
%     aProgMean = aProgMean_1_2;
% elseif rtAcc == 2
%     load ACCstabilityCurveVariables1_initial
%     aProgMean = aProgMean_2;
% elseif rtAcc == 3
%     load IESstabilityCurveVariables1_initial_practice
%     aProgMean = aProgMean_3;
% end

%% Plot the data
% Stability curve
taskStrngfx = {'Flanker','Prime-Probe','Stroop'};
stp = 70;  % 72-2; nSteps
nihTB = 40;  % number of trials of NIH toolbox flanker task
Ebg = 96;  % number of trials of Eisenberg et al.'s (2019) Stroop task
for t = 1:nTask
    figure
    for i = 1:nSubj
        x = conNTi(t,1:stp,i)+incNTi(t,1:stp,i);  % # of trials of the test set sample
        h = plot(x,aProgMean(t,1:stp,i),'Color',cmap(i,:),'LineWidth',1.2); hold on
        set(get(h,'Parent'),'XScale','log')  % log-scale
    end
    xlabel('Number of trials','FontSize',14)
    xlim([0 4000])  % 18session*200trials = 3600
    if rtAcc == 1
        ylabel('Absolute difference from test half (ms)','FontSize',14)
        if t == 1
            ylim([0 25]);
        elseif t == 2
            ylim([0 35]);
        elseif t == 3
            ylim([0 55]);
        end
    elseif rtAcc == 2
        ylabel('Absolute difference from test half','FontSize',14)
        ylim([0 0.08]);
    elseif rtAcc == 3
        ylabel('Absolute difference from test half','FontSize',14)
        if t == 1
            ylim([0 40]);
        elseif t == 2
            ylim([0 60]);
        elseif t == 3
            ylim([0 80]);
        end
    end
    xticks([0 50 100 200 400 800 1600 3200])
    title([taskStrngfx{t} ' Task'],'FontSize',15)
    if x(1) < nihTB
        xline(nihTB,'k','LineWidth',0.7,'LineStyle','--')
    else
        xline(x(1),'k','LineWidth',0.7,'LineStyle','--')
    end
    xline(Ebg,'k','LineWidth',0.7,'LineStyle','-.')
    if t == 1
        legend('Participant 03','Participant 04','Participant 05','Participant 06',...
            'Participant 07','Participant 08','Participant 10','Participant 12',...
            'NIH toolbox','Eisenberg et al.','FontSize',12)
    end
    grid on
end