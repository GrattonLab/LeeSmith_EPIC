%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
% Stroop task
% Facilitation effect and interference effect
%
% To run this script:
% You need the EPIC dataset and "session_numbering.xlsx"
% (download both at https://osf.io/jk9nb)
%
% What this script does:
% Draws stability curves with method 1 (absolute difference between the
% reference set and the test set sample) for RT facilitation effect and
% interference effect
% Facilitation effect: Neutral-congruent
% Interference effect: incongruent-neutral
%
% What this script outputs:
% Stability curves of facilitation effect and interference effect
% Supp. Fig. 14
%
% Created on 12/20/2023 by HJ Lee
% Last modified on 12/20/2023
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
nTask = 1;
taskIndx = 4;  % 1=FL, 3=PP, 4=ST; session_numbering.xlsx
taskStrng = {'STROOP'};
nSession = 18;  % People do it for 9 weeks and each task is used during two sessions each week
nBlocks = 4;
nBlocksT = nSession*nBlocks;  % 72
bD = 2;  % break block into half
nTrials = 108/bD;  % half a block; this amount is added on each step
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
neuMSrt = cell(nSteps*2,nTask,nSubj);  % neutral condition
% (2) Acc
conMSacc = cell(nSteps*2,nTask,nSubj);
incMSacc = cell(nSteps*2,nTask,nSubj);
neuMSacc = cell(nSteps*2,nTask,nSubj);
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
            %if t == 3  % Stroop
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
            neuMSrt{(j-1)*nBlocks*bD+1,t,i} = 1000*T1.rt(and(T1.n0cong==3,T1.acc==1));
            conMSrt{(j-1)*nBlocks*bD+2,t,i} = 1000*T2.rt(and(T2.n0cong==1,T2.acc==1));
            incMSrt{(j-1)*nBlocks*bD+2,t,i} = 1000*T2.rt(and(T2.n0cong==2,T2.acc==1));
            neuMSrt{(j-1)*nBlocks*bD+2,t,i} = 1000*T2.rt(and(T2.n0cong==3,T2.acc==1));
            conMSrt{(j-1)*nBlocks*bD+3,t,i} = 1000*T3.rt(and(T3.n0cong==1,T3.acc==1));
            incMSrt{(j-1)*nBlocks*bD+3,t,i} = 1000*T3.rt(and(T3.n0cong==2,T3.acc==1));
            neuMSrt{(j-1)*nBlocks*bD+3,t,i} = 1000*T3.rt(and(T3.n0cong==3,T3.acc==1));
            conMSrt{(j-1)*nBlocks*bD+4,t,i} = 1000*T4.rt(and(T4.n0cong==1,T4.acc==1));
            incMSrt{(j-1)*nBlocks*bD+4,t,i} = 1000*T4.rt(and(T4.n0cong==2,T4.acc==1));
            neuMSrt{(j-1)*nBlocks*bD+4,t,i} = 1000*T4.rt(and(T4.n0cong==3,T4.acc==1));
            conMSrt{(j-1)*nBlocks*bD+5,t,i} = 1000*T5.rt(and(T5.n0cong==1,T5.acc==1));
            incMSrt{(j-1)*nBlocks*bD+5,t,i} = 1000*T5.rt(and(T5.n0cong==2,T5.acc==1));
            neuMSrt{(j-1)*nBlocks*bD+5,t,i} = 1000*T5.rt(and(T5.n0cong==3,T5.acc==1));
            conMSrt{(j-1)*nBlocks*bD+6,t,i} = 1000*T6.rt(and(T6.n0cong==1,T6.acc==1));
            incMSrt{(j-1)*nBlocks*bD+6,t,i} = 1000*T6.rt(and(T6.n0cong==2,T6.acc==1));
            neuMSrt{(j-1)*nBlocks*bD+6,t,i} = 1000*T6.rt(and(T6.n0cong==3,T6.acc==1));
            conMSrt{(j-1)*nBlocks*bD+7,t,i} = 1000*T7.rt(and(T7.n0cong==1,T7.acc==1));
            incMSrt{(j-1)*nBlocks*bD+7,t,i} = 1000*T7.rt(and(T7.n0cong==2,T7.acc==1));
            neuMSrt{(j-1)*nBlocks*bD+7,t,i} = 1000*T7.rt(and(T7.n0cong==3,T7.acc==1));
            conMSrt{(j-1)*nBlocks*bD+8,t,i} = 1000*T8.rt(and(T8.n0cong==1,T8.acc==1));
            incMSrt{(j-1)*nBlocks*bD+8,t,i} = 1000*T8.rt(and(T8.n0cong==2,T8.acc==1));
            neuMSrt{(j-1)*nBlocks*bD+8,t,i} = 1000*T8.rt(and(T8.n0cong==3,T8.acc==1));
            % Acc
            conMSacc{(j-1)*nBlocks*bD+1,t,i} = T1.acc(T1.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+1,t,i} = T1.acc(T1.n0cong==2);
            neuMSacc{(j-1)*nBlocks*bD+1,t,i} = T1.acc(T1.n0cong==3);
            conMSacc{(j-1)*nBlocks*bD+2,t,i} = T2.acc(T2.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+2,t,i} = T2.acc(T2.n0cong==2);
            neuMSacc{(j-1)*nBlocks*bD+2,t,i} = T2.acc(T2.n0cong==3);
            conMSacc{(j-1)*nBlocks*bD+3,t,i} = T3.acc(T3.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+3,t,i} = T3.acc(T3.n0cong==2);
            neuMSacc{(j-1)*nBlocks*bD+3,t,i} = T3.acc(T3.n0cong==3);
            conMSacc{(j-1)*nBlocks*bD+4,t,i} = T4.acc(T4.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+4,t,i} = T4.acc(T4.n0cong==2);
            neuMSacc{(j-1)*nBlocks*bD+4,t,i} = T4.acc(T4.n0cong==3);
            conMSacc{(j-1)*nBlocks*bD+5,t,i} = T5.acc(T5.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+5,t,i} = T5.acc(T5.n0cong==2);
            neuMSacc{(j-1)*nBlocks*bD+5,t,i} = T5.acc(T5.n0cong==3);
            conMSacc{(j-1)*nBlocks*bD+6,t,i} = T6.acc(T6.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+6,t,i} = T6.acc(T6.n0cong==2);
            neuMSacc{(j-1)*nBlocks*bD+6,t,i} = T6.acc(T6.n0cong==3);
            conMSacc{(j-1)*nBlocks*bD+7,t,i} = T7.acc(T7.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+7,t,i} = T7.acc(T7.n0cong==2);
            neuMSacc{(j-1)*nBlocks*bD+7,t,i} = T7.acc(T7.n0cong==3);
            conMSacc{(j-1)*nBlocks*bD+8,t,i} = T8.acc(T8.n0cong==1);
            incMSacc{(j-1)*nBlocks*bD+8,t,i} = T8.acc(T8.n0cong==2);
            neuMSacc{(j-1)*nBlocks*bD+8,t,i} = T8.acc(T8.n0cong==3);
        end
    end
end

%% Check for practice effects (if exists, remove with simple linear regression)
bFEmatRT = nan(nTask,nSteps*2,nSubj);  % before (nSteps==72)
bFEmatACC = nan(nTask,nSteps*2,nSubj);
bIEmatRT = nan(nTask,nSteps*2,nSubj);
bIEmatACC = nan(nTask,nSteps*2,nSubj);
for i = 1:nSubj
    for s = 1:nSteps*2
        bFEmatRT(1,s,i) = mean(neuMSrt{s,1,i})-mean(conMSrt{s,1,i});
        bFEmatACC(1,s,i) = mean(neuMSacc{s,1,i})-mean(conMSacc{s,1,i});
        bIEmatRT(1,s,i) = mean(incMSrt{s,1,i})-mean(neuMSrt{s,1,i});
        bIEmatACC(1,s,i) = mean(incMSacc{s,1,i})-mean(neuMSacc{s,1,i});
    end
end

% Plot
cmap = turbo(nSubj);
figure
subplot(2,2,1)
for i = 1:nSubj
    plot(1:nSteps*2,bFEmatRT(1,:,i),'Color',cmap(i,:),'LineWidth',1.2); hold on
end
xlabel('Trial units')
ylabel('Reaction time (ms)')
title('Facilitation Effect')
subplot(2,2,2)
for i = 1:nSubj
    plot(1:nSteps*2,bFEmatACC(1,:,i),'Color',cmap(i,:),'LineWidth',1.2); hold on
end
xlabel('Trial units')
ylabel('Accuracy')
title('Facilitation Effect')
subplot(2,2,3)
for i = 1:nSubj
    plot(1:nSteps*2,bIEmatRT(1,:,i),'Color',cmap(i,:),'LineWidth',1.2); hold on
end
xlabel('Trial units')
ylabel('Reaction time (ms)')
title('Interference Effect')
subplot(2,2,4)
for i = 1:nSubj
    plot(1:nSteps*2,bIEmatACC(1,:,i),'Color',cmap(i,:),'LineWidth',1.2); hold on
end
xlabel('Trial units')
ylabel('Accuracy')
title('Interference Effect')

%% Remove the initial 4 steps due to large variance in the beginning of experiment
rmvL = 1:4;  % nTrials*4 = the amount removing 216(ST)
conMSrt(rmvL,:,:) = [];
incMSrt(rmvL,:,:) = [];
neuMSrt(rmvL,:,:) = [];
conMSacc(rmvL,:,:) = [];
incMSacc(rmvL,:,:) = [];
neuMSacc(rmvL,:,:) = [];
bFEmatRT(:,rmvL,:) = [];
bFEmatACC(:,rmvL,:) = [];
bIEmatRT(:,rmvL,:) = [];
bIEmatACC(:,rmvL,:) = [];
nSteps = nSteps-(length(rmvL))/2;

%% Matrix definition and preassignment
aProgMean_FE = nan(nTask,nSteps,nSubj);
aProgMean_IE = nan(nTask,nSteps,nSubj);
conNTi = nan(nTask,nSteps,nSubj);
incNTi = nan(nTask,nSteps,nSubj);
neuNTi = nan(nTask,nSteps,nSubj);

%% RT or accuracy? Run this script TWICE from this line; change rtAcc
rtAcc = 1;  % rt:1, acc:2
if rtAcc == 1
    conMS = conMSrt;
    incMS = incMSrt;
    neuMS = neuMSrt;
elseif rtAcc == 2
    conMS = conMSacc;
    incMS = incMSacc;
    neuMS = neuMSacc;
end

%% Bootstrapping
AsignBase = [zeros(1,nSteps),ones(1,nSteps)];  % 0s for the reference and 1s for the test set
for i = 1:nSubj
    absdiffvec_FE = nan(numTest,nSteps);
    absdiffvec_IE = nan(numTest,nSteps);
    conNT = nan(numTest,nSteps);
    incNT = nan(numTest,nSteps);
    neuNT = nan(numTest,nSteps);
    for k = 1:numTest
        %% Make reference set & test set: Trial unit assignment in random orders
        blockreasign = randperm(nSteps*2);
        Asigner = AsignBase(blockreasign);
        conRef = conMS(Asigner==0,1,i);  % reference set
        incRef = incMS(Asigner==0,1,i);
        neuRef = neuMS(Asigner==0,1,i);
        conTest = conMS(Asigner==1,1,i);  % test set
        incTest = incMS(Asigner==1,1,i);
        neuTest = neuMS(Asigner==1,1,i);

        conTestsm = [];  % Test set sample
        incTestsm = [];
        neuTestsm = [];

        testorder = datasample(1:nSteps,nSteps);
        for s = 1:nSteps
            conTestsm = [conTestsm; conTest{testorder(s)}];
            incTestsm = [incTestsm; incTest{testorder(s)}];
            neuTestsm = [neuTestsm; neuTest{testorder(s)}];
            conNT(k,s) = length(conTestsm);
            incNT(k,s) = length(incTestsm);
            neuNT(k,s) = length(neuTestsm);

            %% Absolute difference
            FEref = mean(cell2mat(neuRef))-mean(cell2mat(conRef));
            FEtest = mean(neuTestsm)-mean(conTestsm);
            absdiffvec_FE(k,s) = abs(FEref-FEtest);
            IEref = mean(cell2mat(incRef))-mean(cell2mat(neuRef));
            IEtest = mean(incTestsm)-mean(neuTestsm);
            absdiffvec_IE(k,s) = abs(IEref-IEtest);
        end
    end
    aProgMean_FE(1,:,i) = mean(absdiffvec_FE);
    aProgMean_IE(1,:,i) = mean(absdiffvec_IE);
    conNTi(1,:,i) = mean(conNT);
    incNTi(1,:,i) = mean(incNT);
    neuNTi(1,:,i) = mean(neuNT);
end

if rtAcc == 1
    save('RTstabilityCurveVariables1_Stroop_FEIE','aProgMean_FE','aProgMean_IE','conNTi','incNTi','neuNTi');
elseif rtAcc == 2
    save('ACCstabilityCurveVariables1_Stroop_FEIE','aProgMean_FE','aProgMean_IE','conNTi','incNTi','neuNTi');
end

%if rtAcc == 1
%    load RTstabilityCurveVariables1_Stroop_FEIE
%elseif rtACC == 2
%    load ACCstabilityCurveVariables1_Stroop_FEIE
%end

%% Plot the data
% Stability curve
taskStrng_print = {'Stroop'};
stp = nSteps;
% Facilitation effect
figure
for i = 1:nSubj
    x = neuNTi(1,1:stp,i)+conNTi(1,1:stp,i);
    h = plot(x,aProgMean_FE(1,1:stp,i),'Color',cmap(i,:),'LineWidth',1.8); hold on
    set(get(h,'Parent'),'XScale','log')
end
set(gca,'FontSize',15)
xlabel('Number of trials','FontSize',17)
xlim([0 4000])  % 18session*200trials = 3600
if rtAcc == 1
    ylabel('Absolute difference from test half (ms)','FontSize',17)
    ylim([0 55]);
elseif rtAcc == 2
    ylabel('Absolute difference from test half','FontSize',17)
    ylim([0 0.08]);
end
xticks([0 50 100 200 400 800 1600 3200])
title('Facilitation Effect','FontSize',19)
grid on
%legend('Participant 03','Participant 04','Participant 05','Participant 06',...
%    'Participant 07','Participant 08','Participant 10','Participant 12','FontSize',12)
% Interference effect
figure
for i = 1:nSubj
    x = incNTi(1,1:stp,i)+neuNTi(1,1:stp,i);
    h = plot(x,aProgMean_IE(1,1:stp,i),'Color',cmap(i,:),'LineWidth',1.8); hold on
    set(get(h,'Parent'),'XScale','log')
end
set(gca,'FontSize',15)
xlabel('Number of trials','FontSize',17)
xlim([0 4000])  % 18session*200trials = 3600
if rtAcc == 1
    ylabel('Absolute difference from test half (ms)','FontSize',17)
    ylim([0 55]);
elseif rtAcc == 2
    ylabel('Absolute difference from test half','FontSize',17)
    ylim([0 0.08]);
end
xticks([0 50 100 200 400 800 1600 3200])
title('Interference Effect','FontSize',19)
grid on
%legend('Participant 03','Participant 04','Participant 05','Participant 06',...
%    'Participant 07','Participant 08','Participant 10','Participant 12','FontSize',12)