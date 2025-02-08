%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC - Flanker, Prime-probe, and Stroop tasks
%
%% Instructions for running the script:
% You will need the EPIC dataset and the "session_numbering.xlsx" file, 
% both available for download at https://osf.io/jk9nb.
% You also need ICC.m (download at https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc).
%
%% Purpose of the script:
% This script generates stability curves using Method 1, but instead of
% calculating the absolute difference between the reference set and the
% test set sample, it calculates the ICC across participants (note our
% dataset has only 9 participants)
%
%% Outputs:
% Supp. Fig. 15
%
% Created on 08/28/2024 by HJ Lee
% Last modified on 02/04/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')

%% Experimental and statistical parameters
% Participant information
tmpSubjID = 1:12;
subjExcld = [1,2,11,9];  % EPIC 9 excluded
subjID = setdiff(tmpSubjID, subjExcld);
nSubj = length(subjID);

% Task information
nTask = 3;  % 1: Flanker, 2: Prime-Probe, 3: Stroop
taskIndx = [1,3,4];  % 1:Flanker, 2:GoNogo, 3:PrimeProbe, 4:STROOP
taskStrng = {'Flanker','PrimeProbe','STROOP'};
nSession = 18;  % number of sessions for each task
nBlocks = 4;
nBlocksT = nSession*nBlocks;  % 72
bD = 2;  % break block into half
nTrials = [100; 96; 108]./bD;  % this amount is added on each step
nSteps = (nBlocksT*bD)/2;  % number of steps; add half the block to the test set on every step
UpB = 97.5;  % upper percentile of 95% confidence interval
LwB = 2.5;  % lower percentile
numTest = 1000;  % number of repetition
%nCond = 2;  % congruent(1) vs. incongruent(2)

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
        dir = ['rawData/EPIC' num2str(subjID(i))];
        %% Load session indexing information
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

%% Regress the practice effects - RT
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

%% Remove the initial 4 steps to get rid of variance in the beginning of the task
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
conNTi = nan(nTask,nSteps);  % number of trials - congruent condition; for plotting
incNTi = nan(nTask,nSteps);  % incongruent

%% RT or accuracy? Run this script FOUR times from this line; change rtAcc and practice
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

%% Calculate ICC
myICC = nan(nSteps,numTest,nTask);
AsignBase = [zeros(1,nSteps),ones(1,nSteps)];
conNT = nan(numTest,nSteps,nTask);
incNT = conNT;
for k = 1:numTest
    ICCbowl_r = nan(nSubj,nSteps,nTask);
    ICCbowl_t = ICCbowl_r;
    for i = 1:nSubj
        %% Make reference set & test set: Trial unit assignment in random orders
        % Randomize the trial units that are added for each step (Test Set)
        % Flanker
        blockreasign = randperm(nSteps*2); % scramble the order
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

        % Make test set storage
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
            %% Record data
            ICCbowl_r(i,s,1) = mean(RefFL);  % should have the same data across steps
            ICCbowl_t(i,s,1) = mean(TestFL);
            ICCbowl_r(i,s,2) = mean(RefPP);
            ICCbowl_t(i,s,2) = mean(TestPP);
            ICCbowl_r(i,s,3) = mean(RefST);
            ICCbowl_t(i,s,3) = mean(TestST);

            if i == 1
                conNT(k,s,1) = length(TestConFL);  % for testing and plotting
                incNT(k,s,1) = length(TestIncFL);
                conNT(k,s,2) = length(TestConPP);
                incNT(k,s,2) = length(TestIncPP);
                conNT(k,s,3) = length(TestConST);
                incNT(k,s,3) = length(TestIncST);
            end
        end
    end
    for s = 1:nSteps
        myICC(s,k,1) = ICC([ICCbowl_r(:,s,1),ICCbowl_t(:,s,1)],'C-1',0.05);  % fl
        myICC(s,k,2) = ICC([ICCbowl_r(:,s,2),ICCbowl_t(:,s,2)],'C-1',0.05);  % pp
        myICC(s,k,3) = ICC([ICCbowl_r(:,s,3),ICCbowl_t(:,s,3)],'C-1',0.05);  % st
    end
end
% Average
plot_myICC = squeeze(mean(myICC,2));
% 95% confidence interval
plot_myICC_ci = squeeze((prctile(myICC,UpB,2)-prctile(myICC,LwB,2))/2);
CIl = plot_myICC-plot_myICC_ci;
CIu = plot_myICC+plot_myICC_ci;

% Trial number
conNTi(1,:) = mean(conNT(:,:,1));
incNTi(1,:) = mean(incNT(:,:,1));
conNTi(2,:) = mean(conNT(:,:,2));
incNTi(2,:) = mean(incNT(:,:,2));
conNTi(3,:) = mean(conNT(:,:,3));
incNTi(3,:) = mean(incNT(:,:,3));

%% Plot results
taskStrngfx = {'Flanker','Prime-Probe','Stroop'};
stp = nSteps;
palette = parula(nTask);
palette(3,:) = [1 0 0];
% Flanker task
figure
x = conNTi(1,1:stp)+incNTi(1,1:stp);
% confidence interval
xconf = [x x(end:-1:1)];
u = CIu(:,1);
yconf = [CIl(:,1); u(end:-1:1)]';
p = fill(xconf,yconf,'k'); hold on
p.FaceColor = [0.6 0.5 0.9];
p.EdgeColor = 'none';
% mean
h1 = plot(x,plot_myICC(:,1)','Color',palette(1,:),'LineWidth',5);
set(get(h1,'Parent'),'XScale','log')  % log-scale
set(gca,'FontSize',16)
xlabel('Number of trials','FontSize',24)
xlim([50 4000])
%ylabel('ICC','FontSize',28)
ylim([0 1.1])
xticks([0 50 100 200 400 800 1600 3200])
grid on

% Prime-probe task
figure
x = conNTi(2,1:stp)+incNTi(2,1:stp);
% confidence interval
xconf = [x x(end:-1:1)];
u = CIu(:,2);
yconf = [CIl(:,2); u(end:-1:1)]';
p = fill(xconf,yconf,'r'); hold on
p.FaceColor = [0.7 0.95 0.95];
p.EdgeColor = 'none';
% mean
h2 = plot(x,plot_myICC(:,2)','Color',palette(2,:),'LineWidth',5);
set(get(h2,'Parent'),'XScale','log')  % log-scale
set(gca,'FontSize',16)
xlabel('Number of trials','FontSize',24)
xlim([50 4000])
%ylabel('ICC','FontSize',28)
ylim([0 1.1])
xticks([0 50 100 200 400 800 1600 3200])
grid on

% Stroop task
figure
x = conNTi(3,1:stp)+incNTi(3,1:stp);
% confidence interval
xconf = [x x(end:-1:1)];
u = CIu(:,3);
yconf = [CIl(:,3); u(end:-1:1)]';
p = fill(xconf,yconf,'b'); hold on
p.FaceColor = [1 0.8 0.8];
p.EdgeColor = 'none';
h3 = plot(x,plot_myICC(:,3)','Color',palette(3,:),'LineWidth',5);
set(get(h3,'Parent'),'XScale','log')  % log-scale
set(gca,'FontSize',16)
xlabel('Number of trials','FontSize',24)
xlim([50 4000])
%ylabel('ICC','FontSize',28)
ylim([0 1.1])
xticks([0 50 100 200 400 800 1600 3200])
grid on