%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC - Public Data 1
%
%% Instructions for running the script:
% You will need Robinson & Steyvers' (2023) flanker task data, available
% for download at https://osf.io/6hjwv/
%
%% Purpose of the script:
% This script plots the corrected between-subject standard deviation by
% adjusting for sampling noise.
% Sample variance = 2*(Mean squared error)/trial size + true variance
% 
%% Outputs:
% Supp. Fig. 18
%
% Created on 01/25/2024 by HJ Lee
% Last modified on 02/04/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')

%% Load Robinson & Steyvers' (2023) flanker task data
load FlankerData_learn

%% Parameters
tmpnSubj = length(d_flanker_tab);  % number of total participants: 495
nCond = 2;  % congruent(1), incongruent(2)
nStd = 3;  % standard deviation from mean for outlier setting
UpB = 97.5;  % upper percentile of 95% confidence interval
LwB = 2.5;  % lower percentile
numTest = 1000;  % number of simulated testing

%% 1. Preprocess data
% Matrix definition and preassignment
nSess = nan(tmpnSubj,1);  % number of sessions each participant had (this differs across participants)
for i = 1:tmpnSubj
    nSess(i) = length(find(d_flanker_tab{i,1}.trial_num==1));
end
[Max,Imax] = max(nSess);
nStep = Max;  % to store them in a unified size of cell; this is the # of steps because one session data (~64 trials) were added to test set sample on each step
conMSrt = cell(nStep,tmpnSubj);
incMSrt = cell(nStep,tmpnSubj);
conMSacc = cell(nStep,tmpnSubj);
incMSacc = cell(nStep,tmpnSubj);

% Break up data based on session and store them in matrices
ageBin = nan(tmpnSubj,1);  % Each participant age
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
% RT matrix is empty when all responses are incorrect
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

%% Calculate standard deviation (SD) and mean squared error (MSE)
nSampling = [50 80 100 200 400 800 1600]./nCond;  % number of trial sampling
% Preassignment
imCE = nan(l_sGrp,length(nSampling),numTest);  % mean congruency effect of each participant
imC = nan(l_sGrp,length(nSampling),numTest);  % congruent trial
imI = nan(l_sGrp,length(nSampling),numTest);  % incongruent trial
fxdC = cell(l_sGrp,length(nSampling),numTest);
fxdI = cell(l_sGrp,length(nSampling),numTest);
for i_t = 1:numTest
    for i = 1:l_sGrp
        for i_s = 1:length(nSampling)
            tmpC = datasample(RTdistC{i},nSampling(i_s));
            tmpI = datasample(RTdistI{i},nSampling(i_s));
            fxdC{i,i_s,i_t} = tmpC;
            fxdI{i,i_s,i_t} = tmpI;
            imCE(i,i_s,i_t) = mean(tmpI)-mean(tmpC);
            imC(i,i_s,i_t) = mean(tmpC);
            imI(i,i_s,i_t) = mean(tmpI);
        end
    end
end
mCE = squeeze(mean(imCE));  % size: nSampling*nTesting
sdCE = squeeze(std(imCE));  % compare it with mean(sqrt(Vs),2); should be the same
% Preassignment
Vs = nan(length(nSampling),numTest);  % Sample variance
iSamVar = nan(l_sGrp,length(nSampling),numTest);  % preassignment for individual sample variance
myMSE = nan(length(nSampling),numTest);
c_Vs = nan(length(nSampling),numTest);
for i_t = 1:numTest
    iMSEc = zeros(nSampling(end),l_sGrp,length(nSampling));  % this may not be filled up as trial size varies across participants but it's okay as you're adding the valules (empty ones=zero)
    iMSEi = zeros(nSampling(end),l_sGrp,length(nSampling));
    for i_s = 1:length(nSampling)
        for i = 1:l_sGrp
            iSamVar(i,i_s,i_t) = (imCE(i,i_s,i_t)-mCE(i_s,i_t))^2;
            for l1 = 1:length(fxdC{i,i_s,i_t})
                iMSEc(l1,i,i_s) = (fxdC{i,i_s,i_t}(l1)-imC(i,i_s,i_t))^2;
            end
            for l2 = 1:length(fxdI{i,i_s,i_t})
                iMSEi(l2,i,i_s) = (fxdI{i,i_s,i_t}(l2)-imI(i,i_s,i_t))^2;
            end
        end
        Vs(i_s,i_t) = sum(iSamVar(:,i_s,i_t))/(l_sGrp-1);
        % check if sqrt(Vs) and sdCE are similar if not the same
        L = nSampling(i_s);  % number of trials PER condition
        tmpa = iMSEc(:,:,i_s);
        tmpa = tmpa(:);
        tmpb = iMSEi(:,:,i_s);
        tmpb = tmpb(:);
        myMSE(i_s,i_t) = (sum(tmpa)+sum(tmpb))/(2*l_sGrp*(L-1));
        mySNR = Vs(i_s,i_t)/myMSE(i_s,i_t)-2/L;
        c_Vs(i_s,i_t) = Vs(i_s,i_t)-2*myMSE(i_s,i_t)/L;
    end
end
save('Robinson_sim_correctedVariance.mat','imCE','imC','imI','fxdC','fxdI','Vs','iSamVar','myMSE','c_Vs')

%% Plot the results (Supp. Fig. 18)
figure
subplot(2,1,1)
x = 1:length(nSampling);
xconf = [x x(end:-1:1)];
myMSE_95CI = (prctile(myMSE',UpB)-prctile(myMSE',LwB))/2;
CIl = mean(myMSE,2)'-myMSE_95CI;
CIu = mean(myMSE,2)'+myMSE_95CI;
u = CIu;
yconf = [CIl u(end:-1:1)];
p = fill(xconf,yconf,'r');
p.FaceColor = [1 0.8 0.8];
p.EdgeColor = 'none';
hold on
plot(x,mean(myMSE,2),'Color','r','LineWidth',2)
set(gca,'FontSize',15)
grid on
xlabel('Number of trials','FontSize',15)
xticks(1:length(nSampling))
xticklabels(nSampling.*nCond)
ylabel('msec','FontSize',15)
title('A) Mean Squared Error','FontSize',15)
ylim([floor(CIl(1)) ceil(CIu(1))])
legend('95% confidence interval','Mean','FontSize',13)
subplot(2,1,2)
% Sample between-subject standard deviation
SD_95CI = (prctile(sqrt(Vs)',UpB)-prctile(sqrt(Vs)',LwB))/2;
CIl = mean(sqrt(Vs),2)'-SD_95CI;
CIu = mean(sqrt(Vs),2)'+SD_95CI;
tmp = CIu;
u = CIu;
yconf = [CIl u(end:-1:1)];
p = fill(xconf,yconf,'k');
p.FaceColor = [0.8 0.8 0.8];
p.EdgeColor = 'none';
hold on
plot(x,mean(sqrt(Vs),2),'Color','k','LineWidth',2)
% Corrected between-subject standard deviation
c_sdCE_95CI = (prctile(real(sqrt(c_Vs))',UpB)-prctile(real(sqrt(c_Vs))',LwB))/2;
CIl = mean(sqrt(c_Vs),2)'-c_sdCE_95CI;
CIu = mean(sqrt(c_Vs),2)'+c_sdCE_95CI;
u = CIu;
yconf = [CIl u(end:-1:1)];
p = fill(xconf,yconf,'b');
p.FaceColor = [0.8 0.8 1];
p.EdgeColor = 'none';
hold on
plot(x,mean(sqrt(c_Vs),2),'Color','b','LineWidth',2)
set(gca,'FontSize',15)
grid on
xlabel('Number of trials','FontSize',15)
xticks(1:length(nSampling))
xticklabels(nSampling.*nCond)
ylabel('msec','FontSize',15)
ylim([real(floor(CIl(1))) ceil(tmp(1))])
legend('','Sample','','Corrected','FontSize',13)
title('B) Between-Subject Standard Deviation','FontSize',15)