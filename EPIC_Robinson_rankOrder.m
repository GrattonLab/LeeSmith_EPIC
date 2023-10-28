%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
%
% To run this sciprt:
% You need Robinson & Steyvers' (2023) flanker task data (download at
% https://osf.io/6hjwv/)
% You also need ICC.m (download at
% https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
% What this script does:
% Plots the rank correlation coefficients between CE and incongruent trial
%
% Exclusion criteria:
% (a) Below 70% accuracy in either experimental conditions (congruent/incongruent)
% (b) Sessions with 0% accuracy
% (c) Below 2500 trials
%
% What this script outputs:
% Supp. Fig. 6A
%
% Created on 05/16/2023 by HJ Lee
% Last modified on 07/03/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')

%% Load Robinson & Steyvers' (2023) flanker task data
load FlankerData_learn

%% Parameter settings
tmpnSubj = length(d_flanker_tab);  % number of participants: 495
nCond = 2;  % number of experimental conditions: congruent/incongruent

%% 1. Preprocess data
% Exclude participants - Have consistent criteria across analyses
% Matrix definition and preassignment
nSess = nan(tmpnSubj,1);  % number of sessions each participant had (because this also differs across participants)
for i = 1:tmpnSubj
    nSess(i) = length(find(d_flanker_tab{i,1}.trial_num==1));
end
[Max,Imax] = max(nSess);
nStep = Max;  % to store them in a unified size of cell
conMSrt = cell(nStep,tmpnSubj);
incMSrt = cell(nStep,tmpnSubj);
conMSacc = cell(nStep,tmpnSubj);
incMSacc = cell(nStep,tmpnSubj);

% Break up data based on session and store them in matrices
for i = 1:tmpnSubj
    tmpT = table(d_flanker_tab{i,1}.trial_num,double(d_flanker_tab{i,1}.compatible),...
        d_flanker_tab{i,1}.accuracy,d_flanker_tab{i,1}.response_time,'VariableNames',["trial_num","n0cong","acc","rt"]);
    indx = find(tmpT.trial_num==1);
    indx = [indx; length(table2array(tmpT))+1];

    for b = 1:length(indx)-1
        tmpTd = tmpT(indx(b):indx(b+1)-1,:);
        conMSrt{b,i} = tmpTd.rt(and(tmpTd.n0cong==1,tmpTd.acc==1));  % congruent RT
        incMSrt{b,i} = tmpTd.rt(and(tmpTd.n0cong==0,tmpTd.acc==1));  % incongruent
        conMSacc{b,i} = tmpTd.acc(tmpTd.n0cong==1);  % accuracy
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
nSubj = length(subjID);  % 448

%% Exclusion criteria (c): Select participants with more than 2500 CORRECT trials
lrg = 2500;
nTrialMatexcld = nan(nSubj,1);
for i = 1:nSubj
    nTrialMatexcld(i) = length(cell2mat(conMSrt(1:stepTL(i),i)))+length(cell2mat(incMSrt(1:stepTL(i),i)));
end
sGrp = find(nTrialMatexcld>lrg);
l_sGrp = length(sGrp);  % 185

%% 2. Descriptive
GrtMatR = nan(nCond,l_sGrp);  % grand mean RT
GaccMatR = nan(nCond,l_sGrp);  % accuracy
for i = 1:l_sGrp
    cMr = cell2mat(conMSrt(1:stepTL(sGrp(i)),sGrp(i)));  % correct congruent
    iMr = cell2mat(incMSrt(1:stepTL(sGrp(i)),sGrp(i)));
    cMr = cMr./1000;  % convert to secs
    iMr = iMr./1000;
    GrtMatR(1,i) = mean(rmoutliers(cMr,"mean"));  % outliers: more than 3 SD from mean
    GrtMatR(2,i) = mean(rmoutliers(iMr,"mean"));
    cMa = cell2mat(conMSacc(1:stepTL(sGrp(i)),sGrp(i)));
    iMa = cell2mat(incMSacc(1:stepTL(sGrp(i)),sGrp(i)));
    GaccMatR(1,i) = mean(cMa);
    GaccMatR(2,i) = mean(iMa);
end

%% Plot the correlation between CE and incongruent trials
% (1) RT
[~,p1] = sort(GrtMatR(2,:)-GrtMatR(1,:),'descend');
r1 = 1:l_sGrp;
r1(p1) = r1;
[~,p2] = sort(GrtMatR(2,:),'descend');
r2 = 1:l_sGrp;
r2(p2) = r2;
% Kendall's tau and ICC
tauH = corr(r1',r2','type','kendall');
iccH = ICC([r1',r2'],'A-k',0.05);
figure
scatter(GrtMatR(2,:),GrtMatR(2,:)-GrtMatR(1,:),'filled')
lsline
% refline
axis square
set(gca,'FontSize',13)
str = sprintf('  Tau = %1.2f\n  ICC = %1.2f', tauH, iccH);
%str = sprtinf('  Tau = %1.2f', tauH);
rText = text(min(get(gca,'xlim')),max(get(gca,'ylim')), str);
set(rText, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
xlabel('Incongruent trial RT (sec)','FontSize',14)
ylabel('CE RT (sec)','FontSize',14)
title('Robinson & Steyvers'' (2023) data','FontSize',15)
% (2) Accuracy
GpeMatR = 100.*(1-GaccMatR);
[~,p1] = sort(GpeMatR(2,:)-GpeMatR(1,:),'descend');
r1 = 1:l_sGrp;
r1(p1) = r1;
[~,p2] = sort(GpeMatR(2,:),'descend');
r2 = 1:l_sGrp;
r2(p2) = r2;
% Kendall's tau and ICC
tauH = corr(r1',r2','type','kendall');
iccH = ICC([r1',r2'],'A-k',0.05);
figure
scatter(GpeMatR(2,:),GpeMatR(2,:)-GpeMatR(1,:),'filled')
lsline
% refline
axis square
%ylim([-0.02 0.15])  % accuracy
set(gca,'FontSize',13)
str = sprintf('  Tau = %1.2f\n  ICC = %1.2f', tauH, iccH);
%str = sprintf('  Tau = %1.2f', tauH);
rText = text(min(get(gca,'xlim')),max(get(gca,'ylim')), str);
set(rText, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
xlabel('Incongruent trial percent correct(%)','FontSize',14)
ylabel('CE percent correct (%)','FontSize',14)
title('Robinson & Steyvers'' (2023) data','FontSize',15)