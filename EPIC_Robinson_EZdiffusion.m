%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
%
% To run this sciprt:
% You need Robinson & Steyvers' (2023) flanker task data (download at
% https://osf.io/6hjwv/),
% a Matlab function named "ezdiffusion.m" (download at https://osf.io/jk9nb/), 
% and ICC.m (download at
% https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
% What this script does:
% Plots the EZ-diffusion model parameters across varying numbers of trials
% (RT needs to be in SECS to run EZ-diffusion modeling)
%
% What this script outputs:
% Supp. Fig. 12A
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
nTdrawn = [50, 100, 200, 400, 800, 1600 3200]/2;  % number of trials drawn
drawLR = 5;  % ~nTdrawn

%% 1. Preprocess data
% (1) Exclude participants - Have consistent criteria across analyses
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

%% 2. EZ-diffusion modeling
HrtR1 = nan(nCond,drawLR,l_sGrp);  % mean RT (first half) in secs
HrtR2 = nan(nCond,drawLR,l_sGrp);
HrtvR1 = nan(nCond,drawLR,l_sGrp);  % RT variance
HrtvR2 = nan(nCond,drawLR,l_sGrp);
HaccR1 = nan(nCond,drawLR,l_sGrp);  % mean accuracy
HaccR2 = nan(nCond,drawLR,l_sGrp);
vMatR1 = nan(nCond,drawLR,l_sGrp);  % drift rate; test 1
aMatR1 = nan(nCond,drawLR,l_sGrp);  % boundary separation
TerMatR1 = nan(nCond,drawLR,l_sGrp);  % nondecision time
vMatR2 = nan(nCond,drawLR,l_sGrp);  % test 2
aMatR2 = nan(nCond,drawLR,l_sGrp);
TerMatR2 = nan(nCond,drawLR,l_sGrp);
for i = 1:l_sGrp
    cMr = cell2mat(conMSrt(1:stepTL(sGrp(i)),sGrp(i)));  % correct congruent
    iMr = cell2mat(incMSrt(1:stepTL(sGrp(i)),sGrp(i)));
    cMr = cMr./1000;  % convert to secs
    iMr = iMr./1000;
    cMa = cell2mat(conMSacc(1:stepTL(sGrp(i)),sGrp(i)));
    iMa = cell2mat(incMSacc(1:stepTL(sGrp(i)),sGrp(i)));
    % split-half option 1: odd / even
%     oddI = 1:2:length(cMr); evenI = 2:2:length(cMr);
%     cMr1 = cMr(oddI);
%     cMr2 = cMr(evenI);
%     oddI = 1:2:length(iMr); evenI = 2:2:length(iMr);
%     iMr1 = iMr(oddI);
%     iMr2 = iMr(evenI);
%     oddI = 1:2:length(cMa); evenI = 2:2:length(cMa);
%     cMa1 = cMa(oddI);
%     cMa2 = cMa(evenI);
%     oddI = 1:2:length(iMa); evenI = 2:2:length(iMa);
%     iMa1 = iMa(oddI);
%     iMa2 = iMa(evenI);
    % split-half option 2: half
    cMr1 = cMr(1:round(length(cMr)/2));
    cMr2 = cMr(round(length(cMr)/2)+1:end);
    iMr1 = iMr(1:round(length(iMr)/2));
    iMr2 = iMr(round(length(iMr)/2)+1:end);
    cMa1 = cMa(1:round(length(cMa)/2));
    cMa2 = cMa(round(length(cMa)/2)+1:end);
    iMa1 = iMa(1:round(length(iMa)/2));
    iMa2 = iMa(round(length(iMa)/2)+1:end);

    cMr1o = rmoutliers(cMr1,"mean");  % outliers: 3SD from mean
    cMr2o = rmoutliers(cMr2,"mean");
    iMr1o = rmoutliers(iMr1,"mean");
    iMr2o = rmoutliers(iMr2,"mean");
    for k = 1:drawLR
        HrtR1(1,k,i) = mean(cMr1o(1:nTdrawn(k)));
        HrtR1(2,k,i) = mean(iMr1o(1:nTdrawn(k)));
        HrtR2(1,k,i) = mean(cMr2o(1:nTdrawn(k)));
        HrtR2(2,k,i) = mean(iMr2o(1:nTdrawn(k)));
        HrtvR1(1,k,i) = var(cMr1o(1:nTdrawn(k)));
        HrtvR1(2,k,i) = var(iMr1o(1:nTdrawn(k)));
        HrtvR2(1,k,i) = var(cMr2o(1:nTdrawn(k)));
        HrtvR2(2,k,i) = var(iMr2o(1:nTdrawn(k)));
        HaccR1(1,k,i) = mean(cMa1(1:nTdrawn(k)));
        HaccR1(2,k,i) = mean(iMa1(1:nTdrawn(k)));
        HaccR2(1,k,i) = mean(cMa2(1:nTdrawn(k)));
        HaccR2(2,k,i) = mean(iMa2(1:nTdrawn(k)));
        % EZ-diffusion modeling
        [v1,a1,Ter1] = ezdiffusion(HaccR1(1,k,i),HrtvR1(1,k,i),HrtR1(1,k,i),nTdrawn(k));
        [v2,a2,Ter2] = ezdiffusion(HaccR1(2,k,i),HrtvR1(2,k,i),HrtR1(2,k,i),nTdrawn(k));
        vMatR1(1,k,i) = v1;
        vMatR1(2,k,i) = v2;
        aMatR1(1,k,i) = a1;
        aMatR1(2,k,i) = a2;
        TerMatR1(1,k,i) = Ter1;
        TerMatR1(2,k,i) = Ter2;
        clear v1 v2 a1 a2 Ter1 Ter2
        [v1,a1,Ter1] = ezdiffusion(HaccR2(1,k,i),HrtvR2(1,k,i),HrtR2(1,k,i),nTdrawn(k));
        [v2,a2,Ter2] = ezdiffusion(HaccR2(2,k,i),HrtvR2(2,k,i),HrtR2(2,k,i),nTdrawn(k));
        vMatR2(1,k,i) = v1;
        vMatR2(2,k,i) = v2;
        aMatR2(1,k,i) = a1;
        aMatR2(2,k,i) = a2;
        TerMatR2(1,k,i) = Ter1;
        TerMatR2(2,k,i) = Ter2;
        clear v1 v2 a1 a2 Ter1 Ter2
    end
end

%% Plot the results: EZD modeling parameters
cmap1 = parula(3);
cmap2 = spring(3);
cmap3 = cool(3);
ceTmp1 = squeeze(HrtR1(2,:,:)-HrtR1(1,:,:));  % inc-con
ceTmp2 = squeeze(HrtR2(2,:,:)-HrtR2(1,:,:));
ceTmpa1 = squeeze(HaccR1(1,:,:)-HaccR1(2,:,:));  % con-inc
ceTmpa2 = squeeze(HaccR2(1,:,:)-HaccR2(2,:,:));
vTmp1 = squeeze(vMatR1(1,:,:)-vMatR1(2,:,:));  % con-inc
vTmp2 = squeeze(vMatR2(1,:,:)-vMatR2(2,:,:));
aTmp1 = squeeze(aMatR1(1,:,:)-aMatR1(2,:,:));
aTmp2 = squeeze(aMatR2(1,:,:)-aMatR2(2,:,:));
TerTmp1  = squeeze(TerMatR1(2,:,:)-TerMatR1(1,:,:));  % inc-con
TerTmp2 = squeeze(TerMatR2(2,:,:)-TerMatR2(1,:,:));
figure
for k = 1:drawLR
    % 1) RT CE
    subplot(drawLR,5,5*(k-1)+1)
    scatter(ceTmp1(k,:),ceTmp2(k,:),8,'filled')
    %lsline
    xlim([min([ceTmp1(:);ceTmp2(:)]) max([ceTmp1(:);ceTmp2(:)])])
    ylim([min([ceTmp1(:);ceTmp2(:)]) max([ceTmp1(:);ceTmp2(:)])])
    [Icc] = ICC([ceTmp1(k,:)',ceTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    ylabel(['N = ' num2str(2*nTdrawn(k))])

    % 2) Acc CE
    subplot(drawLR,5,5*(k-1)+2)
    scatter(ceTmpa1(k,:),ceTmpa2(k,:),8,'filled')
    xlim([min([ceTmpa1(:);ceTmpa2(:)]) max([ceTmpa1(:);ceTmpa2(:)])])
    ylim([min([ceTmpa1(:);ceTmpa2(:)]) max([ceTmpa1(:);ceTmpa2(:)])])
    %lsline
    [Icc] = ICC([ceTmpa1(k,:)',ceTmpa2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square

    % 3) v
    subplot(drawLR,5,5*(k-1)+3)
    scatter(vTmp1(k,:),vTmp2(k,:),8,cmap1(1,:),'filled')
    xlim([min([vTmp1(:);vTmp2(:)]) max([vTmp1(:);vTmp2(:)])])
    ylim([min([vTmp1(:);vTmp2(:)]) max([vTmp1(:);vTmp2(:)])])
    %lsline
    [Icc] = ICC([vTmp1(k,:)',vTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square

    % 4) a
    subplot(drawLR,5,5*(k-1)+4)
    scatter(aTmp1(k,:),aTmp2(k,:),8,cmap1(2,:),'filled')
    xlim([min([aTmp1(:);aTmp2(:)]) max([aTmp1(:);aTmp2(:)])])
    ylim([min([aTmp1(:);aTmp2(:)]) max([aTmp1(:);aTmp2(:)])])
    %lsline
    [Icc] = ICC([aTmp1(k,:)',aTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square

    % 5) Ter
    subplot(drawLR,5,5*(k-1)+5)
    scatter(TerTmp1(k,:),TerTmp2(k,:),8,cmap1(3,:),'filled')
    xlim([min([TerTmp1(:);TerTmp2(:)]) max([TerTmp1(:);TerTmp2(:)])])
    ylim([min([TerTmp1(:);TerTmp2(:)]) max([TerTmp1(:);TerTmp2(:)])])
    %lsline
    [Icc] = ICC([TerTmp1(k,:)',TerTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc);
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
end
subplot(drawLR,5,1); title('CE RT')
subplot(drawLR,5,2); title('CE accuracy')
subplot(drawLR,5,3); title('Drift rate')
subplot(drawLR,5,4); title('Boundary separation')
subplot(drawLR,5,5); title('Nondecision time')
sgtitle('Robinson and Steyvers (2023)')