%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
%
% To run this sciprt:
% You need Robinson & Steyvers' (2023) flanker task data (download at
% https://osf.io/6hjwv/),
% a Matlab function file named "ezdiffusion.m" (download at
% https://osf.io/jk9nb), and ICC.m (download at
% https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
% What this script does:
% Plots the EZ-diffusion model parameters across varying numbers of trials
% (RT needs to be in SECS to run EZ-diffusion modeling)
%
% What this script outputs:
% Fig. 8 and Supp. Fig. 19B
%
% Created on 05/16/2023 by HJ Lee
% Last modified on 01/18/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')

%% Load Robinson & Steyvers' (2023) flanker task data
load FlankerData_learn

%% Parameter settings
nCond = 2;  % number of experimental conditions: congruent, incongruent
UpB = 97.5;
LwB = 2.5;
nBt = 1000;  % number of bootstrapping
tmpnSubj = length(d_flanker_tab);  % number of total participants: 495
nTdrawn = [50, 100, 200, 400, 800, 1600 3200]/2;  % number of trials drawn; divided by 2 as sampling is done for con/inc separately
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
% Descriptive stats
HrtR1 = nan(nCond,drawLR,l_sGrp);  % mean RT (first half) in secs; first half
HrtR2 = nan(nCond,drawLR,l_sGrp);  % second half
HrtvR1 = nan(nCond,drawLR,l_sGrp);  % RT variance
HrtvR2 = nan(nCond,drawLR,l_sGrp);
HaccR1 = nan(nCond,drawLR,l_sGrp);  % mean accuracy
HaccR2 = nan(nCond,drawLR,l_sGrp);
% EZ-diffusion model parameters
vMatR1 = nan(nCond,drawLR,l_sGrp);  % drift rate; test 1
aMatR1 = nan(nCond,drawLR,l_sGrp);  % boundary separation
TerMatR1 = nan(nCond,drawLR,l_sGrp);  % nondecision time
vMatR2 = nan(nCond,drawLR,l_sGrp);  % test 2
aMatR2 = nan(nCond,drawLR,l_sGrp);
TerMatR2 = nan(nCond,drawLR,l_sGrp);
% Standard error (95% confidence interval of the mean)
vMatRse = nan(2,drawLR,l_sGrp);  % 2=first & second half
aMatRse = nan(2,drawLR,l_sGrp);
TerMatRse = nan(2,drawLR,l_sGrp);
% ICC
vICC1 = nan(drawLR,nBt,l_sGrp);  % for calculating 95% confidence interval of ICC
vICC2 = nan(drawLR,nBt,l_sGrp);
aICC1 = nan(drawLR,nBt,l_sGrp);
aICC2 = nan(drawLR,nBt,l_sGrp);
TerICC1 = nan(drawLR,nBt,l_sGrp);
TerICC2 = nan(drawLR,nBt,l_sGrp);
for i = 1:l_sGrp
    % preassignment for bootstrapping
    rt_bt1 = nan(nCond,drawLR,nBt);
    rt_bt2 = nan(nCond,drawLR,nBt);
    rtv_bt1 = nan(nCond,drawLR,nBt);
    rtv_bt2 = nan(nCond,drawLR,nBt);
    acc_bt1 = nan(nCond,drawLR,nBt);
    acc_bt2 = nan(nCond,drawLR,nBt);

    vSE1 = nan(nCond,drawLR,nBt);  % 95% confidence interval
    vSE2 = nan(nCond,drawLR,nBt);
    aSE1 = nan(nCond,drawLR,nBt);
    aSE2 = nan(nCond,drawLR,nBt);
    TerSE1 = nan(nCond,drawLR,nBt);
    TerSE2 = nan(nCond,drawLR,nBt);

    % Preprocess Robinson and Steyvers' data
    cMr = cell2mat(conMSrt(1:stepTL(sGrp(i)),sGrp(i)));  % correct congruent trials RT
    iMr = cell2mat(incMSrt(1:stepTL(sGrp(i)),sGrp(i)));
    cMr = cMr./1000;  % convert to secs
    iMr = iMr./1000;
    cMa = cell2mat(conMSacc(1:stepTL(sGrp(i)),sGrp(i)));  % accuracy
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
        [v1,a1,Ter1] = ezdiffusion(HaccR1(1,k,i),HrtvR1(1,k,i),HrtR1(1,k,i),nTdrawn(k));  % first half
        [v2,a2,Ter2] = ezdiffusion(HaccR1(2,k,i),HrtvR1(2,k,i),HrtR1(2,k,i),nTdrawn(k));
        vMatR1(1,k,i) = v1;
        vMatR1(2,k,i) = v2;
        aMatR1(1,k,i) = a1;
        aMatR1(2,k,i) = a2;
        TerMatR1(1,k,i) = Ter1;
        TerMatR1(2,k,i) = Ter2;
        clear v1 v2 a1 a2 Ter1 Ter2
        [v1,a1,Ter1] = ezdiffusion(HaccR2(1,k,i),HrtvR2(1,k,i),HrtR2(1,k,i),nTdrawn(k));  % second half
        [v2,a2,Ter2] = ezdiffusion(HaccR2(2,k,i),HrtvR2(2,k,i),HrtR2(2,k,i),nTdrawn(k));
        vMatR2(1,k,i) = v1;
        vMatR2(2,k,i) = v2;
        aMatR2(1,k,i) = a1;
        aMatR2(2,k,i) = a2;
        TerMatR2(1,k,i) = Ter1;
        TerMatR2(2,k,i) = Ter2;
        clear v1 v2 a1 a2 Ter1 Ter2

        % Bootstrap
        simrtC1 = [];
        simrtC2 = [];
        simrtI1 = [];
        simrtI2 = [];
        simrtvC1 = [];
        simrtvC2 = [];
        simrtvI1 = [];
        simrtvI2 = [];
        simaccC1 = [];
        simaccC2 = [];
        simaccI1 = [];
        simaccI2 = [];
        for i_b = 1:nBt
            tc1 = datasample(cMr1o,nTdrawn(k));
            tc2 = datasample(cMr2o,nTdrawn(k));
            ti1 = datasample(iMr1o,nTdrawn(k));
            ti2 = datasample(iMr2o,nTdrawn(k));
            simrtC1 = [simrtC1, mean(tc1)];
            simrtC2 = [simrtC2, mean(tc2)];
            simrtI1 = [simrtI1, mean(ti1)];
            simrtI2 = [simrtI2, mean(ti2)];
            simrtvC1 = [simrtvC1, var(tc1)];
            simrtvC2 = [simrtvC2, var(tc2)];
            simrtvI1 = [simrtvI1, var(ti1)];
            simrtvI2 = [simrtvI2, var(ti2)];
            simaccC1 = [simaccC1, mean(datasample(cMa1,nTdrawn(k)))];
            simaccC2 = [simaccC2, mean(datasample(cMa2,nTdrawn(k)))];
            simaccI1 = [simaccI1, mean(datasample(iMa1,nTdrawn(k)))];
            simaccI2 = [simaccI2, mean(datasample(iMa2,nTdrawn(k)))];
        end
        rt_bt1(1,k,:) = simrtC1;
        rt_bt2(1,k,:) = simrtC2;
        rt_bt1(2,k,:) = simrtI1;
        rt_bt2(2,k,:) = simrtI2;
        rtv_bt1(1,k,:) = simrtvC1;
        rtv_bt2(1,k,:) = simrtvC2;
        rtv_bt1(2,k,:) = simrtvI1;
        rtv_bt2(2,k,:) = simrtvI2;
        acc_bt1(1,k,:) = simaccC1;
        acc_bt2(1,k,:) = simaccC2;
        acc_bt1(2,k,:) = simaccI1;
        acc_bt2(2,k,:) = simaccI2;
        % EZ-diffusion modeling of the bootstrapped data
        for ib = 1:nBt
            [vSE1(1,k,ib),aSE1(1,k,ib),TerSE1(1,k,ib)] = ezdiffusion(acc_bt1(1,k,ib),rtv_bt1(1,k,ib),rt_bt1(1,k,ib),nTdrawn(k));
            [vSE2(1,k,ib),aSE2(1,k,ib),TerSE2(1,k,ib)] = ezdiffusion(acc_bt2(1,k,ib),rtv_bt2(1,k,ib),rt_bt2(1,k,ib),nTdrawn(k));
            [vSE1(2,k,ib),aSE1(2,k,ib),TerSE1(2,k,ib)] = ezdiffusion(acc_bt1(2,k,ib),rtv_bt1(2,k,ib),rt_bt1(2,k,ib),nTdrawn(k));
            [vSE2(2,k,ib),aSE2(2,k,ib),TerSE2(2,k,ib)] = ezdiffusion(acc_bt2(2,k,ib),rtv_bt2(2,k,ib),rt_bt2(2,k,ib),nTdrawn(k));
        end
        % Standard error
        VceSE1 = squeeze(vSE1(1,k,:)-vSE1(2,k,:));  % con-inc
        vICC1(k,:,i) = VceSE1;  % to later calculate 95% confidence interval of the ICC
        vMatRse(1,k,i) = (prctile(VceSE1,UpB)-prctile(VceSE1,LwB))/2;
        VceSE2 = squeeze(vSE2(1,k,:)-vSE2(2,k,:));
        vICC2(k,:,i) = VceSE2;
        vMatRse(2,k,i) = (prctile(VceSE2,UpB)-prctile(VceSE2,LwB))/2;
        AceSE1 = squeeze(aSE1(1,k,:)-aSE1(2,k,:));  % con-inc
        aICC1(k,:,i) = AceSE1;
        aMatRse(1,k,i) = (prctile(AceSE1,UpB)-prctile(AceSE1,LwB))/2;
        AceSE2 = squeeze(aSE2(1,k,:)-aSE2(2,k,:));
        aICC2(k,:,i) = AceSE2;
        aMatRse(2,k,i) = (prctile(AceSE2,UpB)-prctile(AceSE2,LwB))/2;
        TerceSE1 = squeeze(TerSE1(2,k,:)-TerSE1(1,k,:));  % inc-con
        TerICC1(k,:,i) = TerceSE1;
        TerMatRse(1,k,i) = (prctile(TerceSE1,UpB)-prctile(TerceSE1,LwB))/2;
        TerceSE2 = squeeze(TerSE2(2,k,:)-TerSE2(1,k,:));
        TerICC2(k,:,i) = TerceSE2;
        TerMatRse(2,k,i) = (prctile(TerceSE2,UpB)-prctile(TerceSE2,LwB))/2;
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

% x,y lim
xL = nan(drawLR,2,3);  % 2:min,max; 3:v,a,Ter
yL = nan(drawLR,2,3);
for k = 1:drawLR
    xL(k,1,1) = min(vTmp1(k,:)-squeeze(vMatRse(1,k,:))');
    xL(k,2,1) = max(vTmp1(k,:)+squeeze(vMatRse(1,k,:))');
    xL(k,1,2) = min(aTmp1(k,:)-squeeze(aMatRse(1,k,:))');
    xL(k,2,2) = max(aTmp1(k,:)+squeeze(aMatRse(1,k,:))');
    xL(k,1,3) = min(TerTmp1(k,:)-squeeze(TerMatRse(1,k,:))');
    xL(k,2,3) = max(TerTmp1(k,:)+squeeze(TerMatRse(1,k,:))');
    yL(k,1,1) = min(vTmp2(k,:)-squeeze(vMatRse(2,k,:))');
    yL(k,2,1) = max(vTmp2(k,:)+squeeze(vMatRse(2,k,:))');
    yL(k,1,2) = min(aTmp2(k,:)-squeeze(aMatRse(2,k,:))');
    yL(k,2,2) = max(aTmp2(k,:)+squeeze(aMatRse(2,k,:))');
    yL(k,1,3) = min(TerTmp2(k,:)-squeeze(TerMatRse(2,k,:))');
    yL(k,2,3) = max(TerTmp2(k,:)+squeeze(TerMatRse(2,k,:))');
end
vminX = min(xL(:,1,1));
vmaxX = max(xL(:,2,1));
vminY = min(yL(:,1,1));
vmaxY = max(yL(:,2,1));
aminX = min(xL(:,1,2));
amaxX = max(xL(:,2,2));
aminY = min(yL(:,1,2));
amaxY = max(yL(:,2,2));
TerminX = min(xL(:,1,3));
TermaxX = max(xL(:,2,3));
TerminY = min(yL(:,1,3));
TermaxY = max(yL(:,2,3));

% WITHOUT confidence intervals for the modeling parameters
for k = 1:drawLR
    % 1) RT CE
    figure
    scatter(ceTmp1(k,:),ceTmp2(k,:),60,[0 0.4 0.9],'filled')
    xlim([min([ceTmp1(:);ceTmp2(:)]) max([ceTmp1(:);ceTmp2(:)])])
    ylim([min([ceTmp1(:);ceTmp2(:)]) max([ceTmp1(:);ceTmp2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [Icc] = ICC([ceTmp1(k,:)',ceTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc);
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'FontSize',40,'FontWeight','bold','verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontSize',35,'FontWeight','bold')
    %title('Reaction time','FontSize',30)

    % 2) Acc CE
    figure
    scatter(ceTmpa1(k,:),ceTmpa2(k,:),60,[0 0.2 1],'filled')
    xlim([min([ceTmpa1(:);ceTmpa2(:)]) max([ceTmpa1(:);ceTmpa2(:)])])
    ylim([min([ceTmpa1(:);ceTmpa2(:)]) max([ceTmpa1(:);ceTmpa2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [Icc] = ICC([ceTmpa1(k,:)',ceTmpa2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'FontSize', 40,'FontWeight','bold','verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontSize',35,'FontWeight','bold')
    %title('Accuracy','FontSize',30)

    % 3) v
    figure
    scatter(vTmp1(k,:),vTmp2(k,:),60,cmap1(1,:),'filled')
    xlim([min([vTmp1(:);vTmp2(:)]) max([vTmp1(:);vTmp2(:)])])
    ylim([min([vTmp1(:);vTmp2(:)]) max([vTmp1(:);vTmp2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [Icc] = ICC([vTmp1(k,:)',vTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'FontSize',40,'FontWeight','bold','verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontSize',35,'FontWeight','bold')
    %title('Drift rate','FontSize',30)

    % 4) a
    figure
    scatter(aTmp1(k,:),aTmp2(k,:),60,cmap1(2,:),'filled')
    xlim([min([aTmp1(:);aTmp2(:)]) max([aTmp1(:);aTmp2(:)])])
    ylim([min([aTmp1(:);aTmp2(:)]) max([aTmp1(:);aTmp2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [Icc] = ICC([aTmp1(k,:)',aTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'FontSize',40,'FontWeight','bold','verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontSize',35,'FontWeight','bold')
    %title('Boundary separation','FontSize',30)

    % 5) Ter
    figure
    scatter(TerTmp1(k,:),TerTmp2(k,:),60,cmap1(3,:),'filled')
    xlim([min([TerTmp1(:);TerTmp2(:)]) max([TerTmp1(:);TerTmp2(:)])])
    ylim([min([TerTmp1(:);TerTmp2(:)]) max([TerTmp1(:);TerTmp2(:)])])
    ax = gca;
    ax.FontSize = 36;
    [Icc] = ICC([TerTmp1(k,:)',TerTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText,'FontSize',40,'FontWeight','bold','verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    %ylabel(['N = ' num2str(2*nTdrawn(k))],'FontSize',35,'FontWeight','bold')
    %title('Nondecision time','FontSize',30)
end

%% Plot ICC and its 95% confidence interval as a function of number of trials (Supp. Fig. 19B)
myICC = nan(drawLR,3);  % 3: EZ-diffusion model parameters
ICC_ci = nan(drawLR,3);
for k = 1:drawLR
    myICC(k,1) = ICC([vTmp1(k,:)',vTmp2(k,:)'],'A-k',0.05);
    myICC(k,2) = ICC([aTmp1(k,:)',aTmp2(k,:)'],'A-k',0.05);
    myICC(k,3) = ICC([TerTmp1(k,:)',TerTmp2(k,:)'],'A-k',0.05);
    tmp1 = nan(drawLR,nBt);  % drift rate
    tmp2 = nan(drawLR,nBt);  % boundary separation
    tmp3 = nan(drawLR,nBt);  % nondecision time
    for i_b = 1:nBt
        tmp1(k,i_b) = ICC([squeeze(vICC1(k,i_b,:)),squeeze(vICC2(k,i_b,:))],'A-k',0.05);
        tmp2(k,i_b) = ICC([squeeze(aICC1(k,i_b,:)),squeeze(aICC2(k,i_b,:))],'A-k',0.05);
        tmp3(k,i_b) = ICC([squeeze(TerICC1(k,i_b,:)),squeeze(TerICC2(k,i_b,:))],'A-k',0.05);
    end
    ICC_ci(k,1) = (prctile(tmp1(k,:),UpB)-prctile(tmp1(k,:),LwB))/2;
    ICC_ci(k,2) = (prctile(tmp2(k,:),UpB)-prctile(tmp2(k,:),LwB))/2;
    ICC_ci(k,3) = (prctile(tmp3(k,:),UpB)-prctile(tmp3(k,:),LwB))/2;
end
% 95% confidence interval
x = 1:drawLR;
xconf = [x x(end:-1:1)];
%figure
for pa = 1:3  % 3 EZD parameters
    %subplot(3,1,pa)
    figure
    CIl = myICC(:,pa)-ICC_ci(:,pa);
    CIu = myICC(:,pa)+ICC_ci(:,pa);
    u = CIu;
    yconf = [CIl u(end:-1:1)];
    pl = fill(xconf,yconf(:),'b');
    %pl.FaceColor = [1 0.5 0];
    pl.FaceColor = cmap1(pa,:);
    pl.EdgeColor = 'none';
    hold on
    %p2 = plot(x,myICC(:,pa),'Color',cmap1(pa,:),'LineWidth',2.4);
    p2 = plot(x,myICC(:,pa),'Color',[1 0.5 0],'LineWidth',2.4);
    set(gca,'FontSize',18)
    xlabel('Number of trials','FontSize',19)
    xticks(1:drawLR)
    xticklabels([50 100 200 400 800])
    legend('95% confidence interval','ICC','Location','southeast','FontSize',15)
    if pa == 1
        title('Drift rate','FontSize',20)
    elseif pa == 2
        title('Boundary separation','FontSize',20)
    elseif pa == 3
        title('Nondecision time','FontSize',20)
    end
    ylim([-0.2 1])
end