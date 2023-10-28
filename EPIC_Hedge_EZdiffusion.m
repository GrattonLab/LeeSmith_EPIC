%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
%
% To run this sciprt:
% You need Hedge et al.'s (2018) flanker and Stroop task data 
% (download at https://osf.io/cwzds/),
% a Matlab function named "ezdiffusion.m" (download at https://osf.io/jk9nb/),
% and ICC.m (download at
% https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc)
%
% What this script does:
% Plots the EZ-diffusion model parameters across varying numbers of trials
% (RT needs to be in SECS to run EZ-diffusion modeling)
%
% What this script outputs:
% Supp. Fig. 12B
%
% Created on 05/16/2023 by HJ Lee
% Last modified on 07/03/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

rng('shuffle')

%% Parameter settings
nCond = 2;  % number of experimental conditions: congruent/incongruent 
nStd = 3;  % STD from mean for excluding outliers

tmpSubjID1 = 1:50;
tmpSubjID2 = 1:62;
subjExcld1 = [6, 17, 37, 21, 33];  % based on Hedge et al.'s Readme.txt; also removed 21 33
subjExcld2 = [25, 28, 34, 38, 56, 32];  % based on Readme.txt; also removed 32
subjID1 = setdiff(tmpSubjID1, subjExcld1);
subjID2 = setdiff(tmpSubjID2, subjExcld2);
nSubj1 = length(subjID1);  % study 1
nSubj2 = length(subjID2);  % study 2
nSubjT = nSubj1+nSubj2;
nTask = 1;  % just use the flanker task data for EZdiffusion; otherwise 2
taskStrng = {'Flanker','Stroop'};
%nSession = 2;
nBlocks = 5;
nTrials = 144;  % per block
%nTrialsCE = 96;  % number of trials for calculating CE (excluding neutral)

nTdrawn = [50, 100, 200, 400, 800, 1600 3200]/2;  % number of trials drawn (per condition)
drawLH = 4;  % number of trials drawn for Hedge et al.'s data; max 480 trials (cf. nTdrawn)

%% 1. Load and preprocess Hedge et al.'s (2018) flanker task data
HrtMat1 = nan(nCond,length(nTdrawn),nSubjT);  % mean RT (first half) in secs
HrtMat2 = nan(nCond,length(nTdrawn),nSubjT);  % second half
HrtvMat1 = nan(nCond,length(nTdrawn),nSubjT);  % RT variance
HrtvMat2 = nan(nCond,length(nTdrawn),nSubjT);
HaccMat1 = nan(nCond,length(nTdrawn),nSubjT);  % mean accuracy
HaccMat2 = nan(nCond,length(nTdrawn),nSubjT);
% Study 1 data
for i = 1:nSubj1
    for t = 1:nTask  % 1
        dir = ['HedgeData/Study1-' taskStrng{t}];
        cd(dir)
        % session 1
        fileID  = ['Study1_P' num2str(subjID1(i)) taskStrng{t} num2str(1) '.csv'];
        rawData = readmatrix(fileID);
        T1 = table((1:nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T1(T1.n0cong==1,:) = [];  % exclude neutral condition
        % session 2
        fileID  = ['Study1_P' num2str(subjID1(i)) taskStrng{t} num2str(2) '.csv'];
        rawData = readmatrix(fileID);
        T2 = table((1+nBlocks*nTrials:2*nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T2(T2.n0cong==1,:) = [];  % exclude neutral condition
        %T = [T1;T2];  % combine session 1&2
        cd('../../')  % home directory
        % Exclude outlier trials
        % T1
        T1_rm = T1(T1.acc==1,:);
        TGM1 = varfun(@mean,T1_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS1 = varfun(@std,T1_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        inbound = TGM1.mean_rt-(TGS1.std_rt*nStd);  % 3SD
        for i_inbound = 1:size(inbound,1)
            if inbound(i_inbound) < 0
                inbound(i_inbound) = 0;
            end
        end
        %inbound = 0.15;  % alternatively
        outbound = TGM1.mean_rt+(TGS1.std_rt*nStd);
        out_ids = [];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==0 & T1_rm.rt<inbound(1))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==2 & T1_rm.rt<inbound(2))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==0 & T1_rm.rt>outbound(1))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==2 & T1_rm.rt>outbound(2))];
        T1_rm.acc(ismember(T1_rm.trialNum,out_ids)) = 99;
        T1.acc(ismember(T1.trialNum,out_ids)) = 99;
        % T2
        T2_rm = T2(T2.acc==1,:);
        TGM2 = varfun(@mean,T2_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS2 = varfun(@std,T2_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        inbound = TGM2.mean_rt-(TGS2.std_rt*nStd);
        for i_inbound = 1:size(inbound,1)
            if inbound(i_inbound) < 0
                inbound(i_inbound) = 0;
            end
        end
        %inbound = 0.15;  % alternatively
        outbound = TGM2.mean_rt+(TGS2.std_rt*nStd);  % upper bound
        out_ids = [];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==0 & T2_rm.rt<inbound(1))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==2 & T2_rm.rt<inbound(2))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==0 & T2_rm.rt>outbound(1))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==2 & T2_rm.rt>outbound(2))];
        T2_rm.acc(ismember(T2_rm.trialNum,out_ids)) = 99;
        T2.acc(ismember(T2.trialNum,out_ids)) = 99;

        % descriptive
        rtC1 = T1.rt(and(T1.n0cong==0,T1.acc==1));  % con
        rtI1 = T1.rt(and(T1.n0cong==2,T1.acc==1));  % inc
        rtC2 = T2.rt(and(T2.n0cong==0,T2.acc==1));
        rtI2 = T2.rt(and(T2.n0cong==2,T2.acc==1));
        accC1 = T1.acc(T1.n0cong==0);
        accI1 = T1.acc(T1.n0cong==2);
        accC2 = T2.acc(T2.n0cong==0);
        accI2 = T2.acc(T2.n0cong==2);
        accC1(find(accC1==99)) = [];  % removing outlier trials
        accI1(find(accI1==99)) = [];
        accC2(find(accC2==99)) = [];
        accI2(find(accI2==99)) = [];
        for k = 1:drawLH
            if length(rtC1) < nTdrawn(k)  % this happens when k == 4
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtC1))])
                HrtMat1(1,k,i) = mean(rtC1);  % don't need trimmean because outlier trials were marked as 99 and removed
                HrtvMat1(1,k,i) = var(rtC1);
                HaccMat1(1,k,i) = mean(accC1);
            else
                HrtMat1(1,k,i) = mean(rtC1(1:nTdrawn(k)));
                HrtvMat1(1,k,i) = var(rtC1(1:nTdrawn(k)));
                HaccMat1(1,k,i) = mean(accC1(1:nTdrawn(k)));
            end
            if length(rtI1) < nTdrawn(k)
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtI1))])
                HrtMat1(2,k,i) = mean(rtI1); 
                HrtvMat1(2,k,i) = var(rtI1);
                HaccMat1(2,k,i) = mean(accI1);
            else
                HrtMat1(2,k,i) = mean(rtI1(1:nTdrawn(k)));
                HrtvMat1(2,k,i) = var(rtI1(1:nTdrawn(k)));
                HaccMat1(2,k,i) = mean(accI1(1:nTdrawn(k)));
            end
            if length(rtC2) < nTdrawn(k)
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtC2))])
                HrtMat2(1,k,i) = mean(rtC2);
                HrtvMat2(1,k,i) = var(rtC2);
                HaccMat2(1,k,i) = mean(accC2);
            else
                HrtMat2(1,k,i) = mean(rtC2(1:nTdrawn(k)));
                HrtvMat2(1,k,i) = var(rtC2(1:nTdrawn(k)));
                HaccMat2(1,k,i) = mean(accC2(1:nTdrawn(k)));
            end
            if length(rtI2) < nTdrawn(k)
                disp(['Subj' num2str(i) ' step' num2str(k) ' #ofT:' num2str(length(rtI2))])
                HrtMat2(2,k,i) = mean(rtI2);
                HrtvMat2(2,k,i) = var(rtI2);
                HaccMat2(2,k,i) = mean(accI2);
            else
                HrtMat2(2,k,i) = mean(rtI2(1:nTdrawn(k)));
                HrtvMat2(2,k,i) = var(rtI2(1:nTdrawn(k)));
                HaccMat2(2,k,i) = mean(accI2(1:nTdrawn(k)));
            end
        end
    end
end
% Study 2 data
for i = 1:nSubj2
    for t = 1:nTask
        dir = ['HedgeData/Study2-' taskStrng{t}];
        cd(dir)
        % session 1
        fileID  = ['Study2_P' num2str(subjID2(i)) taskStrng{t} num2str(1) '.csv'];
        rawData = readmatrix(fileID);
        T1 = table((1:nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T1(T1.n0cong==1,:) = [];  % exclude neutral condition
        % session 2
        fileID  = ['Study2_P' num2str(subjID2(i)) taskStrng{t} num2str(2) '.csv'];
        rawData = readmatrix(fileID);
        T2 = table((1+nBlocks*nTrials:2*nBlocks*nTrials)',rawData(:,4),rawData(:,5),rawData(:,6),'VariableNames',["trialNum","n0cong","acc","rt"]);
        T2(T2.n0cong==1,:) = [];  % exclude neutral condition
        %T = [T1;T2];
        cd('../../')  % home directory
        % Exclude outlier trials
        % T1
        T1_rm = T1(T1.acc==1,:);
        TGM1 = varfun(@mean,T1_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS1 = varfun(@std,T1_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        inbound = TGM1.mean_rt-(TGS1.std_rt*nStd);
        for i_inbound = 1:size(inbound,1)
            if inbound(i_inbound) < 0
                inbound(i_inbound) = 0;
            end
        end
        %inbound = 0.15;  % alternatively
        outbound = TGM1.mean_rt+(TGS1.std_rt*nStd);  % upper bound
        out_ids = [];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==0 & T1_rm.rt<inbound(1))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==2 & T1_rm.rt<inbound(2))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==0 & T1_rm.rt>outbound(1))];
        out_ids = [out_ids; T1_rm.trialNum(T1_rm.n0cong==2 & T1_rm.rt>outbound(2))];
        T1_rm.acc(ismember(T1_rm.trialNum,out_ids)) = 99;
        T1.acc(ismember(T1.trialNum,out_ids)) = 99;
        % T2
        T2_rm = T2(T2.acc==1,:);
        TGM2 = varfun(@mean,T2_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        TGS2 = varfun(@std,T2_rm,'GroupingVariables',{'n0cong'},'OutputFormat','table');
        inbound = TGM2.mean_rt-(TGS2.std_rt*nStd);
        for i_inbound = 1:size(inbound,1)
            if inbound(i_inbound) < 0
                inbound(i_inbound) = 0;
            end
        end
        %inbound = 0.15;  % alternatively
        outbound = TGM2.mean_rt+(TGS2.std_rt*nStd);  % upper bound
        out_ids = [];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==0 & T2_rm.rt<inbound(1))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==2 & T2_rm.rt<inbound(2))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==0 & T2_rm.rt>outbound(1))];
        out_ids = [out_ids; T2_rm.trialNum(T2_rm.n0cong==2 & T2_rm.rt>outbound(2))];
        T2_rm.acc(ismember(T2_rm.trialNum,out_ids)) = 99;
        T2.acc(ismember(T2.trialNum,out_ids)) = 99;

        % descriptive (cont.)
        rtC1 = T1.rt(and(T1.n0cong==0,T1.acc==1));
        rtI1 = T1.rt(and(T1.n0cong==2,T1.acc==1));
        rtC2 = T2.rt(and(T2.n0cong==0,T2.acc==1));
        rtI2 = T2.rt(and(T2.n0cong==2,T2.acc==1));
        accC1 = T1.acc(T1.n0cong==0);
        accI1 = T1.acc(T1.n0cong==2);
        accC2 = T2.acc(T2.n0cong==0);
        accI2 = T2.acc(T2.n0cong==2);
        accC1(find(accC1==99)) = [];
        accI1(find(accI1==99)) = [];
        accC2(find(accC2==99)) = [];
        accI2(find(accI2==99)) = [];
        for k = 1:drawLH
            if length(rtC1) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtC1))])
                HrtMat1(1,k,nSubj1+i) = mean(rtC1);  % don't need trimmean because outlier trials are marked as 99 above and removed
                HrtvMat1(1,k,nSubj1+i) = var(rtC1);
                HaccMat1(1,k,nSubj1+i) = mean(accC1);
            else
                HrtMat1(1,k,nSubj1+i) = mean(rtC1(1:nTdrawn(k)));
                HrtvMat1(1,k,nSubj1+i) = var(rtC1(1:nTdrawn(k)));
                HaccMat1(1,k,nSubj1+i) = mean(accC1(1:nTdrawn(k)));
            end
            if length(rtI1) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtI1))])
                HrtMat1(2,k,nSubj1+i) = mean(rtI1);
                HrtvMat1(2,k,nSubj1+i) = var(rtI1);
                HaccMat1(2,k,nSubj1+i) = mean(accI1);
            else
                HrtMat1(2,k,nSubj1+i) = mean(rtI1(1:nTdrawn(k)));
                HrtvMat1(2,k,nSubj1+i) = var(rtI1(1:nTdrawn(k)));
                HaccMat1(2,k,nSubj1+i) = mean(accI1(1:nTdrawn(k)));
            end
            if length(rtC2) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtC2))])
                HrtMat2(1,k,nSubj1+i) = mean(rtC2);
                HrtvMat2(1,k,nSubj1+i) = var(rtC2);
                HaccMat2(1,k,nSubj1+i) = mean(accC2);
            else
                HrtMat2(1,k,nSubj1+i) = mean(rtC2(1:nTdrawn(k)));
                HrtvMat2(1,k,nSubj1+i) = var(rtC2(1:nTdrawn(k)));
                HaccMat2(1,k,nSubj1+i) = mean(accC2(1:nTdrawn(k)));
            end
            if length(rtI2) < nTdrawn(k)
                disp(['Subj' num2str(nSubj1+i) ' step' num2str(k) ' #ofT:' num2str(length(rtI2))])
                HrtMat2(2,k,nSubj1+i) = mean(rtI2);
                HrtvMat2(2,k,nSubj1+i) = var(rtI2);
                HaccMat2(2,k,nSubj1+i) = mean(accI2);
            else
                HrtMat2(2,k,nSubj1+i) = mean(rtI2(1:nTdrawn(k)));
                HrtvMat2(2,k,nSubj1+i) = var(rtI2(1:nTdrawn(k)));
                HaccMat2(2,k,nSubj1+i) = mean(accI2(1:nTdrawn(k)));
            end
        end
    end
end

%% EZ-diffusion modeling using Hedge et al.'s data - To compare it with simulated data for evaulating the model
% run EZ-diffusion modeling across different number of trial draws
vMatH1 = nan(nCond,drawLH,nSubjT);
aMatH1 = nan(nCond,drawLH,nSubjT);
TerMatH1 = nan(nCond,drawLH,nSubjT);
vMatH2 = nan(nCond,drawLH,nSubjT);
aMatH2 = nan(nCond,drawLH,nSubjT);
TerMatH2 = nan(nCond,drawLH,nSubjT);
for i = 1:nSubjT
    for k = 1:drawLH
        [v1, a1, Ter1] = ezdiffusion(HaccMat1(1,k,i),HrtvMat1(1,k,i),HrtMat1(1,k,i),nTdrawn(k));  % con
        [v2, a2, Ter2] = ezdiffusion(HaccMat1(2,k,i),HrtvMat1(2,k,i),HrtMat1(2,k,i),nTdrawn(k));  % inc
        vMatH1(1,k,i) = v1;  % con
        vMatH1(2,k,i) = v2;  % inc
        aMatH1(1,k,i) = a1;
        aMatH1(2,k,i) = a2;
        TerMatH1(1,k,i) = Ter1;
        TerMatH1(2,k,i) = Ter2;
        clear v1 v2 a1 a2 Ter1 Ter2
        [v1, a1, Ter1] = ezdiffusion(HaccMat2(1,k,i),HrtvMat2(1,k,i),HrtMat2(1,k,i),nTdrawn(k));
        [v2, a2, Ter2] = ezdiffusion(HaccMat2(2,k,i),HrtvMat2(2,k,i),HrtMat2(2,k,i),nTdrawn(k));
        vMatH2(1,k,i) = v1;
        vMatH2(2,k,i) = v2;
        aMatH2(1,k,i) = a1;
        aMatH2(2,k,i) = a2;
        TerMatH2(1,k,i) = Ter1;
        TerMatH2(2,k,i) = Ter2;
        clear v1 v2 a1 a2 Ter1 Ter2
    end
end
% Plot EZD modeling results
cmap1 = parula(3);
cmap2 = spring(3);
cmap3 = cool(3);
ceTmp1 = squeeze(HrtMat1(2,:,:)-HrtMat1(1,:,:));  % inc-con
ceTmp2 = squeeze(HrtMat2(2,:,:)-HrtMat2(1,:,:));
ceTmpa1 = squeeze(HaccMat1(1,:,:)-HaccMat1(2,:,:));  % con-inc
ceTmpa2 = squeeze(HaccMat2(1,:,:)-HaccMat2(2,:,:));
vTmp1 = squeeze(vMatH1(1,:,:)-vMatH1(2,:,:));  % con-inc
vTmp2 = squeeze(vMatH2(1,:,:)-vMatH2(2,:,:));
aTmp1 = squeeze(aMatH1(1,:,:)-aMatH1(2,:,:));
aTmp2 = squeeze(aMatH2(1,:,:)-aMatH2(2,:,:));
TerTmp1  = squeeze(TerMatH1(2,:,:)-TerMatH1(1,:,:));  % inc-con
TerTmp2 = squeeze(TerMatH2(2,:,:)-TerMatH2(1,:,:));
figure
for k = 1:drawLH
    % 1) RT CE
    subplot(drawLH+1,5,5*(k-1)+1)
    scatter(ceTmp1(k,:),ceTmp2(k,:),8,'filled')
    %lsline
    xlim([min([ceTmp1(:);ceTmp2(:)]) max([ceTmp1(:);ceTmp2(:)])])
    ylim([min([ceTmp1(:);ceTmp2(:)]) max([ceTmp1(:);ceTmp2(:)])])
    [Icc] = ICC([ceTmp1(k,:)',ceTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc);
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
    ylabel(['N = ' num2str(2*nTdrawn(k))])

    % 2) Acc CE
    subplot(drawLH+1,5,5*(k-1)+2)
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
    subplot(drawLH+1,5,5*(k-1)+3)
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
    subplot(drawLH+1,5,5*(k-1)+4)
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
    subplot(drawLH+1,5,5*(k-1)+5)
    scatter(TerTmp1(k,:),TerTmp2(k,:),8,cmap1(3,:),'filled')
    xlim([min([TerTmp1(:);TerTmp2(:)]) max([TerTmp1(:);TerTmp2(:)])])
    ylim([min([TerTmp1(:);TerTmp2(:)]) max([TerTmp1(:);TerTmp2(:)])])
    %lsline
    [Icc] = ICC([TerTmp1(k,:)',TerTmp2(k,:)'],'A-k',0.05);
    str = sprintf('  ICC = %1.2f', Icc); clear Icc
    rText = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(rText, 'fontsize', 8, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    refline
    axis square
end
subplot(drawLH+1,5,1); title('CE RT')
subplot(drawLH+1,5,2); title('CE accuracy')
subplot(drawLH+1,5,3); title('Drift rate')
subplot(drawLH+1,5,4); title('Boundary separation')
subplot(drawLH+1,5,5); title('Nondecision time')
sgtitle('Hedge et al. (2018)')