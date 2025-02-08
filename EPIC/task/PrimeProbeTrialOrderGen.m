numsubs=30;
subtrials=cell(61,numsubs);
for thesubs=1:numsubs;
 subID=['IndivVar00' num2str(thesubs)];
 for ivsess=1:60;
%Script Makes Trial Order for the ST-SI Blocks%
%1=II
%2=IC
%3=CC
%4=CI
%My cells%
ListIncAB=10;
ListIncYZ=20;
ListConAB=30;
ListConYZ=40;
NumBlock=8;
BrokenVector=cell(49,NumBlock);
for i=1:NumBlock;
    BaseVector=listgen1(4,3);
    for ii=1:length(BaseVector);
        if BaseVector(ii)==1;
            BrokenVector{ii,i}=[ListIncAB;ListIncYZ];
        elseif BaseVector(ii)==2;
            BrokenVector{ii,i}=[ListIncAB;ListConYZ];
        elseif BaseVector(ii)==3;
            BrokenVector{ii,i}=[ListConAB;ListConYZ];
        elseif BaseVector(ii)==4;
            BrokenVector{ii,i}=[ListConAB;ListIncYZ];
        else
            disp('problem')
        end 
    end
end
NumMat=[];
for w=1:size(BrokenVector,2);
    Colholder=[];
    for ww=1:size(BrokenVector,1);
        Unit=BrokenVector{ww,w};
        Colholder=[Colholder;Unit];
    end
    NumMat(:,w)=Colholder;
end
TrimNumMat=NumMat(3:size(NumMat,1),:);
CheckConA=sum(TrimNumMat(:,1) == 30);
CheckConB=sum(TrimNumMat(:,1) == 40);
CheckIncA=sum(TrimNumMat(:,1) == 10);
CheckIncB=sum(TrimNumMat(:,1) == 20);
%Test II IC CC CI Balance: CC=1, CI=2, II=3, IC=4%
BalanceCheck=[];
for q=2:length(TrimNumMat(:,1));
    if TrimNumMat(q,1)==30 & (TrimNumMat(q-1,1)==30|TrimNumMat(q-1,1)==40);
        BalanceCheck(q)=1;
    elseif TrimNumMat(q,1)==40 & (TrimNumMat(q-1,1)==30|TrimNumMat(q-1,1)==40); %don't the or statement can remove the on that correpsonds to B or vice versa in the case of an a trial%
        BalanceCheck(q)=1;
    elseif TrimNumMat(q,1)==10 & (TrimNumMat(q-1,1)==30|TrimNumMat(q-1,1)==40);
        BalanceCheck(q)=2;
    elseif TrimNumMat(q,1)==20 & (TrimNumMat(q-1,1)==30|TrimNumMat(q-1,1)==40);
        BalanceCheck(q)=2;
    elseif TrimNumMat(q,1)==10 & (TrimNumMat(q-1,1)==10|TrimNumMat(q-1,1)==20);
        BalanceCheck(q)=3;
    elseif TrimNumMat(q,1)==20 & (TrimNumMat(q-1,1)==10|TrimNumMat(q-1,1)==20);
        BalanceCheck(q)=3;
    elseif TrimNumMat(q,1)==30 & (TrimNumMat(q-1,1)==10|TrimNumMat(q-1,1)==20);
        BalanceCheck(q)=4;
    elseif TrimNumMat(q,1)==40 & (TrimNumMat(q-1,1)==10|TrimNumMat(q-1,1)==20);
        BalanceCheck(q)=4;
    else
        disp('problem')
    end
end
CheckCC=sum(BalanceCheck(:) == 1);
CheckCI=sum(BalanceCheck(:) == 2);
CheckII=sum(BalanceCheck(:) == 3);
CheckIC=sum(BalanceCheck(:) == 4);
Eprimestack=[];
for t=1:size(TrimNumMat,2);
    grab=TrimNumMat(:,t);
    Eprimestack=[Eprimestack;grab];
end
global EprimeCell;
EprimeCell=cell(length(Eprimestack),1);
for qu=1:length(Eprimestack);
    if Eprimestack(qu)==10;
        EprimeCell{qu,1}='ListIncAB';
    elseif Eprimestack(qu)==20;
        EprimeCell{qu,1}='ListIncYZ';
    elseif Eprimestack(qu)==30;
        EprimeCell{qu,1}='ListConAB';
    elseif Eprimestack(qu)==40;
        EprimeCell{qu,1}='ListConYZ';
    else
        disp('problem')
    end
end
%BreakUp into Two Segs%
STandSIOne=cell(384,1);
STandSITwo=cell(384,1);
LengthofSeg=(length(EprimeCell)/2);
for BOne=1:LengthofSeg
STandSIOne{BOne,1}=EprimeCell{BOne,1};
end
for BTwo=1:LengthofSeg
STandSITwo{BTwo,1}=EprimeCell{(BTwo+384),1};
end
EprimestackUsed=Eprimestack(1:384);
stimorder=[];
Eprimeblocks=[EprimestackUsed(1:96),EprimestackUsed(97:192),EprimestackUsed(193:288),EprimestackUsed(289:384)];
for theblocks=1:4;
EprimeTarget=Eprimeblocks(:,theblocks);
YZcon=Shuffle([[ones(1,12)],[ones(1,12)*3]]);
YZinc=Shuffle([[ones(1,12)*5],[ones(1,12)*7]]);
ABcon=Shuffle([[ones(1,12)*2],[ones(1,12)*4]]);
ABinc=Shuffle([[ones(1,12)*6],[ones(1,12)*8]]);
stimblock=[];
for eachstim=1:length(EprimeTarget);
    %10 & 30 are the AB; 20 & 40 are the YZ; 10 & 20 are Inc; 30 & 40 are Con%
    %2=AA 4=BB 6=AB 8=BA; 1=YY 3=ZZ 5=YZ 7=ZY%
    if EprimeTarget(eachstim)==10
        stimblock=[stimblock;ABinc(1)];
        ABinc(1)=[];
    elseif EprimeTarget(eachstim)==20
        stimblock=[stimblock;YZinc(1)];
        YZinc(1)=[];
    elseif EprimeTarget(eachstim)==30
        stimblock=[stimblock;ABcon(1)];
        ABcon(1)=[];
    elseif EprimeTarget(eachstim)==40
        stimblock=[stimblock;YZcon(1)];
        YZcon(1)=[];
    end
end
stimorder=[stimorder;stimblock];
end
subtrials{ivsess,thesubs}=stimorder;
 end
subtrials{61,thesubs}=subID;
end