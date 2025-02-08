function tt_rc_StroopTask
%By DMS based on template code from MD
sca;
close all;
clearvars;

rng('shuffle');

%--------------------------------------------------------------------------
%                        Participant Setup
%--------------------------------------------------------------------------

prompt = 'Participant #: ';
subID = input(prompt, 's');

current = pwd();
subFolder = [current '/ParticipantInfo/' subID];
%Create Participant folder (if new)
if ~isfolder(subFolder)
    mkdir(subFolder)
end
taskFolder = [subFolder '/STROOP'];

if ~isfolder(taskFolder)
    mkdir(taskFolder)
end


overwrite = 0;
while overwrite == 0
    prompt = 'Run #: ';
    run = input(prompt);
    cd(taskFolder)
    if isfile(['run' num2str(run) '.mat'])
        prompt = 'Run already exists. Overwirte? (y/n): ';
        ow = input(prompt, 's');
        if strcmpi(ow, 'y')
            overwrite = 1;
        end
    else
        break
    end
end

cd(current)

prompt = 'Practice? (y/n): ';
practice = input(prompt, 's');

if practice == 'y'
    practice = 1;
else
    practice = 0;
end
%--------------------------------------------------------------------------
%                       Change Variables here
%--------------------------------------------------------------------------

blocks = 4; 
trialsperblock = 108;
trials = blocks*trialsperblock; 
trialBreakdown = [(1/3), (1/3), (1/3)]; 
stimTime=1;
ITI = 1; %Secs between words
timeLimit = 2; %Time limit to answer 

numColors = 4;

allRGB = zeros(numColors, 3);
allRGB(1, :) = [1 0 0]; %Red
allRGB(2, :) = [1 1 0]; %Yellow
allRGB(3, :) = [0 1 0]; %Green
allRGB(4, :) = [0 0 1]; %Blue


allNames = {'RED', 'YELLOW', 'GREEN', 'BLUE'};

otherWords = {'HORSE', 'BIRD', 'CAT','DOG'}; 


oneKey = KbName('1!');
twoKey = KbName('2@');
threeKey = KbName('3#');
fourKey = KbName('4$');


allData = cell(trials, 8);
oldBlockRT = 2;
block = 1;

%--------------------------------------------------------------------------
%                              Screen Setup
%--------------------------------------------------------------------------

% Skip sync tests to avoid error
Screen('Preference', 'SkipSyncTests', 1);

PsychDefaultSetup(2);

screens = Screen('Screens'); 
screenNumber = max(screens);

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% Center of Screen
[xCenter, yCenter] = RectCenter(windowRect);

HideCursor()

Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%Screen Text 
Screen('TextFont', window, 'Helvetica');
Screen('TextSize', window, 36); %Font size for instructions 

%Fixation cross 
fixCrossDimPix = 15;
lineWidthPix = 2;
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];



%--------------------------------------------------------------------------
%                          Word Setup
%--------------------------------------------------------------------------
asterisks = zeros(trialBreakdown(1)*trials, 1);
differents = ones(trialBreakdown(2)*trials, 1);
sames = ones(trialBreakdown(3)*trials, 1) *2;

trialOrder = Shuffle(vertcat(asterisks, differents, sames));

%No condition more than 3 in a row 
fixed = 0;
while ~fixed
    for i = 4:trials
        if trialOrder(i) == trialOrder(i-1) && trialOrder(i) == trialOrder(i-2) && trialOrder(i) == trialOrder(i-3)
            trialOrder = Shuffle(trialOrder);
            break
        elseif i == trials
            fixed = 1;
            break
        end
    end
end

% ---For only 4 colors ------
% colorOrder = datasample(1:numColors, trials);
wordOrderCon=[];
colorOrderCon=[];
wordOrderInc=[];
colorOrderInc=[];
Neutral=[];
for bl=1:blocks;
thewordc=[];
thecolorc=[];
thewordi=[];
thecolori=[];
convec=Shuffle(repmat([1,2,3,4],1,9));
incvec=Shuffle(repmat([1,2,3,4,5,6,7,8,9,10,11,12],1,3));
nevec=Shuffle(repmat([1,2,3,4],1,9));
%Mappings incvec to colorOrder 1,2,3 are 1; 4,5,6 are 2; 7,8,9 are 3; 10,11,12 are 4;
%2,3,4 ; 1,3,4 ; 1,2,4; 1,2,3;
%Mappings incvec to wordOrder 4,7,10 are 1; 1,8,11 are 2; 2,5,12 are 3; 3,6,9;
for tr=1:36;
    if convec(tr)==1;
        thewordc=[thewordc,1];
        thecolorc=[thecolorc,1];
    elseif convec(tr)==2;
        thewordc=[thewordc,2];
        thecolorc=[thecolorc,2];
    elseif convec(tr)==3;
        thewordc=[thewordc,3];
        thecolorc=[thecolorc,3];
    elseif convec(tr)==4;
        thewordc=[thewordc,4];
        thecolorc=[thecolorc,4];
    end
    if incvec(tr)==4 || incvec(tr)==7 || incvec(tr)==10;
        thewordi=[thewordi,1];
    elseif incvec(tr)==1 || incvec(tr)==8 || incvec(tr)==11;
        thewordi=[thewordi,2];
    elseif incvec(tr)==2 || incvec(tr)==5 || incvec(tr)==12;
        thewordi=[thewordi,3];
    elseif incvec(tr)==3 || incvec(tr)==6 || incvec(tr)==9;
        thewordi=[thewordi,4];
    end
    if incvec(tr)==1 || incvec(tr)==2 || incvec(tr)==3;
        thecolori=[thecolori,1];
    elseif incvec(tr)==4 || incvec(tr)==5 || incvec(tr)==6;
        thecolori=[thecolori,2];
    elseif incvec(tr)==7 || incvec(tr)==8 || incvec(tr)==9;
        thecolori=[thecolori,3];
    elseif incvec(tr)==10 || incvec(tr)==11 || incvec(tr)==12;
        thecolori=[thecolori,4];
    end
end
wordOrderCon=[wordOrderCon,thewordc];
colorOrderCon=[colorOrderCon,thecolorc];
wordOrderInc=[wordOrderInc,thewordi];
colorOrderInc=[colorOrderInc,thecolori];
Neutral=[Neutral,nevec];
end

NextC=1;
NextI=1;
NextN=1;
wordOrder = NaN(trials, 1);
colorOrder = [];
for word=1:trials;
    if trialOrder(word) == 0 %asterisks/other word
        wordOrder(word) = 0;
        colorOrder(word)=Neutral(NextN);
        NextN=NextN+1;
    elseif trialOrder(word) == 1 %Different
        wordOrder(word)=wordOrderInc(NextI);
        colorOrder(word)=colorOrderInc(NextI);
        NextI=NextI+1;
    else %Same
        wordOrder(word)=wordOrderCon(NextC);
        colorOrder(word)=colorOrderCon(NextC);
        NextC=NextC+1;
    end
end

%--------------------------------------------------------------------------       
%                        Instructions/Key Bindings 
%--------------------------------------------------------------------------

DrawFormattedText(window, 'You will see a word appear in a color that may or may not match the word.', 'center', yCenter-300, white);
DrawFormattedText(window, 'You should respond with the COLOR the word is written in, NOT the word itself.', 'center', yCenter-200, white);
DrawFormattedText(window, 'Respond using the number of the color.', 'center', yCenter-100, white);
DrawFormattedText(window, 'You should work as quickly and as accurately as possible.', 'center', yCenter, white);

DrawFormattedText(window, '1=RED      2=YELLOW      3=GREEN      4=BLUE', 'center', yCenter+100, white);
DrawFormattedText(window, 'Note: Use Right Hand, 1=Index      2=Middle      3=Ring      4=Pinky', 'center', yCenter+200, white);
DrawFormattedText(window, 'Press ANY KEY to start', 'center', yCenter+300, white);

Screen('Flip', window)
WaitSecs(.5);
KbWait;

Screen('DrawLines', window, allCoords,...
    lineWidthPix, white, [xCenter, yCenter], 2);
Screen('TextSize', window, 28);
DrawFormattedText(window, '1=RED      2=YELLOW      3=GREEN      4=BLUE', 'center', yCenter+400, white);
Screen('Flip', window);

WaitSecs(2); %Wait for time 

if practice == 1 %check for practice 
    
%--------------------------------------------------------------------------
%                          Practice
%--------------------------------------------------------------------------
pracTrialOrder = [0 0 0 0 0 0 2 2 0 2 0 2 0 0 2 1 1 1 2 0 2 1 0 1 1];
numPracTrials = 25;
pracColors = randi(numColors, numPracTrials, 1);
pracWords = NaN(numPracTrials, 1);

for word = 1:numPracTrials
    if pracTrialOrder(word) == 0 %Other word
        pracWords(word) = 0;
    elseif pracTrialOrder(word) == 1 % different 
        pracWords(word) = randi(numColors);
        while pracWords(word) == pracColors(word)
            pracWords(word) = randi(numColors);
        end
    else %Same 
        pracWords(word) = randi(numColors);
    end
end


for pracTrial = 1:numPracTrials
    trialColor = allRGB(pracColors(pracTrial), :);
    Screen('TextSize', window, 48);
        trialCenter = yCenter+24;
    if pracTrialOrder(pracTrial) == 0
        trialWord = cell2mat(datasample(otherWords, 1));

    else
        trialWord = allNames{pracWords(pracTrial)};
    end
    
    DrawFormattedText(window, trialWord, 'center', trialCenter, trialColor);
    Screen('TextSize', window, 28)
    DrawFormattedText(window, '1=RED      2=YELLOW      3=GREEN      4=BLUE', 'center', yCenter+400, white);
    Screen('Flip', window);
    
    while 1
       [~ ,~, keyCode] = KbCheck;
        if keyCode(oneKey)
            response = 1; %Red
            break
        elseif keyCode(twoKey)
            response = 2; %Yellow
            break
        elseif keyCode(threeKey)
            response = 3; %Green
            break
        elseif keyCode(fourKey)
            response = 4; %Blue
            break
        end
    end
    
    
    %Accuracy 
    acc = response == pracColors(pracTrial);
    
    
    if acc == 1
        message = ['Correct! The word is printed in ' allNames{pracColors(pracTrial)}];
    else
        message = ['Incorrect. The word is printed in ' allNames{pracColors(pracTrial)}];
    end
    
    
    Screen('TextSize', window, 36);
    DrawFormattedText(window, message, 'center', yCenter-100, white);
    DrawFormattedText(window, trialWord, 'center', trialCenter, trialColor);
    DrawFormattedText(window, 'Press ANY KEY to continue', 'center', yCenter+300, white);
    Screen('TextSize', window, 28)
    DrawFormattedText(window, '1=RED      2=YELLOW      3=GREEN      4=BLUE', 'center', yCenter+400, white);
    Screen('Flip', window);
    WaitSecs(.5);
    KbWait;
    DrawFormattedText(window, '1=RED      2=YELLOW      3=GREEN      4=BLUE', 'center', yCenter+400, white);
    Screen('DrawLines', window, allCoords,...
        lineWidthPix, white, [xCenter, yCenter], 2);
    Screen('Flip', window);
    WaitSecs(ITI);
    
end

Screen('TextSize', window, 36);
DrawFormattedText(window, 'Practice complete!', 'center', yCenter-100, white)
DrawFormattedText(window, 'During the actual task, you will not receive feedback until the end of a block.',...
                        'center', 'center', white)
DrawFormattedText(window, 'Press ANY KEY to begin actual task', 'center', yCenter+300, white);
Screen('Flip', window)
WaitSecs(.5);
KbWait
Screen('Flip', window);
WaitSecs(2);


end
%--------------------------------------------------------------------------       
%                           Experimental Loop
%--------------------------------------------------------------------------
startTime = clock;
for itrial = 1:trials
    trialColor = allRGB(colorOrder(itrial), :);
    
    trialCenter = yCenter+24;
    Screen('TextSize', window, 56);
    
    if trialOrder(itrial) == 0
        trialWord = cell2mat(datasample(otherWords, 1));

    else
        trialWord = allNames{wordOrder(itrial)};
    end
    
    DrawFormattedText(window, trialWord, 'center', trialCenter, trialColor);
    Screen('TextSize', window, 28)
    DrawFormattedText(window, '1=RED      2=YELLOW      3=GREEN      4=BLUE', 'center', yCenter+400, white);
    timeStamp(itrial,1) = Screen('Flip', window);
    
% Record response, 
    
fixflickerB=0;
    %Number - color matching 
    response = NaN;
    RT = NaN;
    Rcounter=0;
    tic;
    while toc < timeLimit
        if Rcounter==0;
        [~,~, keyCode] = KbCheck;
        if keyCode(oneKey)
            response = 1; %Red
            RT = toc;
            Rcounter=1;
        elseif keyCode(twoKey)
            response = 2; %Yellow
            RT = toc;
            Rcounter=1;
        elseif keyCode(threeKey)
            response = 3; %Green
            RT = toc;
            Rcounter=1;
        elseif keyCode(fourKey)
            response = 4; %Blue
            RT = toc;
            Rcounter=1;
        end
        end
        if toc > stimTime 
            if fixflickerB == 0
            Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter, yCenter], 2);
            Screen('TextSize', window, 28)
            DrawFormattedText(window, '1=RED      2=YELLOW      3=GREEN      4=BLUE', 'center', yCenter+400, white);
            timeStamp(itrial,2) = Screen('Flip', window);
            fixflickerB = 1;
            end 
        end
    end
    
    %Accuracy 
    acc = response == colorOrder(itrial);
    
    allData(itrial, :) = {trialWord, colorOrder(itrial), trialOrder(itrial), response, RT, acc, wordOrder(itrial),timeStamp};
    
    %Block Break 
       if mod(itrial, trialsperblock) == 0
        blockRT = nanmean(cell2mat(allData(trialsperblock*(block-1)+1:trialsperblock*block, 5)));
        blockAcc = nanmean(cell2mat(allData(trialsperblock*(block-1)+1:trialsperblock*block, 6)));
        Accmessage = ['This block, you got ' num2str(blockAcc*100) '% correct.'];
        RTmessage = ['Your average reaction time was ' num2str(blockRT) ' seconds.'];
        DrawFormattedText(window, Accmessage, 'center', yCenter-100, white);
        DrawFormattedText(window, RTmessage, 'center', yCenter+100, white);
        Screen(window, 'Flip')
        WaitSecs(5);
        block = block + 1;
        
        %Fixation
Screen('DrawLines', window, allCoords,...
            lineWidthPix, white, [xCenter, yCenter], 2);
            Screen('TextSize', window, 28)
            DrawFormattedText(window, '1=RED      2=YELLOW      3=GREEN      4=BLUE', 'center', yCenter+400, white);
            Screen('Flip', window);
        WaitSecs(ITI);
   end
    
end
endTime= clock;


%--------------------------------------------------------------------------
%                          Save Data
%--------------------------------------------------------------------------
allData{1, 9} = startTime(4:6);%Start time
allData{1, 10} = startTime(1:3); %Date
allData{trials, 9} = endTime(4:6); %End time 
        
save(fullfile(taskFolder, ['run' num2str(run) '.mat']), 'allData');


sca;
end


