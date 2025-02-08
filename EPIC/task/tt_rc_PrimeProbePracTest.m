function tt_rc_PrimeProbePracTest
%By DMS based on template code from MD
%384 trials with 96 trials per block each trial being 1,500ms Prime = 200ms blank = 33ms Probe = 200ms fix = 1067ms% 
sca; %Closes any open Psychtoolbox windows 
close all; %closes any other matlab windows
clearvars; %clear all variables 

rng('shuffle'); %shuffle RNG generator
%2=AA 4=BB 6=AB 8=BA; 1=YY 3=ZZ 5=YZ 7=ZY%
load('subtrialsGood.mat'); %Get trial orders
%--------------------------------------------------------------------------
%                        Participant Setup
%--------------------------------------------------------------------------
%Enter participant ID in command window
prompt = 'Participant #: ';
subID = input(prompt, 's');

current = pwd();
subFolder = [current '/ParticipantInfo/' subID];
%Create Participant folder (if new)
if ~isfolder(subFolder)
    mkdir(subFolder)
end
taskFolder = [subFolder '/PrimeProbe'];
%Make Flanker folder (if 1st run)
if ~isfolder(taskFolder)
    mkdir(taskFolder)
end

%Check if run exists to prevent accidental overwrite
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



%--------------------------------------------------------------------------
%                       Change Variables here
%--------------------------------------------------------------------------
trials = 384; %Number of total trials
blocks = 4;
trialsperblock = trials/blocks; %Number of tirals in each block

stimTime = .2; %Time stimli on screne
blankTime = .033; %Blank Screen Time
ITI = 1.067; %Secs between stimili
responseTime = stimTime+blankTime+stimTime+ITI; %Secs to respond
stimblank=stimTime+blankTime;
totstimTime=stimTime+blankTime+stimTime; %Total stimulus time

numTrialStims = 5; %[A B Y Z & blank]
stimNames={'A','B','Y','Z','blank'};


AKey = KbName('1!');
BKey = KbName('2@');
YKey = KbName('3#');
ZKey = KbName('4$');

allData = cell(trials, 7);

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

HideCursor();

Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%Screen Text 
Screen('TextFont', window, 'Helvetica');
Screen('TextSize', window, 36); %Font size for instructions

%Fixation cross
fixCrossDimPix = 15;
lineWidthPix = 3;
fixXCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
fixYCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
fixAllCoords = [fixXCoords; fixYCoords];
%--------------------------------------------------------------------------
%                            Stimuli Setup
%--------------------------------------------------------------------------

%Load arrrow pics and make textures 
pictures = zeros(trials, 1);
for img = 1:numTrialStims
    tmpImg = imread([char(stimNames(img)) '.png']);
    pictures(img) = Screen('MakeTexture', window, tmpImg);
end

[imgX, imgY, ~] = size(tmpImg);
arrowSize = .2; %mulitplier for img size


%Set positions

arrowRect = [xCenter-(imgY*arrowSize)/2; yCenter-(imgX*arrowSize)/2; ...
             xCenter+(imgY*arrowSize)/2; yCenter+(imgX*arrowSize)/2];

%Get the subject's stimorder%
for orderchecker=1:size(subtrials,2);
    if strcmp(subID,subtrials{61,orderchecker});
       stimOrder=subtrials{run,orderchecker};
    end
end
%Note the order of stimNames {'A','B','Y','Z','blank'} 1 is A, 2 is B, 3 is Y, 4 is Z
%2=AA 4=BB 6=AB 8=BA; 1=YY 3=ZZ 5=YZ 7=ZY (These numbers from subtrials map to the letter pairs)%
% {'A','B','Y','Z','blank'} (FirstStimOrder & SecondStimOrder are used to make the letter pair for a trial)%
%The numbers (1 to 8) from subtrials map to columns in FirstStimOrder and SecondStimOrder
% Note how in stimNames Y is 3 [A,B,Y,Z,blank] and in the stimOrder which comes from substrials YY=1. Element one for both FirstStimOrder and SecondStimOrder is 3
FirstStimOrder=[3,1,4,2,3,1,4,2];   %{'Y','A','Z','B','Y','A','Z','B'};
SecondStimOrder=[3,1,4,2,4,2,3,1];  %{'Y','A','Z','B','Z','B','Y','A'};
%--------------------------------------------------------------------------       
%                        Instructions/Key Bindings 
%--------------------------------------------------------------------------
DrawFormattedText(window, 'In this task, on each trial you will be shown a letter followed by another letter.', 'center', yCenter-100, white);
DrawFormattedText(window, 'Ignore the first letter and respond to the second letter:', 'center', yCenter, white);
DrawFormattedText(window, 'A = hit "1", B = hit "2"', 'center', yCenter+100, white);
DrawFormattedText(window, 'Y = hit "3", Z = hit "4"', 'center', yCenter+200, white);
DrawFormattedText(window, 'Respond with your Right hand: 1=Index, 2=Middle, 3=Ring, 4=Pinky.', 'center', yCenter+300, white);
DrawFormattedText(window, 'You should work as quickly and as accurately as possible.', 'center', yCenter+400, white);
DrawFormattedText(window, 'Press ANY KEY to start the practice trials', 'center', yCenter+500, white);


Screen('Flip', window)
WaitSecs(.5);
KbWait;
Screen('DrawLines', window, fixAllCoords,...
        lineWidthPix, white, [xCenter, yCenter], 2); %Draw fixation
Screen('Flip', window);
WaitSecs(ITI);

%--------------------------------------------------------------------------
%                             Practice
%--------------------------------------------------------------------------
pracTrialOrder=[6;5;2;7;4;3;8;1;2;3;4;5;6;7;8;1;2;7;8;5;6;1;4;3];
numPracTrials = 24;
startTime = clock;
block = 1;
for ptrial = 1:numPracTrials
    if pracTrialOrder(ptrial) == 1 || pracTrialOrder(ptrial) == 2 || pracTrialOrder(ptrial) == 3 || pracTrialOrder(ptrial) == 4;
        trialType = 'CONGRUENT';
    else
        trialType = 'INCONGRUENT';
    end
    
    if pracTrialOrder(ptrial) == 2 || pracTrialOrder(ptrial) == 4 || pracTrialOrder(ptrial) == 6 || pracTrialOrder(ptrial) == 8;
        trialSet = 'AB';
    else
        trialSet = 'YZ';
    end
    
    %Draw arrows to correct spots 
    Screen('DrawTexture', window, pictures(FirstStimOrder(pracTrialOrder(ptrial))), [], arrowRect); %The loop number itrial acts as the index of stimOrder and the value of stimOrder at that index is to get FirststimOrder value at an index equal to the value of stimOrder and this value is used as the index of pictures%
    Screen('Flip', window);
    
    %Response
    fixflicker = 0;
    fixflickerB = 0;
    fixflickerS = 0;
    response = NaN;
    RT = NaN;
    Rcounter=0;
    tic;
    while toc < responseTime
        if Rcounter==0;
        [~,~, keyCode] = KbCheck;
        if keyCode(AKey) 
            response = 'A'; 
            RT = toc;
            Rcounter=1;
        elseif keyCode(BKey)
            response = 'B'; 
            RT = toc;
            Rcounter=1;
        elseif keyCode(YKey)
            response = 'Y'; 
            RT = toc;
            Rcounter=1;
        elseif keyCode(ZKey)
            response = 'Z'; 
            RT = toc;
            Rcounter=1;
        end
        end
        if toc > stimTime 
            if fixflickerB == 0
                Screen('DrawTexture', window, pictures(5), [], arrowRect); %Blank Period%
                Screen('Flip', window);
                fixflickerB = 1;
            end 
        end
        if toc > stimblank
            if fixflickerS == 0
                Screen('DrawTexture', window, pictures(SecondStimOrder(pracTrialOrder(ptrial))), [], arrowRect); %Blank Period%
                Screen('Flip', window);
                fixflickerS = 1;
            end
        end
        if toc > totstimTime %DMS: Made this so that the fix comes on after the 3 stim events
            if fixflicker == 0
                Screen('DrawLines', window, fixAllCoords,...
                lineWidthPix, white, [xCenter, yCenter], 2); %Draw fixation
                Screen('Flip', window);
                fixflicker = 1; %Stops cross from flickering 
            end
        end
    end
    
    
   if pracTrialOrder(ptrial) == 2 || pracTrialOrder(ptrial) == 8 % AA or BA
       if strcmp(response, 'A') %A response 
           trialacc = 1;
       else 
           trialacc =0;
       end
   elseif pracTrialOrder(ptrial) == 4 || pracTrialOrder(ptrial) == 6 % BB or AB
       if strcmp(response, 'B') %B response 
           trialacc = 1;
       else
           trialacc = 0;
       end
   elseif pracTrialOrder(ptrial) == 1 || pracTrialOrder(ptrial) == 7 % YY or ZY
       if strcmp(response, 'Y') %Y response 
           trialacc = 1;
       else
           trialacc = 0;
       end
   elseif pracTrialOrder(ptrial) == 3 || pracTrialOrder(ptrial) == 5 % ZZ or YZ
       if strcmp(response, 'Z') %Z response 
           trialacc = 1;
       else
           trialacc = 0;
       end
   end
   %Accuracy 
    acc = trialacc==1;
    
    
    if acc == 1
        message = 'Correct!';
    else
        message = 'Incorrect';
    end
    
    
    Screen('TextSize', window, 36);
    DrawFormattedText(window, message, 'center', yCenter-100, white);
    DrawFormattedText(window, 'Press ANY KEY to continue', 'center', yCenter+300, white);
    Screen('TextSize', window, 28)
    DrawFormattedText(window, 'A=Key1 B=Key2 Y=Key3 Z=Key4', 'center', yCenter+400, white);
    Screen('Flip', window);
    WaitSecs(.5);
    KbWait;
    DrawFormattedText(window, 'A=Key1 B=Key2 Y=Key3 Z=Key4', 'center', yCenter+400, white);
    Screen('Flip', window);
    WaitSecs(1); 
    %Back to fix%
    Screen('DrawLines', window, fixAllCoords,...
    lineWidthPix, white, [xCenter, yCenter], 2); %Draw fixation
    Screen('Flip', window);
    WaitSecs(.5); 
end
    %Ending Prac Below%
    Screen('TextSize', window, 36);
    DrawFormattedText(window, 'Practice complete!', 'center', yCenter-100, white)
    DrawFormattedText(window, 'During the actual task, you will receive feedback at the end of a block.',...
                        'center', 'center', white)
    DrawFormattedText(window, 'Press ANY KEY to begin actual task', 'center', yCenter+300, white);
    Screen('Flip', window)
    WaitSecs(.5);
    KbWait
    Screen('Flip', window);
    WaitSecs(2);
%--------------------------------------------------------------------------       
%                           Experimental Loop
%--------------------------------------------------------------------------
startTime = clock;
block = 1;
for itrial = 1:trials
    if stimOrder(itrial) == 1 || stimOrder(itrial) == 2 || stimOrder(itrial) == 3 || stimOrder(itrial) == 4;
        trialType = 'CONGRUENT';
    else
        trialType = 'INCONGRUENT';
    end
    
    if stimOrder(itrial) == 2 || stimOrder(itrial) == 4 || stimOrder(itrial) == 6 || stimOrder(itrial) == 8;
        trialSet = 'AB';
    else
        trialSet = 'YZ';
    end
    
    %Draw pcis to correct spots 
    Screen('DrawTexture', window, pictures(FirstStimOrder(stimOrder(itrial))), [], arrowRect); %The loop number itrial acts as the index of stimOrder and the value of stimOrder at that index is to get FirststimOrder value at an index equal to the value of stimOrder and this value is used as the index of pictures%
    timeStamp(itrial,1) = Screen('Flip', window);
    
    %Response
    fixflicker = 0;
    fixflickerB = 0;
    fixflickerS = 0;
    response = NaN;
    RT = NaN;
    Rcounter=0;
    tic;
    while toc < responseTime
        if Rcounter==0;
        [~,~, keyCode] = KbCheck;
        if keyCode(AKey) 
            response = 'A'; 
            RT = toc;
            Rcounter=1;
        elseif keyCode(BKey)
            response = 'B'; 
            RT = toc;
            Rcounter=1;
        elseif keyCode(YKey)
            response = 'Y'; 
            RT = toc;
            Rcounter=1;
        elseif keyCode(ZKey)
            response = 'Z'; 
            RT = toc;
            Rcounter=1;
        end
        end
        if toc > stimTime 
            if fixflickerB == 0
                Screen('DrawTexture', window, pictures(5), [], arrowRect); %Blank Period%
                timeStamp(itrial,2) = Screen('Flip', window);
                fixflickerB = 1;
            end 
        end
        if toc > stimblank
            if fixflickerS == 0
                Screen('DrawTexture', window, pictures(SecondStimOrder(stimOrder(itrial))), [], arrowRect); %Blank Period%
               timeStamp(itrial,3) =  Screen('Flip', window);
                fixflickerS = 1;
            end
        end
        if toc > totstimTime %DMS: Made this so that the fix comes on after the 3 stim events
            if fixflicker == 0
                Screen('DrawLines', window, fixAllCoords,...
                lineWidthPix, white, [xCenter, yCenter], 2); %Draw fixation
                timeStamp(itrial,4) = Screen('Flip', window);
                fixflicker = 1; %Stops cross from flickering 
            end
        end
    end
    
    
   if stimOrder(itrial) == 2 || stimOrder(itrial) == 8 % AA or BA
       if strcmp(response, 'A') %A response 
           trialacc = 1;
       else 
           trialacc =0;
       end
   elseif stimOrder(itrial) == 4 || stimOrder(itrial) == 6 % BB or AB
       if strcmp(response, 'B') %B response 
           trialacc = 1;
       else
           trialacc = 0;
       end
   elseif stimOrder(itrial) == 1 || stimOrder(itrial) == 7 % YY or ZY
       if strcmp(response, 'Y') %Y response 
           trialacc = 1;
       else
           trialacc = 0;
       end
   elseif stimOrder(itrial) == 3 || stimOrder(itrial) == 5 % ZZ or YZ
       if strcmp(response, 'Z') %Z response 
           trialacc = 1;
       else
           trialacc = 0;
       end
   end
    
   allData(itrial, :) = {trialType, stimOrder(itrial), response, trialacc, RT,trialSet,timeStamp};
   
   %Feedback
   if mod(itrial, trialsperblock) == 0
        blockRT = nanmean(cell2mat(allData(trialsperblock*(block-1)+1:trialsperblock*block, 5)));
        blockAcc = nanmean(cell2mat(allData(trialsperblock*(block-1)+1:trialsperblock*block, 4)));
        Accmessage = ['This block, you got ' num2str(blockAcc*100) '% correct.'];
        RTmessage = ['Your average reaction time was ' num2str(blockRT) ' seconds.'];
        DrawFormattedText(window, Accmessage, 'center', yCenter-100, white);
        DrawFormattedText(window, RTmessage, 'center', yCenter+100, white);
        Screen(window, 'Flip')
        WaitSecs(5);
        block = block + 1;
        
        %Fixation
        Screen('DrawLines', window, fixAllCoords,...
            lineWidthPix, white, [xCenter, yCenter], 2); %Draw fixation
        Screen('Flip', window);
        WaitSecs(ITI);
   end
   
end  
endTime = clock;

%--------------------------------------------------------------------------
%                              Save Data
%--------------------------------------------------------------------------
allData{1, 8} = startTime(4:6); %Start time
allData{2, 8} = startTime(1:3); %Date
allData{trials, 6} = endTime(4:6); %#ok<NASGU> %End time

save(fullfile(taskFolder, ['run' num2str(run) '.mat']), 'allData');

sca;
end

