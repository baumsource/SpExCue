function bonus = SpExCue_BaDaGa(subID,Nrep,fnext,flags,kv)
% SpExCue_BaDaGa - experimental procedure for syllable segregation with
% different spatial cues
% 
% Usage:  SpExCue_BaDaGa(subID,Nrep,fnext,flags,kv)
% 
% Input parameters:
%   subID :   subject ID
%   Nrep  :   # repetitions per condition
%   fnext :   filename extension (use: train or test)
%   flags :   flags, options see below
%   kv    :   key-value pairs, options see below
%
%   flags:
%     flags.do_debugMode = 0; % set to 1 in order to show subject screen only in a small window 
%     flags.do_bonus = 1;     % set to 1 in order to include bonus payment of 1 cent per correct trial
%     flags.do_feedback = 1;  % set to 1 in order to display feedback screen
%
%   kv:
%     kv.screenNumber = 1; % choose screen number for psychtoolbox
%
%   Requirements:
%     1)  SOFA API

% AUTHOR: Robert Baumgarnter


hrtfPath = []; % if empty genStimuBaDaGa will use 'C:\Experiments\Robert\Experiments_GIT\HRTFs\';
trialsPerBlock = 50; % lasts 6.25 minutes
SPL = 75; % dB
respDuration = 4000e-3; % sec
feedbackDuration = 400e-3; % sec
angularEccentricity = 30; %deg

fs=2*24414; 

trgVal.right = 251;
trgVal.wrong = 252;
trgVal.na = 253;
trgVal.unpause = 254; % per my .cfg file for the Biosemi, a trigger of 254 unpauses and 255 pauses the EEG recording
trgVal.pause = 255;

savename = fullfile('data',[subID,'_',num2str(angularEccentricity),'deg_',fnext]);
ntest = 0;
while exist([savename,'.mat'],'file') % if it already exists, increase test number  
  ntest = ntest+1;
  savename = [savename,'_',num2str(ntest)];
end

if IsOSX
  flags.do_TDTon = false;
else
  flags.do_TDTon = true;
end

%% Stimulus creation

stimu = genStimuBaDaGa(fs,angularEccentricity*[-1 1],Nrep,subID,hrtfPath,SPL );
targetAngles = unique([stimu.targetlocation]);
trialnumber = length(stimu); 
spat = fliplr(unique({stimu.spatialization}));
% randomization
neworder=randperm(trialnumber);
randstimu = stimu(neworder);

%% Initialize the TDT and playback system
fsVal = floor(fs/1e3); % sampling rate, 48 stands for 48828.125 Hz (-> max buffer length 85 seconds)
minmaxChVolt = 1;%db2mag(3*3)*2.53;  % min/max voltage for the output channels (scaling), 2.53 yields 78db SPL at -3dB headphone amplification
trigDuration = 0.005;  % duration of trigger in seconds

if flags.do_TDTon
  myTDT = tdt('playback_2channel', fsVal, minmaxChVolt, trigDuration );
  if floor(myTDT.sampleRate/1e3) ~= fsVal
    error('RB: Sampling rate issue. Reboot PC and TDT!')
  end
end

%% Screen setup
PsychDefaultSetup(1); % makes sure Screen is functional and unifies keyCodes across OS
HideCursor;

screens = Screen('Screens');

white = [255 255 255]; %WhiteIndex(screenNumber);
black = [0 0 0]; %BlackIndex(screenNumber);
blue = [0 0 255];
green = [0 255 0];
red = [255 0 0];

if flags.do_debugMode
    debugRect = [10 10 200 200];
    [win,winRect] = Screen('OpenWindow',kv.screenNumber,black,debugRect);
    ShowCursor;
else
    [win,winRect] = Screen('OpenWindow',kv.screenNumber, black);
end

[screenXpix,screenYpix] = Screen('WindowSize', win);
[x_center,y_center] = RectCenter(winRect);


%% Listener instructions
instruction = [...
'You will hear two simultaneous sequences from different voices, one to the left and one to the right.\n',... 
'The different voices can switch locations from trial to trial. \n',... 
'Before each trial, a different third voice will cue the position you should attend to.\n',... 
'The syllables /ba/, /da/, and /ga/ are combined to a random sequence of these syllables.\n',... 
'Repetitions of the same syllable can occur.\n',...
'Your task is to report back the sequence presented at the cued position by\n\n',... 
'     pressing key *1* for /ba/, *2* for /da/, and *3* for /ga/.\n\n',... 
'Enter your response within the time you see a second circle around the center dot.\n',... 
'After your response, the center dot will turn green if your response is correct and red if it is incorrect.\n',... 
'Keep your eyes focused on the center dot and try not to move especially during sound presentation.\n',...
'Press any key to continue.\n']; 

DrawFormattedText(win,instruction,.2*x_center,'center',white,120,0,0,1.5);
Screen('Flip',win);
KbStrokeWait;

%% Trial presentation
correctcount=zeros(length(spat),1);
NblocksTot = ceil(trialnumber/trialsPerBlock);
bonus = 0;
for itrial=1:trialnumber
  
    disp(itrial)
    nblock = ceil(itrial/trialsPerBlock);

    if nblock > ceil((itrial-1)/trialsPerBlock)
      if flags.do_TDTon
        send_event(myTDT, trgVal.unpause) 
        pause(1);
      end
      iSessionSentence = ['Block # ' num2str(nblock) ' of ' num2str(NblocksTot),'\n\n',...
        'Press any key to start this block.'];
      DrawFormattedText(win,iSessionSentence,.2*x_center,'center',white,120,0,0,1.5);
      Screen('Flip',win);
      KbStrokeWait;
    end

    %% Stimulation frame

    % Fixation point
    Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
    Screen('Flip',win);
    pause(1)

    %   Sound presentation
    if flags.do_TDTon
      myTDT.load_stimulus(randstimu(itrial).sound, randstimu(itrial).triggers);
      myTDT.play_blocking()
    else
      sound(randstimu(itrial).sound,fs)% Response via keyboard 
      pause(length(randstimu(itrial).sound)/fs)
    end
    if flags.do_debugMode
      disp((randstimu(itrial).triggers))
      disp(randstimu(itrial).spatialization)
    end
  

    %% Response frame
    Screen('DrawDots',win, [x_center,y_center], 14, blue, [], 2);
    Screen('Flip',win);
    
    answer = NaN * ones(1,3);
    t = 0;
    btncount = 0;
    status = 0;
    t0resp = tic;
    while ( t < respDuration)
        [keyIsDown, resptime, keyCode] = KbCheck;
        if keyIsDown
            if status == 0
                status = 1;
            else
                continue
            end
            btncount = btncount + 1;
            cc = KbName(keyCode);  % find out which key was pressed
            answer(btncount) = str2double(cc(1));
        else
            status = 0;
        end
        if btncount >= 3
            t = toc(t0resp);
            break
        end
        t = toc(t0resp);
    end
    randstimu(itrial).response = answer(1:3); 
    randstimu(itrial).responsetime=t;
    
    %% feedback frame
    hrtfFlag = find(ismember(spat,randstimu(itrial).spatialization));
    points = ( randstimu(itrial).response == randstimu(itrial).seqtar);

    if isnan(answer)
%             renderVisFrame(PS,'NORESP');
        Screen('DrawDots',win, [x_center,y_center], 14, black, [], 2);
        ifbtrg = trgVal.na;
    else
        if ( sum(points) ==3 )
%                 renderVisFrame(PS,'CORRECT');
            Screen('DrawDots',win, [x_center,y_center], 14, green, [], 2);
            correctcount(hrtfFlag)=correctcount(hrtfFlag)+1;
            ifbtrg = trgVal.right;
        else
%                 renderVisFrame(PS,'WRONG');
            Screen('DrawDots',win, [x_center,y_center], 14, red, [], 2);
            ifbtrg = trgVal.wrong;
        end
    end
        
    if flags.do_TDTon
      send_event(myTDT, ifbtrg)
    end
    if flags.do_debugMode
      disp(ifbtrg)
    end
    if flags.do_feedback
      Screen('Flip',win);
      pause(feedbackDuration);
    end
    
    %% clean up

  if nblock < ceil((itrial+1)/trialsPerBlock) || itrial == trialnumber
    
    endSentence=['Block #' num2str(nblock),' completed.\n'];
    if flags.do_bonus
      bonus = sum(correctcount)/100;
      endSentence=[endSentence,'\n','Bonus earned: $',num2str(bonus,'%1.2f')];
    end
    DrawFormattedText(win,endSentence,.2*x_center,'center',white,120,0,0,1.5);
    Screen('Flip',win);
    disp(endSentence)
    
    % Pause TDT
    if flags.do_TDTon
      pause(1)
      send_event(myTDT, trgVal.pause) 
      myTDT.reset();  % rewind and clear buffers
    end
    
    KbStrokeWait;

    lastBlock = nblock;
    
    result=randstimu(1:itrial);
    save(savename, 'result');
  end


end

sca;

NperSpat = zeros(length(spat),1);
pc = zeros(length(spat),1);
for ispat = 1:length(spat)
  NperSpat(ispat) = sum(ismember({randstimu(1:itrial).spatialization},spat{ispat}));
  pc(ispat) = 100*correctcount(ispat)/NperSpat(ispat);
end
pcTab = array2table(pc,'VariableNames',{'percentCorrect'},'RowNames',spat);
disp(pcTab)

end