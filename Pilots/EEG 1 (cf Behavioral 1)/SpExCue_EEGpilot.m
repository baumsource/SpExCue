function SpExCue_EEGpilot(ID,varargin)
%SpExCue_runEEGpilot - Experimental routine for EEG pilot experiment of the 
%SpExCue project testing the effect of spectral magnitude compression of HRTFs
%on sound externalization 
%
% Key-value pairs:
% azi = [30,0,-30]; % set of azimuths
% ele = [0,30,0]; % set of elevations
% M = [1,0.5,0]; % set of spectral saliences
% Nrep = 3; % repetitions
% flow = 700; % lower cut-off frequency
% fhigh = 18000; % upper cut-off frequency
% dur = 1.5; % duration of stimulus pair in sec
% fadeDur = 0.05; % duration of fade-in/out in sec
% jitter = 0.1; % temporal jitter in sec
% SPL = 70; % SPL in dB
% SPLrove = 10; % roving of SPL in dB
% rdepth = 10; % depth (peak-to-peak) in dB of spectral magnitude ripple
% rdensity = 0.25; % density of spectral magnitude ripple
% screenNumber = 1; % For 3rd floor lab (when using dual screens) choose #1, otherwise use #0
%
% Flags:
% stimulation = {'binaural','monaural'};
% roving = {'componentRove','pairRove'};
% feedback = {'TbTfeedback','blockedFeedback','noFeedback'}; % TbTfeedback:
% feedback on trial-by-trial basis by changing fixation dot color;
% blockedFeedback: provide only summary scores after every block
% tdt = {'TDTon','TDToff'};
% psychtoolbox = {'','debugMode'}; % to generate a very small (useless) Psychtoolbox window but allow access to Matlab command window while executing script
% HRTFs = {'HRCeq','HRC3ms','ARI'};

% AUTHOR: Robert Baumgartner

saveFileNamePrefix = 'SpExCue_EEGpilotResult';

%% Experimental variables
definput.keyvals.fnExtension = ''; % filename extension
definput.keyvals.azi = 30;%[30,0,-30]; % set of azimuths
definput.keyvals.ele = 0;%[0,30,0]; % set of elevations
definput.keyvals.M = [1,0.5,0]; % set of spectral saliences
definput.keyvals.Nrep = 100; % repetitions
definput.keyvals.flow = 700; % lower cut-off frequency
definput.keyvals.fhigh = 18000; % upper cut-off frequency
definput.keyvals.dur = 1.2; % duration of stimulus pair in sec
definput.keyvals.fadeDur = 0.05; % duration of fade-in/out in sec
definput.keyvals.jitter = 0.1; % temporal jitter in sec
definput.keyvals.SPL = 80; % SPL in dB
definput.keyvals.SPLrove = 0; % roving of SPL in dB
definput.keyvals.rdepth = 20; % depth (peak-to-peak) in dB of spectral magnitude ripple
definput.keyvals.rdensity = 0.25; % density of spectral magnitude ripple
definput.keyvals.screenNumber = 1; % For 3rd floor lab (when using dual screens) choose #1, otherwise use #0
definput.flags.stimulation = {'binaural','monaural'};
definput.flags.roving = {'componentRove','pairRove'};
definput.flags.familiarization = {'familiarize','skipFamiliarization'};
definput.flags.feedback = {'blockedFeedback','TbTfeedback','noFeedback'}; % Give feedback on trial-by-trial basis by changing fixation dot color
definput.flags.tdt = {'TDTon','TDToff'};
definput.flags.psychtoolbox = {'','debugMode'}; % to generate a very small (useless) Psychtoolbox window but allow access to Matlab command window while executing script
definput.flags.HRTFs = {'HRCeq','HRC3ms','ARI'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

%% Enter listener ID
if not(exist('ID','var'))
  subj.ID = upper(input('Enter subject ID: ','s'));
else
  subj.ID = ID;
end

%% Save path
savename = fullfile('Results',[saveFileNamePrefix '_' subj.ID,kv.fnExtension]);
if not(exist('Results','dir'))
  mkdir('Results')
end

%% Internal seetings
if flags.do_TDTon
  closerKey = 67; % C key
  fartherKey = 70; % F key
else
  closerKey = 6; % C key
  fartherKey = 9; % F key
end

trigVals = struct(... 
  'correctResponse',8,...
  'wrongResponse',9,...
  'left',000,...
  'top',100,...
  'right',200,...
  'M1_1',30,... % M1 (at stimulus onset) of 1
  'M1_p5',20,...
  'M1_0',10,...
  'M2_1',3,... % M2 (at stimulus change) of 1
  'M2_p5',2,...
  'M2_0',1,...
  'stimulusOffset',9);

% direction = round(trig/100)
% M1 = round(mod(trig,100)/10)
% M2 = mod(trig,10)
% specific D -> combination of one or two from above

%% Check availability of dependent functions
if not(exist('SpExCue_stim','file'))
  addpath(fullfile('..','MATLAB_general'))
end
if not(exist('fftreal','file'))
  amtstart
end

%% Initialize the TDT and playback system
fsVal = 48; % sampling rate, 48 stands for 48828.125 Hz (-> max buffer length 85 seconds)
minmaxChVolt = 1.0;  % min/max voltage for the output channels (scaling)
trigDuration = 0.001;  % duration of trigger in seconds

if flags.do_TDTon
  myTDT = tdt('playback_2channel', fsVal, minmaxChVolt, trigDuration );
end

%% Initialize the graphical interface for Psychtoolbox
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

%% Listener instruction
instruction1 = [...
  'In this experiment you will hear pairs of two different sounds.\n',... 
  'Only focus on the SPATIAL LOCATION of the sounds\n',...
  'and try to ignore other differences like intensity or timbre.\n',...
  'Sometimes it might sound as coming from a specific location outside your head and\n',...
  'sometimes it might sound as coming from close to your skin or even inside your head.\n',...
  '\n',...
  'During the experiment you are asked to perform the following task:\n',...
  'Press the *C* key if the SECOND sound appears to be CLOSER than the first sound or\n',...
  'press the *F* key if the SECOND sound appears to be FARTHER than the first sound.\n',...
  'Keep your eyes focused on the centered dot and respond when the dot turned blue (after sound turned off)!\n',...
  '\n',...
  'Please try not to move during sound presentation!\n',...
  'You will have the opportunity to take breaks and move after blocks of less than one minute.\n',...
  '\n',...
  'Press any key for some examples before starting the actual experiment!\n',...
  'If these exemplary sounds are incomfortably loud, please tell the experimenter immediately!\n'];
DrawFormattedText(win,instruction1,.2*x_center,'center',white,120,0,0,1.5);
Screen('Flip',win);

%% Stimulus generation and position permutation
if length(kv.azi) > 1 && length(kv.ele) == 1
  kv.ele = kv.ele*ones(length(kv.azi),1);
end
pos = [kv.azi(:),kv.ele(:)];
Npos = length(kv.azi);
pos = pos(randperm(Npos),:);
if flags.do_TDTon
  fs = myTDT.sampleRate;
else
  fs = 48.8e3;
end
subj.stim = SpExCue_stim( kv.M,subj.ID,pos,round(fs),kv.flow,kv.fhigh,kv.SPL,flags.HRTFs );

% Monaural?
if flags.do_monaural
  for ii = 1:numel(subj.stim.sig)
    ILD(ii) = -diff(dbspl(subj.stim.sig{ii}));
    if ILD(ii) > 0
      subj.stim.sig{ii}(:,2) = 0;
    else
      subj.stim.sig{ii}(:,1) = 0;
    end
  end
  savename = [savename '_mon'];
end

% define all possible combinations of two different M
C = nchoosek(kv.M,2);
Mcomb = [C;fliplr(C)]; % combinations of M
C = nchoosek(1:length(kv.M),2);
iM = [C;fliplr(C)]; % corresponding combinations of indices of M
Ncomb = size(Mcomb,1);
Ntotal = Npos*kv.Nrep*Ncomb;

%% Listener familiarization & SPL check
if flags.do_familiarize
    KbStrokeWait;
    % Fixation point
    Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
    Screen('Flip',win);
    pause(1)
    for pp = 1:size(subj.stim.sig,2)
      for ii = 1:size(subj.stim.sig,1)
        i1 = ii;
        i2 = mod(ii+1,length(kv.M))+1;
        sigpair = SpExCue_crossfade(subj.stim.sig{i1,pp},subj.stim.sig{i2,pp},...
          subj.stim.fs,kv.dur,kv.dur/2,kv.fadeDur);
        sigpair = 10^(kv.SPLrove/2/20)*sigpair; % present max level
        if flags.do_TDTon
          myTDT.load_stimulus(sigpair);
          myTDT.play()
        else
          sound(sigpair,fs)
        end
        
        % response after stimulus offset
        pause(kv.dur)
        Screen('DrawDots',win, [x_center,y_center], 14, blue, [], 2);
        Screen('Flip',win);
      
        % response via keyboard 
        keyCodeVal = 0;
        while not(keyCodeVal==closerKey || keyCodeVal==fartherKey) % 67...C, 70...F
            [tmp,keyCode] = KbWait([],2);
            keyCodeVal = find(keyCode,1);
        end
        E = sign(keyCodeVal-(closerKey+fartherKey)/2); 
        D = kv.M(i2)-kv.M(i1);
        hit = E*D > 0;

        % feedback
        if flags.do_TbTfeedback
          if hit
            Screen('DrawDots',win, [x_center,y_center], 14, green, [], 2);
          else
            Screen('DrawDots',win, [x_center,y_center], 14, red, [], 2);
          end
        else
          Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
        end
        Screen('Flip',win);
        pause(0.4)
        Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
        Screen('Flip',win);
        
        pause(1)
      end
    end
    instruction2 = [...
      'These were the test examples.\n',...
      'Do you have any questions?\n',...
      'If not press any key to start!\n'];
    DrawFormattedText(win,instruction2,'center','center',white,120,0,0,1.5);
    Screen('Flip',win);
end

%% Test procedure
KbStrokeWait;

% initialize variables
ii = 0; % incremental counter
subj.Mcomb = nan(Ntotal,2); % M combination
subj.D = nan(Ntotal,1); % D = M2-M1
if flags.do_componentRove
  subj.SPL = nan(Ntotal,2); % SPL
  subj.rphase = nan(Ntotal,2); % spectral ripple phase
else 
  subj.SPL = nan(Ntotal,1); % SPL
  subj.rphase = nan(Ntotal,1); % spectral ripple phase
end
subj.E = nan(Ntotal,1); % relative externalization response: -1...closer, +1...farther
subj.RT = nan(Ntotal,1); % reaction time
NperBlock = 2; % runs per block
skipBreak = false; % flag for skipping intermediate breaks within blocks
for irep = 1:kv.Nrep
  
  if not(skipBreak)
    if flags.do_TDTon
      send_event(myTDT, 254) % per my .cfg file for the Biosemi, a trigger of 254 unpauses and 255 pauses the EEG recording
    end
    pause(1)  % add a pause to ensure the pause is active before the sound 
  end
  
  for jj = 1:Npos
    
    % Show position
    if pos(jj,1) > 0
      aziLabel = 'left';
      aziTrigVal = trigVals.left;
    elseif pos(jj,1) < 0
      aziLabel = 'right';
      aziTrigVal = trigVals.right;
    else
      aziLabel = 'front';
      aziTrigVal = trigVals.top;
    end
    if Npos > 1
      DrawFormattedText(win,['Sounds from ' aziLabel '.'],'center','center',white);
      Screen('Flip',win);
      pause(2)
    end

    % Fixation point
    Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
    Screen('Flip',win);
    pause(1)
    
    % Randomized presentation order
    for iC = randperm(Ncomb);
      ii = ii+1;
      subj.pos(ii,:) = [kv.azi(jj),kv.ele(jj)];
      subj.Mcomb(ii,:) = Mcomb(iC,:);
      subj.D(ii) = diff(Mcomb(iC,:));
      sig1 = subj.stim.sig{iM(iC,1),jj};
      sig2 = subj.stim.sig{iM(iC,2),jj};
      M1TrigVal = (2*Mcomb(iC,1)+1)*10;
      M2TrigVal = 2*Mcomb(iC,2)+1;

      % roving of stimulus components
      if flags.do_componentRove
        % spectral roving
        subj.rphase(ii,:) = 2*pi*rand(1,2);
        sig1 = SpExCue_spectralRipple(sig1,subj.stim.fs,...
          kv.rdepth,kv.rdensity,subj.rphase(ii,1),kv.flow,kv.fhigh,'oct');
        sig2 = SpExCue_spectralRipple(sig2,subj.stim.fs,...
          kv.rdepth,kv.rdensity,subj.rphase(ii,2),kv.flow,kv.fhigh,'oct');
        % level roving
        dSPL = kv.SPLrove*(rand(1,2)-0.5);
        subj.SPL(ii,:) = kv.SPL + dSPL;
        sig1 = 10^(dSPL(1)/20)*sig1;
        sig2 = 10^(dSPL(2)/20)*sig2;
      end

      % combine stimulus pairs with temporal jitter of crossfade
      dt = kv.jitter*(rand-0.5);
      [sigpair,nM2] = SpExCue_crossfade(sig1,sig2,...
        subj.stim.fs,kv.dur,kv.dur/2+dt,kv.fadeDur);

      % roving of stimulus pair
      if flags.do_pairRove
        % spectral roving
        subj.rphase(ii) = 2*pi*rand;
        sigpair = SpExCue_spectralRipple(sigpair,subj.stim.fs,...
          kv.rdepth,kv.rdensity,subj.rphase(ii),kv.flow,kv.fhigh,'oct');
        % level roving
        dSPL = kv.SPLrove*(rand-0.5);
        subj.SPL(ii) = kv.SPL + dSPL;
        sigpair = 10^(dSPL/20)*sigpair;
      end
      
      if flags.do_debugMode
        figSgram = figure;
        subplot(1,2,1)
        sgram(sigpair(:,1),fs,'dynrange',60,'db')
        title('Left')
        subplot(1,2,2)
        sgram(sigpair(:,2),fs,'dynrange',60,'db')
        title('Right')
      end
      
      % playback
      if flags.do_TDTon
        triggerInfo = [   1,aziTrigVal+M1TrigVal; ...
                        nM2,aziTrigVal+M1TrigVal+M2TrigVal;...
                        size(sigpair,1),trigVals.stimulusOffset];
        myTDT.load_stimulus(sigpair, triggerInfo);
        tic;
        myTDT.play()
      else
        tic;
        sound(sigpair,fs)
      end
      
      % response after stimulus offset
      pause(kv.dur)
      Screen('DrawDots',win, [x_center,y_center], 14, blue, [], 2);
      Screen('Flip',win);
      
      % response via keyboard 
      keyCodeVal = 0;
      while not(keyCodeVal==closerKey || keyCodeVal==fartherKey) % 67...C, 70...F
          [tmp,keyCode] = KbWait([],2);
          keyCodeVal = find(keyCode,1);
      end
      subj.RT(ii) = toc-kv.dur/2+dt;
      subj.E(ii) = sign(keyCodeVal-(closerKey+fartherKey)/2); 
      subj.hit(ii) = subj.E(ii)*subj.D(ii) > 0;
      
      if flags.do_TDTon
        if subj.hit(ii)
          send_event(myTDT, trigVals.correctResponse)
        else
          send_event(myTDT, trigVals.wrongResponse)
        end
      end
      
      % feedback
      if flags.do_TbTfeedback
        if subj.hit(ii)
          Screen('DrawDots',win, [x_center,y_center], 14, green, [], 2);
        else
          Screen('DrawDots',win, [x_center,y_center], 14, red, [], 2);
        end
      else
        Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
      end
      Screen('Flip',win);
      pause(0.4)
      Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
      Screen('Flip',win);

      pause(0.5 + kv.jitter*(rand-0.5))
      
      if flags.do_debugMode
        close(figSgram)
      end
      
    end

  end
  
  if not(mod(irep,NperBlock)) || irep == kv.Nrep
    % Pause TDT
    if flags.do_TDTon
      send_event(myTDT, 255)  % per my .cfg file for the Biosemi, a trigger of 254 unpauses and 255 pauses the EEG recording
      myTDT.reset();  % rewind and clear buffers
    end

    % Save results
    save(savename,'subj','kv','flags')

    % Display time course and intermediate score
    infotext = [num2str(irep/NperBlock) ' of ' num2str(ceil(kv.Nrep/NperBlock)) ' blocks completed.'];
    if not(flags.do_noFeedback)
      pcorrect = 100*sum(subj.hit(1:ii))/ii;
      infotext = [infotext,'\n\n\n',num2str(pcorrect,'%3.2f'),'% correct.'];
    end
    if irep == kv.Nrep/2
      infotext = [infotext,'\n\n\n','Great! You already finished the first half.\n'];
      infotext = [infotext,'Please knock on the door of the booth!'];
    else
      infotext = [infotext,'\n\n\n','Press any key to continue!'];
    end
    disp(infotext)
    DrawFormattedText(win,infotext,'center','center',white);
    Screen('Flip',win);
    KbStrokeWait;
    skipBreak = false;
  else
    skipBreak = true;
  end
  
end

% calculate performance statistics

%% Inform listener that experiment is completed
i1 = 'Thank you! The experiment is completed.';
i2 = '\n\n\n\n Press any key to close this window!';

DrawFormattedText(win,[i1, i2],'center','center',white);
Screen('Flip',win);
KbStrokeWait;
Screen('CloseAll');

end