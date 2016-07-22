function SpExCue_EEGpilot3(ID,varargin)
%SpExCue_EEGpilot3 - Experimental routine for EEG pilot experiment of the 
%SpExCue project testing the effect of spectral magnitude compression of HRTFs
%on sound externalization 
%
% Key-value pairs:
% keyvals.fnExtension = ''; % filename extension
% keyvals.azi = 90;%[30,0,-30]; % set of azimuths
% keyvals.ele = 0;%[0,30,0]; % set of elevations
% keyvals.M = [1,1/3,0]; % set of spectral saliences
% keyvals.Nrep = 120; % repetitions
% keyvals.flow = 1e3; % lower cut-off frequency
% keyvals.fhigh = 16e3; % upper cut-off frequency
% keyvals.dur = 1.2; % duration of stimulus pair in sec
% keyvals.fadeDur = 0.05; % duration of fade-in/out in sec
% keyvals.xFadeDur = 0.01; % duration of cross-fade in sec
% keyvals.jitter = 0.1; % temporal jitter in sec
% keyvals.SPL = 70; % SPL in dB
% keyvals.SPLrove = 0; % roving of SPL in dB
% keyvals.rdepth = 10; % depth (peak-to-peak) in dB of spectral magnitude ripple
% keyvals.rdensity = 0.25; % density of spectral magnitude ripple
% keyvals.screenNumber = 1; % For 3rd floor lab (when using dual screens) choose #1, otherwise use #0
% keyvals.NperBlock = 20; % trials per block
%
% Flags:
% flags.HRTFs = {'HRC3ms','HRCeq','ARI'};
% flags.stimulation = {'binaural','monaural'};
% flags.Mrepetition = {'repeateM','changeM'};
% flags.positionOrderPermutation = {'completePosOrderPermutation','halfPosOrderPermutation'}; % randomized permutation of sequential order of positions either tested completely or only the first half
% flags.roving = {'componentRove','pairRove','noRoving'};
% flags.familiarization = {'familiarize','skipFamiliarization'};
% flags.feedback = {'blockedFeedback','TbTfeedback','noFeedback','consistencyFeedback','D0detectionFeedback'}; % Give feedback on trial-by-trial basis by changing fixation dot color
% flags.tdt = {'TDTon','TDToff'};
% flags.psychtoolbox = {'','debugMode'}; % to generate a very small (useless) Psychtoolbox window but allow access to Matlab command window while executing script

% AUTHOR: Robert Baumgartner

saveFileNamePrefix = 'EEG3pilot';

%% Experimental variables
definput = arg_SpExcue;
[flags,kv]  = ltfatarghelper({},definput,varargin);

%% Enter listener ID
if not(exist('ID','var'))
  subj.ID = upper(input('Enter subject ID: ','s'));
else
  subj.ID = ID;
end

%% Save path
savename = fullfile('data',[saveFileNamePrefix '_' subj.ID,kv.fnExtension]);
if not(exist('./data','dir'))
  mkdir('data')
end

%% Internal seetings
KbName('UnifyKeyNames');
closerKey = KbName('c');
fartherKey = KbName('f');
equalKey = KbName('space');

trigVals = struct(...
  'left',000,...
  'top',100,...
  'right',200,...
  'M1_1',30,... % M1 (at stimulus onset) of 1
  'M1_p5',20,...
  'M1_0',10,...
  'M2_1',3,... % M2 (at stimulus change) of 1
  'M2_p5',2,...
  'M2_0',1,...
  'closer',5,...
  'farther',6,...
  'equal',7,...
  'correctResponse',8,...
  'wrongResponse',9,...
  'stimulusOffset',0);

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
minmaxChVolt = db2mag(5*3)*2.53;  % min/max voltage for the output channels (scaling), 2.53 yields 78db SPL at -3dB headphone amplification
trigDuration = 0.005;  % duration of trigger in seconds

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
  'In this experiment you will hear pairs of sounds. Only focus on their changes in \n',... 
  'SPATIAL LOCATION and try to ignore other differences like intensity or timbre. \n',...
  'Sometimes it might sound as coming from a specific location outside your head and \n',...
  'sometimes it might sound as coming from close to your scalp or even inside your head. \n',...
  '\n',...
  'During the experiment your task is to:\n',...
  'Press the *C* key if the SECOND sound appears to be CLOSER than the first sound or\n'];
if flags.do_repeateM
  instruction1 = [instruction1,...
  'press the *F* key if the SECOND sound appears to be FARTHER than the first sound or\n',...
  'press the *SPACEBAR* if the sound does not change at all (catch trial).\n'];
else
  instruction1 = [instruction1,...,...
  'press the *F* key if the SECOND sound appears to be FARTHER than the first sound.\n'];
end
instruction1 = [instruction1,...
  '\n',...
  'For good EEG recordings, it is very important that you keep \n',...
  'your eyes focused on the centered white dot during sound presentation\n',...
  'and respond not before the dot turned blue, which happens after sound off-set.\n',...
  'Also, please try not to move during sound presentation! You will have the \n',...
  'opportunity to take breaks and move after blocks of less than one minute.\n',...
  '\n'];
if flags.do_consistencyFeedback
  instruction1 = [instruction1,...
    'After each block you will receive feedback on your consistency ratio definded as\n',... 
    'the absolute difference between the number of closer vs farther judgments divided by \n',... 
    'the total number of trials testing the same stimulus pair in both possible orders.\n',...
    '\n'];
elseif flags.do_D0detectionFeedback
  instruction1 = [instruction1,...
    'After each block you will receive feedback on how many catch trials you correctly detected\n',...
    'and how many trials you erroneously marked as catch trials (false alarms).\n',... 
    '\n'];
end
instruction1 = [instruction1,...
  'If the sounds are incomfortably loud, please tell the experimenter immediately!\n\n'];
if flags.do_skipFamiliarization
  instruction1 = [instruction1,'Press any key to start!\n'];
else
  instruction1 = [instruction1,'Press any key for some examples before starting the actual experiment!\n'];
end
DrawFormattedText(win,instruction1,.2*x_center,'center',white,120,0,0,1.5);
Screen('Flip',win);

%% Stimulus generation
if length(kv.azi) > 1 && length(kv.ele) == 1
  kv.ele = kv.ele*ones(length(kv.azi),1);
end
pos = [kv.azi(:),kv.ele(:)];
Npos = length(kv.azi);
% pos = pos(randperm(Npos),:);
if flags.do_TDTon
  fs = myTDT.sampleRate;
else
  fs = 48.8e3;
end
subj.stim = SpExCue_stim( kv.M,subj.ID,pos,round(fs),kv.flow,kv.fhigh,kv.SPL,flags.HRTFs,'signalDuration',kv.dur );

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

%% Define permutations

% define combinations of M
NM = length(kv.M);
C = nchoosek(1:NM,2);
iMvar = [C;fliplr(C)]; % set of all index pairs without repetition
iMrep = repmat(1:NM,[NM,1]);
iMtransp = transpose(iMrep);
iMrep = [iMrep(:),iMtransp(:)]; % set of all index pairs with repetition

if flags.do_changeM
  iM = iMvar;
  iM = repmat(iM,[kv.Nrep,1]); % repeat Nrep times
elseif flags.do_repeateM
  iM = [iMrep;repmat(iMvar,[NM-1,1])]; % combine sets such that M repetitions are as frequent as other M changes (D equally distributed)
  iM = repmat(iM,[round(kv.Nrep/NM),1]); % repeat Nrep/NM times
end
NperPos = length(iM); % 
Ntotal = NperPos*Npos;

% Randomization of M combinations
idrandM = randperm(NperPos);
iM = iM(idrandM,:);

% Randomization of positions
iPos = perms(1:length(kv.azi))'; % all permutations of position orders
NposPermOrder = size(iPos,2); % # order permutations possible with set of positions
iPos = iPos(:,randperm(NposPermOrder)); % randomize across experimental runs
if flags.do_halfPosOrderPermutation
  NposPermOrder = NposPermOrder/2;
  iPos = iPos(:,1:NposPermOrder);
end
NposPermTotal = numel(iPos); % length of position sequence (min # of blocks)
if kv.Nrep/NposPermOrder ~= round(kv.Nrep/NposPermOrder)
  error(['Number of repetitions must be a multiple of ',num2str(NposPermOrder),' to assure complete randomization of positions.'])
end

% Repetitions for various positions
Nblocks = Ntotal/kv.NperBlock;
while Nblocks ~= round(Nblocks) || ...
    Nblocks/NposPermTotal ~= round(Nblocks/NposPermTotal)
  kv.NperBlock = kv.NperBlock-1;
  Nblocks = Ntotal/kv.NperBlock;
  dispFlag = true; 
end
if exist('dispFlag','var')
  disp(['Number of trials per block reduced to ',num2str(kv.NperBlock),'.'])
end
% repeat iM sequence within blocks Npos times
idPosRep = reshape(1:NperPos,kv.NperBlock,[]); % incrementing index arranged in columns (blocks)
idPosRep = repmat(idPosRep,[Npos,1]); % order within each block is repeated for every position
iM = iM(idPosRep(:),:);
iPos = reshape(iPos,1,[]); % turn permutation matrix to row vector
iPos = repmat(iPos,[kv.NperBlock,1]); % position constant within block
iPos = iPos(:);
iPos = repmat(iPos,[Nblocks/NposPermTotal,1]);

% Set values
subj.Mcomb = kv.M(iM); % M combination
subj.D = diff(subj.Mcomb,1,2); % D = M2-M1
subj.pos = pos(iPos,:);
subj.SPL = kv.SPL + kv.SPLrove*(rand(Ntotal,2)-0.5); % SPL
subj.rphase = 2*pi*rand(Ntotal,2); % spectral ripple phase
if flags.do_pairRove 
  subj.SPL = subj.SPL(:,1);
  subj.rphase = subj.rphase(:,1);
end

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
          subj.stim.fs,kv.dur,kv.dur/2,kv.xFadeDur);
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

if Ntotal > 300 % force breaks if testing takes longer than 15 min
  flags.do_forceBreaks = true;
else
  flags.do_forceBreaks = false;
end

% initialize response variables
subj.E = nan(Ntotal,1); % relative externalization response: -1...closer, +1...farther
subj.RT = nan(Ntotal,1); % reaction time
subj.hit = nan(Ntotal,1); % hits

if flags.do_consistencyFeedback
  U = unique(abs(subj.D)); % absolute difference in M
  consistencyCounter = zeros(length(U),Npos); % consistency counter measureing test-retest reliability with collapsed order
end

ii = 0; % incremental counter
for bb = 1:Nblocks
  
  if flags.do_TDTon
    send_event(myTDT, 254) % per my .cfg file for the Biosemi, a trigger of 254 unpauses and 255 pauses the EEG recording
  end
  pause(1.5)  % add a pause to ensure the pause is active before the sound 
    
  % Positioning
  if subj.pos(ii+1,1) > 0
    aziLabel = 'left';
    aziTrigVal = trigVals.left;
  elseif subj.pos(ii+1,1) < 0
    aziLabel = 'right';
    aziTrigVal = trigVals.right;
  else
    aziLabel = 'front';
    aziTrigVal = trigVals.top;
  end
  if Npos > 1 % display position
    DrawFormattedText(win,['Sounds from ' aziLabel '.'],'center','center',white);
    Screen('Flip',win);
    pause(2)
  end

  % Fixation point
  Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
  Screen('Flip',win);
  pause(1)

  % Randomized presentation order
  for iC = 1:kv.NperBlock
    ii = ii+1;
    sig1 = subj.stim.sig{iM(ii,1),iPos(ii)};
    sig2 = subj.stim.sig{iM(ii,2),iPos(ii)};
    M1TrigVal = ceil(2*subj.Mcomb(ii,1)+1)*10;
    M2TrigVal = ceil(2*subj.Mcomb(ii,2)+1);
    TrigValOnset = aziTrigVal+M1TrigVal;
    TrigValChange = aziTrigVal+M1TrigVal+M2TrigVal;

    % roving of stimulus components
    if flags.do_componentRove
      % spectral roving
      if subj.D(ii) == 0 % make sure that noises are spectrally not too similar
        phaseDiff = abs(diff(subj.rphase(ii,:)));
        if phaseDiff < pi/4 || phaseDiff > 7/4*pi
          subj.rphase(ii,1) = mod(subj.rphase(ii,1) + pi,2*pi);
        end
      end
      sig1 = SpExCue_spectralRipple(sig1,subj.stim.fs,...
        kv.rdepth,kv.rdensity,subj.rphase(ii,1),kv.flow,kv.fhigh,'oct');
      sig2 = SpExCue_spectralRipple(sig2,subj.stim.fs,...
        kv.rdepth,kv.rdensity,subj.rphase(ii,2),kv.flow,kv.fhigh,'oct');
      % level roving
      dSPL = subj.SPL(ii,:) - kv.SPL;
      sig1 = 10^(dSPL(1)/20)*sig1;
      sig2 = 10^(dSPL(2)/20)*sig2;
    end

    % combine stimulus pairs with temporal jitter of crossfade
    dt = kv.jitter*(rand-0.5);
    [sigpair,nM2] = SpExCue_crossfade(sig1,sig2,...
      subj.stim.fs,kv.dur,kv.dur/2+dt,kv.xFadeDur);

    % roving of stimulus pair
    if flags.do_pairRove
      % spectral roving
      sigpair = SpExCue_spectralRipple(sigpair,subj.stim.fs,...
        kv.rdepth,kv.rdensity,subj.rphase(ii),kv.flow,kv.fhigh,'oct');
      % level roving
      dSPL = subj.SPL(ii) - kv.SPL;
      sigpair = 10^(dSPL/20)*sigpair;
    end

    if flags.do_debugMode
      % Plot spectrograms
      figSgram = figure;
      subplot(1,2,1)
      sgram(sigpair(:,1),fs,'dynrange',60,'db')
      title('Left')
      subplot(1,2,2)
      sgram(sigpair(:,2),fs,'dynrange',60,'db')
      title('Right')
      % Display trial values
      M = [subj.Mcomb(ii,1);subj.Mcomb(ii,2)];
      Trigger = [TrigValOnset;TrigValChange];
      RipplePhase = subj.rphase(ii,:)';
      SPL = subj.SPL(ii,:)';
      RowNames = {'Onset','Change'};
      table(M,Trigger,RipplePhase,SPL,'RowNames',RowNames)
    end

    % playback
    if flags.do_TDTon
      triggerInfo = [   1,TrigValOnset; ...
                      nM2,TrigValChange;...
                      size(sigpair,1),trigVals.stimulusOffset];
      myTDT.load_stimulus(sigpair, triggerInfo);
      tic;
      myTDT.play_blocking()
    else
      tic;
      sound(sigpair,fs)
      pause(kv.dur)
    end

    % Visual marker at stimulus offset for response request
    Screen('DrawDots',win, [x_center,y_center], 14, blue, [], 2);
    Screen('Flip',win);

    % Response via keyboard 
    keyCodeVal = 0;
    while not(any(keyCodeVal==[closerKey,fartherKey,equalKey]))
        [tmp,keyCode] = KbWait([],2);
        keyCodeVal = find(keyCode,1);
    end
    subj.RT(ii) = toc-kv.dur/2+dt;
    
    % Externalization response
    if keyCodeVal == closerKey
      subj.E(ii) = -1;
      Etrig = trigVals.closer;
    elseif keyCodeVal == fartherKey
      subj.E(ii) = 1;
      Etrig = trigVals.farther;
    elseif keyCodeVal == equalKey
      subj.E(ii) = 0;
      Etrig = trigVals.equal;
    else
      subj.E(ii) = nan;
      Etrig = [];
    end
    if flags.do_TDTon
      send_event(myTDT, Etrig)
    end
    
    % Relationship between D and E 
    subj.hit(ii) = subj.E(ii) == sign(subj.D(ii));
    if flags.do_TDTon % Send response evaluation trigger
      pause(2*trigDuration)
      if subj.hit(ii)
        send_event(myTDT, trigVals.correctResponse)
      else
        send_event(myTDT, trigVals.wrongResponse)
      end
    end

    % Feedback
    if flags.do_TbTfeedback
      if isnan(subj.hit(ii)) || subj.hit(ii) % positive feedback for catch trials
        Screen('DrawDots',win, [x_center,y_center], 14, green, [], 2);
      else
        Screen('DrawDots',win, [x_center,y_center], 14, red, [], 2);
      end
    elseif flags.do_consistencyFeedback
      iU = U==abs(subj.D(ii));
      consistencyCounter(iU,iPos(ii)) = consistencyCounter(iU,iPos(ii)) + sign(subj.D(ii))*subj.E(ii);
      Nchange = sum(subj.D(1:ii)~=0);
      subj.consistency = 100*sum(abs(consistencyCounter(:)))/Nchange;
    elseif flags.do_D0detectionFeedback
      iD0 = subj.D(1:ii) == 0;
      D0hitRate = 100* sum(subj.hit(iD0)) / sum(iD0);
      iDx = subj.D(1:ii) ~= 0;
      D0faRate = 100* sum(subj.E(iDx)==0) / sum(iDx);
    else
      Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
    end
    Screen('Flip',win);
    pause(0.4)
    Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
    Screen('Flip',win);

    % ITI jitter
    pause(0.5 + kv.jitter*(rand-0.5))

    if flags.do_debugMode
      close(figSgram)
    end
  end
  
  % Intermediate score
  subj.pcorrect = 100* nansum(subj.hit(1:ii)) / ii;

  if kv.NperBlock > 6
    % Save results
    save(savename,'subj','kv','flags')

    % Pause TDT
    if flags.do_TDTon
      pause(0.5)
      send_event(myTDT, 255)  % per my .cfg file for the Biosemi, a trigger of 254 unpauses and 255 pauses the EEG recording
      myTDT.reset();  % rewind and clear buffers
    end

    % Display time course and intermediate score
    infotext = [num2str(bb) ' of ' num2str(Nblocks) ' blocks completed.'];
    if flags.do_consistencyFeedback
      infotext = [infotext,'\n\n\n Consistency ratio: ',num2str(subj.consistency,'%3.2f'),'%'];
    elseif flags.do_D0detectionFeedback
      if not(isnan(D0hitRate))
        infotext = [infotext,'\n\n\n Catch trials: ',num2str(D0hitRate,'%3.2f'),...
          '% correctly detected with ',num2str(D0faRate,'%3.2f'),'% false alarms.'];
      end
    elseif not(flags.do_noFeedback)
      infotext = [infotext,'\n\n\n',num2str(subj.pcorrect,'%3.2f'),'% correct.'];
    end
    if flags.do_forceBreaks
      forcedBreak = 60; % seconds
      if bb == round(Nblocks/4)
        infotext = [infotext,'\n\n\n','First quarter done! Please knock on the door and take a break!\n'];
      elseif bb == round(Nblocks/2)
        infotext = [infotext,'\n\n\n','Great! You already finished the first half.\n'];
        infotext = [infotext,'Please knock on the door and take a break!\n'];
      elseif bb == round(Nblocks*3/4)
        infotext = [infotext,'\n\n\n','Almost done! Please knock on the door and take a break!\n'];
      else
        forcedBreak = 5;
      end
      disp(infotext)
      DrawFormattedText(win,infotext,'center','center',white);
      Screen('Flip',win);
      pause(forcedBreak)
    end
    infotext = [infotext,'\n\n\n','Press any key to continue!'];
    DrawFormattedText(win,infotext,'center','center',white);
    Screen('Flip',win);
    % Experimenter monitoring info
    disp('Info given to listener:')
    disp(infotext)
    KbStrokeWait;
  end
  
  % Experimenter monitoring info
  disp('Experimenter-only info:')
  disp(['Percent correct: ',num2str(subj.pcorrect),'%'])
  
end


%% Inform listener that experiment is completed
i1 = 'Thank you! The experiment is completed.';
i2 = '\n\n\n\n Press any key to close this window!';

DrawFormattedText(win,[i1, i2],'center','center',white);
Screen('Flip',win);
KbStrokeWait;
Screen('CloseAll');

end