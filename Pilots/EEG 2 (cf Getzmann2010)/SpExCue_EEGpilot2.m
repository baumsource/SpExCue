function SpExCue_EEGpilot2(ID,varargin)
%SpExCue_runEEGpilot - Experimental routine for EEG pilot experiment of the 
%SpExCue project testing the effect of spectral magnitude compression of HRTFs
%on sound externalization 
%
% Pilot 2: replication of Getzmann & Lewald (2010) with various spectral
% compressions (M) and discontinous spatial changes -> task: left vs
% right sound at various M, onset stimulus always from front 
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
% feedback = {'feedback',''}; % Give feedback on trial-by-trial basis by changing fixation dot color
% tdt = {'TDTon','TDToff'};
% psychtoolbox = {'','debugMode'}; % to generate a very small (useless) Psychtoolbox window but allow access to Matlab command window while executing script

% AUTHOR: Robert Baumgartner

saveFileNamePrefix = 'SpExCue_EEGpilot2Result';

%% Experimental variables
definput.keyvals.azi = 0; % set of azimuths
definput.keyvals.ele = 0; % set of elevations
definput.keyvals.M = [1,.5,.25,0]; % set of spectral saliences
definput.keyvals.Nrep = 125; % repetitions
definput.keyvals.flow = 250; % lower cut-off frequency
definput.keyvals.fhigh = 18000; % upper cut-off frequency
definput.keyvals.signalDuration = 1.2; % duration of stimulus pair in sec
definput.keyvals.fadeDur = 0.015; % duration of fade-in/out in sec
definput.keyvals.jitter = 0; % temporal jitter in sec
definput.keyvals.changeTime = 0.7; % onset of stimulus change in sec
definput.keyvals.SPL = 70; % SPL in dB
definput.keyvals.SPLrove = 0; % roving of SPL in dB
definput.keyvals.rdepth = 00; % depth (peak-to-peak) in dB of spectral magnitude ripple
definput.keyvals.rdensity = 0.25; % density of spectral magnitude ripple
definput.keyvals.screenNumber = 1; % For 3rd floor lab (when using dual screens) choose #1, otherwise use #0
definput.keyvals.responseTrialRatio = 0.2;
definput.flags.stimulation = {'binaural','monaural'};
definput.flags.roving = {'componentRove','pairRove'};
definput.flags.familiarization = {'familiarize','skipFamiliarization'};
definput.flags.excludeStaticSounds = {'staticSounds','noStaticSounds'};
definput.flags.feedback = {'feedback',''}; % Give feedback on trial-by-trial basis by changing fixation dot color
definput.flags.tdt = {'TDTon','TDToff'};
definput.flags.psychtoolbox = {'','debugMode'}; % to generate a very small (useless) Psychtoolbox window but allow access to Matlab command window while executing script

[flags,kv]  = ltfatarghelper({},definput,varargin);

%% Internal seetings
% if flags.do_TDTon % lab PC
%   leftKey = 67;
%   rightKey = 70;
% else % own PC
%   leftKey = 6;
%   rightKey = 9;
% end
KbName('UnifyKeyNames');
lkey = KbName('LeftArrow');
rkey = KbName('RightArrow');

trigVals = struct(... 
  'center',130,...
  'left',110,...
  'right',120,...
  'M1',5,... % M of 1
  'Mp5',3,...
  'Mp25',2,...
  'M0',1,...
  'responseTrial',100,...
  'precursorOnset',6,...
  'correctResponse',7,...
  'wrongResponse',8,...
  'stimulusOffset',9,...
  'pause',255,...
  'unpause',254);

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

%% Enter listener ID
if not(exist('ID','var'))
  subj.ID = upper(input('Enter subject ID: ','s'));
else
  subj.ID = ID;
end

%% Create filename for saving
fn = fullfile('Results',saveFileNamePrefix);
if not(exist('Results','dir')) 
  mkdir('Results')
end
if flags.do_monaural
  fn = [fn '_mon'];
end
if kv.responseTrialRatio > 0
  fn = [fn '_RT'];
end
fn = [fn,'_',subj.ID];

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
instruction = [...
  'In this experiment you will hear noise stimuli.\n',... 
  'At first they will always come from the center and \n',...
  'then they will either move to the left or right hand side.\n',...
  '\n',...
  'During the experiment you are asked to determine the directional change\n',...
  'whenever you hear a short tone preceding the pairs of noise bursts:\n',...
  'If the second burst appears to the LEFT, press the LEFT button.\n',...
  'If the second burst appears to the RIGHT, press the RIGHT button.\n',...
  '\n',...
  'It is very important to respond as FAST and ACCURATE as you can.\n',...
  '\n',...
  'Press any key to hear some exemplary sounds before starting the experiment!\n',...
  'Please inform the experimenter if the loudness appears uncomfortable!\n'];
DrawFormattedText(win,instruction,'center','center',white);
Screen('Flip',win);

%% Stimulus generation and position permutation
pos = [kv.azi(:),kv.ele(:)];
if flags.do_TDTon
  fs = myTDT.sampleRate;
else
  fs = 48.8e3;
end
% subj.stim = SpExCue_stim( kv.M,subj.ID,pos,round(fs),kv.flow,kv.fhigh,kv.SPL,...
%   'fadeDuration',0,'ARI','AMnoise' ); % dimensions of sig field: M x pos
subj.stim = SpExCue_stim( kv.M,subj.ID,pos,round(fs),kv.flow,kv.fhigh,kv.SPL,...
  'fadeDuration',0,'ARI','AMnoise','signalDuration',kv.signalDuration,'moveLeftRight' ); % dimensions of sig field: M x L/R

% 100-ms 500-Hz tone precursor
t = 0:1/fs:.1-1/fs;
precursor = sin(2*pi*500*t);
precursor(end-10:end) = precursor(end-10:end).*cos(0:pi/20:pi/2).^2;
precursor = setdbspl(precursor,kv.SPL);
precursor = repmat(precursor(:),[1,2]);

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
end

%% Listener familiarization & SPL check
if flags.do_familiarize
    KbStrokeWait;
    % Fixation point
    Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
    Screen('Flip',win);
    pause(1)
    for ii = 1:2:numel(subj.stim.sig)
      if flags.do_TDTon
        myTDT.load_stimulus(subj.stim.sig{ii});
        myTDT.play()
      else
        sound(subj.stim.sig{ii},fs)
      end
      pause(2)
    end
    DrawFormattedText(win,...
      ' Stimuli precued by a tone require your response. \n Press any key to hear examples!',...
      'center','center',white);
    Screen('Flip',win);
    KbStrokeWait;
    Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2); % Fixation point
    Screen('Flip',win);
    pause(1)
    for ii = [1,5]
      if flags.do_TDTon
        myTDT.load_stimulus(precursor);
        myTDT.play()
        pause(0.5)
        myTDT.load_stimulus(subj.stim.sig{ii});
        myTDT.play()
      else
        sound(precursor,fs)
        pause(0.5)
        sound(subj.stim.sig{ii},fs)
      end
      pause(2)
    end
    DrawFormattedText(win,'Press any key to start the experiment!','center','center',white);
    Screen('Flip',win);
end

%% Test procedure
KbStrokeWait;

% define all possible stimulus combinations
NM = length(kv.M);
if flags.do_noStaticSounds
  Ndir = 2; % only L/R
else
  Ndir = 3; % also static
end
id_M = repmat(transpose(1:NM),[1,Ndir]);
id_M = repmat(id_M(:),[1,kv.Nrep]);
id_dirChange = repmat(1:Ndir,[NM,1]);
id_dirChange = repmat(id_dirChange(:),[1,kv.Nrep]);
iRT = false(NM*Ndir,kv.Nrep);
iRT(:,1:round(kv.Nrep*kv.responseTrialRatio)) = true;
if flags.do_staticSounds
  iRT(2*NM+1:3*NM,:) = false; % no static response trials
end

% randomizations
iC = nan(NM*Ndir,kv.Nrep); % index for M-direction combinations
for irep = 1:kv.Nrep % randomize within blocks
  iC(:,irep) = randperm(NM*Ndir)+(irep-1)*NM*Ndir;
end
for ii = 1:NM*Ndir % randomize across blocks
  iRT(ii,:) = iRT(ii,randperm(kv.Nrep));
end
iRT = iRT(iC);
id_M = id_M(iC);
id_dirChange = id_dirChange(iC);

% initialize variables
subj.M = kv.M(id_M); % M combination
subj.dirChange = id_dirChange; % M combination
subj.SPL = nan(NM*Ndir,kv.Nrep); % SPL
subj.resp = nan(NM*Ndir,kv.Nrep); % L/R response: 1...left, -1...right
subj.RT = nan(NM*Ndir,kv.Nrep); % reaction time
subj.hit = nan(NM*Ndir,kv.Nrep); % reaction time

for irep = 1:kv.Nrep
  
  if mod(irep,2)==1 % is odd
    if flags.do_TDTon
      send_event(myTDT, trigVals.unpause) % per my .cfg file for the Biosemi, a trigger of 254 unpauses and 255 pauses the EEG recording
    end
    pause(1)  % add a pause to ensure the pause is active before the sound 
  end
  
  % Fixation point
  Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
  Screen('Flip',win);
  pause(1)

  % Randomized presentation order
  for ii = 1:NM*Ndir
%     idx = iC(ii,irep);
    sigpair = subj.stim.sig{id_M(ii,irep),id_dirChange(ii,irep)};
    changeTime = 0.7;
    nM2 = changeTime*fs;
    
%     sig1 = subj.stim.sig{id_M(idx),1};
%     sig2 = subj.stimLR.sig{id_M(idx),id_dirChange(idx)};
% 
%     % combine stimulus pairs with temporal jitter of crossfade
%     dt = kv.jitter*(rand-0.5);
%     changeTime = kv.changeTime+dt;
%     [sigpair,nM2] = SpExCue_crossfade(sig1,sig2,...
%       subj.stim.fs,kv.signalDuration,changeTime,kv.fadeDur,1);
    
    % level roving
    dSPL = kv.SPLrove*(rand-0.5);
    subj.SPL(ii,irep) = kv.SPL + dSPL;
    sigpair = 10^(dSPL/20)*sigpair;
    
    % precursor
    if iRT(ii,irep)
%       if  id_dirChange(idx)==3 % skip static response trials
%         subj.RT(ii,irep) = false;
%         continue
%       end
      if flags.do_TDTon
        myTDT.load_stimulus(precursor, [1,trigVals.precursorOnset]);
        myTDT.play()
      else
        sound(precursor,fs)
      end
      pause(0.5)
    else
      pause(0.6)
    end
    
    % playback triggers
    TrigValM = 4*subj.M(ii,irep)+1;
    TrigValDir = 100 + 10*subj.dirChange(ii,irep);
    TrigValRT = iRT(ii,irep)*trigVals.responseTrial;
    triggerInfo = [   1,TrigValM; ...
                    nM2,TrigValM+TrigValDir+TrigValRT;...
                    size(sigpair,1),trigVals.stimulusOffset];

    % prepare for response via keyboard 
    keyisdown = 1;
    while(keyisdown) % first wait until all keys are released
      [keyisdown,secs,keycode] = KbCheck;
      WaitSecs(0.001); % delay to prevent CPU hogging 
    end
    
    if ii>1; disp(toc); end; tic
    
    % playback
    if flags.do_TDTon
      myTDT.load_stimulus(sigpair, triggerInfo);
      myTDT.play()
    else
      sound(sigpair,fs)
    end
    
    if iRT(ii,irep)
      % record reaction time
  %     tic;
      start_time = GetSecs + changeTime;
  %     wait_time = rand * (tmax-tmin) + tmin;
      while not(keycode(rkey) || keycode(lkey))
  %       if GetSecs-start_time > wait_time
  %         wait_time = Inf; % so as not to repeat this part Screen('FillOval',w,[255 255 255],[x0-r,y0-r,x0+r,y0+r]); Screen('Flip',w);
  %         time0=GetSecs;
  %       end
        [keyisdown,secs,keycode] = KbCheck;
        WaitSecs(0.001); % delay to prevent CPU hogging 
      end
      subj.RT(ii,irep) = secs - start_time;
  %     subj.RT2(ii,irep) = toc - changeTime; % for debugging

      % response evaluation
      if keycode(lkey)
        subj.resp(ii,irep) = -1;
      elseif keycode(rkey)
        subj.resp(ii,irep) = 1;
      end
      subj.hit(ii,irep) = subj.resp(ii,irep) == sign(subj.dirChange(ii,irep)-1.5);

%       if flags.do_TDTon
%         if subj.hit(ii)
%           send_event(myTDT, trigVals.correctResponse)
%         else
%           send_event(myTDT, trigVals.wrongResponse)
%         end
%       end

    % feedback
%     if flags.do_feedback
%       if subj.hit(ii)
%         Screen('DrawDots',win, [x_center,y_center], 14, green, [], 2);
%       else
%         Screen('DrawDots',win, [x_center,y_center], 14, red, [], 2);
%       end
%       Screen('Flip',win);
%       pause(0.4)
%       Screen('DrawDots',win, [x_center,y_center], 14, white, [], 2);
%       Screen('Flip',win);
%     else
%       pause(0.4)
%     end
      pause(max(1.2-subj.RT(ii,irep)-changeTime,0))
    else
      pause(1.2)
    end
    pause(0.2 + kv.jitter*(rand-0.5))

  end
  
  if mod(irep,2)==0 || irep==kv.Nrep % if is even - allow break every second round
    % Pause TDT
    if flags.do_TDTon
      send_event(myTDT, trigVals.pause)  % per my .cfg file for the Biosemi, a trigger of 254 unpauses and 255 pauses the EEG recording
      myTDT.reset();  % rewind and clear buffers
    end

    % Display time course and intermediate score
    infotext = [num2str(irep) ' of ' num2str(kv.Nrep) ' blocks completed.'];
    if flags.do_feedback
      if sum(iRT(1:ii*irep)) > 0
        pcorrect = 100*sum(subj.hit(iRT(1:ii*irep)))/sum(iRT(1:ii*irep));
        infotext = [infotext,'\n\n\n',num2str(pcorrect,'%3.2f'),'% correct.'];
      end
      infotext = [infotext,'\n\n\n','Press a key to continue!'];
    end
    DrawFormattedText(win,infotext,'center','center',white);
    Screen('Flip',win);
    
    % Save results
    save(fn,'subj','kv','flags')
    
    KbStrokeWait;
  end
  
end

%% Inform listener that experiment is completed
i1 = 'Experiment completed.';
i2 = '\n\n\n\n Press any key to close this window';

DrawFormattedText(win,[i1, i2],'center','center',white);
Screen('Flip',win);
KbStrokeWait;
Screen('CloseAll');

end