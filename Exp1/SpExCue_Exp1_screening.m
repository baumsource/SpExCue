function tab = SpExCue_Exp1_screening( ID,varargin )
%SpExCue_Exp1_screening - Screening to assess affected dimensionality by changes
%in spectral contrast
%
% Key-value pairs:
% keyvals.Mcomb = [1,0];
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

%% Check availability of dependent functions
if not(exist('SpExCue_stim','file'))
  addpath(fullfile('..','MATLAB_general'))
end
if not(exist('fftreal','file'))
  amtstart
end

saveFileNamePrefix = 'Exp1behav';
saveFileNamePostfix = 'screen';

%% Experimental variables
definput = arg_SpExcue;
definput.keyvals.Mcomb = [1,0];
[flags,kv]  = ltfatarghelper({},definput,varargin);

if length(kv.azi) > 1 && length(kv.ele) == 1
  kv.ele = kv.ele*ones(length(kv.azi),1);
end

%% Enter listener ID
if not(exist('ID','var'))
  subj.ID = upper(input('Enter subject ID: ','s'));
else
  subj.ID = ID;
end

%% Initialize the TDT and playback system
fsVal = 48; % sampling rate, 48 stands for 48828.125 Hz (-> max buffer length 85 seconds)
minmaxChVolt = db2mag(3*3)*2.53;  % min/max voltage for the output channels (scaling), 2.53 yields 78db SPL at -3dB headphone amplification
trigDuration = 0.005;  % duration of trigger in seconds

if flags.do_TDTon
  myTDT = tdt('playback_2channel', fsVal, minmaxChVolt, trigDuration );
end

%% Psychtoolbox seetings and initialization

%     Screen('Preference', 'SkipSyncTests', 1);

KbName('UnifyKeyNames');
key.none = [KbName('0'),KbName('0)')];
key.elevation = [KbName('1'),KbName('1!')];
key.distance = [KbName('2'),KbName('2@')];
key.yes = KbName('y');
key.no = KbName('n');
key.play = KbName('space');

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

Screen('WindowSize', win);
[x_center,y_center] = RectCenter(winRect);


%% Stimulus generation
pos = [kv.azi(:),kv.ele(:)];
Npos = length(kv.azi);
% pos = pos(randperm(Npos),:);
if flags.do_TDTon
  fs = myTDT.sampleRate;
else
  fs = 48.8e3;
end

subj.stim = SpExCue_stim( kv.Mcomb,subj.ID,pos,round(fs),'argimport',flags,kv );

for pp = 1:Npos
  sigpair{pp} = SpExCue_crossfade(subj.stim.sig{1,pp},subj.stim.sig{2,pp},...
            subj.stim.fs,kv.dur,kv.dur/2,kv.xFadeDur);
end

%% Movement?
% Instruction
infotext = [...
  'You will hear three different sounds. The first one is coming from the right,\n\n',...
  'the second from the front, and the third from the left.\n\n',...
  '   Press *spacebar* to play each sound as often as you want and answer the following question! \n\n \n\n \n\n',...
  'Does the sound appear to move? If yes, does it move mainly in elevation or distance?\n\n',...
  '(Elevation and distance are considered relative to the center of your head.)\n\n',...
  '   Press *0* for no movement or change in location!\n\n', ...
  '   Press *1* for movement in elevation (upwards/downwards)! \n\n',...
  '   Press *2* for movement in distance (closer/farther)! \n\n \n\n \n\n',...
  'This screen remains displayed during stimulus presentation.'];
DrawFormattedText(win,infotext,.2*x_center,'center',white);
Screen('Flip',win);

movement = cell(Npos,1);
for pp = 1:Npos
  keyCodeVal = 0;
  played = false;
  while not(any(keyCodeVal==[key.none,key.distance,key.elevation]) && played)
    [tmp,keyCode] = KbWait([],2);
    keyCodeVal = find(keyCode,1);
    disp(num2str(keyCodeVal))
    if keyCodeVal == key.play
      pause(.1)
      if flags.do_TDTon
        myTDT.load_stimulus(sigpair{pp});
        myTDT.play_blocking()
      else
        sound(sigpair{pp},fs)
      end
      played = true;
    end
  end
  
  if any(keyCodeVal==key.elevation)
    movement{pp} = 'elevation';
  elseif any(keyCodeVal==key.distance)
    movement{pp} = 'distance';
  elseif keyCodeVal==key.none
    movement{pp} = 'none';
  else
    movement{pp} = 'NA';
  end
  
  if pp < Npos
    DrawFormattedText(win,'Next stimulus...',.2*x_center,'center',white);
    Screen('Flip',win);
    pause(1)
    DrawFormattedText(win,infotext,.2*x_center,'center',white);
    Screen('Flip',win);
  else
    DrawFormattedText(win,'Thank you!',.2*x_center,'center',white);
    Screen('Flip',win);
    pause(1)
  end
  
end


%% Ceck results and go on with distance if distance was reported for one or two directions
distance = cell(Npos,1);

distfalg = ismember(movement,'distance');
if any(distfalg) && not(all(distfalg))

  % Instruction 2: distance?
    infotext2 = [...
      'Now, listen to the same sounds again, but only focus on changes in distance and ignore other changes!\n\n',...
      'Press *spacebar* to play the sound as often as you want and answer the following question!\n\n\n\n\n\n',...
      'Does the sound change in distance?\n\n',...
      '   Press *y* for YES!\n\n', ...
      '   Press *n* for NO!'];
    DrawFormattedText(win,infotext2,.2*x_center,'center',white);
    Screen('Flip',win);
    
  % play sounds
  for pp = 1:Npos
      keyCodeVal = 0;
      played = false;
      while not(any(keyCodeVal==[key.yes,key.no]) && played)
        [tmp,keyCode] = KbWait([],2);
        keyCodeVal = find(keyCode,1);
        disp(num2str(keyCodeVal))
        if keyCodeVal == key.play
          pause(.1)
          if flags.do_TDTon
            myTDT.load_stimulus(sigpair{pp});
            myTDT.play_blocking()
          else
            sound(sigpair{pp},fs)
          end
          played = true;
        end
      end

      if any(keyCodeVal==key.yes)
        distance{pp} = 'yes';
      elseif any(keyCodeVal==key.no)
        distance{pp} = 'no';
      else
        distance{pp} = 'NA';
      end

      if pp < Npos
        DrawFormattedText(win,'Next stimulus...',.2*x_center,'center',white);
        Screen('Flip',win);
        pause(1)
        DrawFormattedText(win,infotext2,.2*x_center,'center',white);
        Screen('Flip',win);
      else
        DrawFormattedText(win,'Thank you!',.2*x_center,'center',white);
        Screen('Flip',win);
        pause(1)
      end

  end
  
end
%% Save results
azimuth = kv.azi(:);
tab = table(azimuth,movement,distance);

savename = fullfile('data',[saveFileNamePrefix,'_',subj.ID,'_',saveFileNamePostfix]);
if not(exist('./data','dir'))
  mkdir('data')
end

% introduce counter in filename for multiple runs of screening
% tmp = dir([savename,'*']);
% savename = [savename,num2str(length(tmp)+1)];

save(savename,'tab')

sca
end

