function response = SpExCue_Screening( ID,varargin )
%SpExCue_Screening - Screening to assess affected dimensionality by changes
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

saveFileNamePrefix = 'EEG3pilotScreening';

%% Experimental variables
definput = arg_SpExcue;
definput.keyvals.Mcomb = [1,0];
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

%% Initialize the TDT and playback system
fsVal = 48; % sampling rate, 48 stands for 48828.125 Hz (-> max buffer length 85 seconds)
minmaxChVolt = db2mag(3*3)*2.53;  % min/max voltage for the output channels (scaling), 2.53 yields 78db SPL at -3dB headphone amplification
trigDuration = 0.005;  % duration of trigger in seconds

if flags.do_TDTon
  myTDT = tdt('playback_2channel', fsVal, minmaxChVolt, trigDuration );
end

%% Psychtoolbox seetings and initialization
if flags.do_TDToff
    Screen('Preference', 'SkipSyncTests', 1);
end

KbName('UnifyKeyNames');
key.none = KbName('0');
key.LR = [KbName('1'),KbName('1!')];
key.UD = [KbName('2'),KbName('2@')];
key.FB = [KbName('3'),KbName('3#')];
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
subj.stim = SpExCue_stim( kv.Mcomb,subj.ID,pos,round(fs),kv.flow,kv.fhigh,kv.SPL,flags.HRTFs,'signalDuration',kv.dur );

%% Instruction
infotext = ['Press *spacebar* to play the sound as often as you want and answer the following question(s)!\n\n\n',...
	'Does the sound appear to move? \n',...
  'If not, press *0*! \n', ...
  'If it does, what is the dominant direction of movement? \n',...
  'Press *1* if it is left/right! \n',...
  'Press *2* if it is up/down! \n',...
  'Press *3* if it is front/back!'];
DrawFormattedText(win,infotext,'center','center',white);
Screen('Flip',win);

%% Play sounds
response = cell(Npos,1);
for pp = 1:Npos
  sigpair = SpExCue_crossfade(subj.stim.sig{1,pp},subj.stim.sig{2,pp},...
            subj.stim.fs,kv.dur,kv.dur/2,kv.xFadeDur);
  keyCodeVal = 0;
  while not(any(keyCodeVal==[key.none,key.LR,key.UD,key.FB]))
    [tmp,keyCode] = KbWait([],2);
    keyCodeVal = find(keyCode,1);
    disp(num2str(keyCodeVal))
    if keyCodeVal == key.play
      if flags.do_TDTon
        myTDT.load_stimulus(sigpair);
        myTDT.play_blocking()
      else
        sound(sigpair,fs)
      end
    end
  end
  if any(keyCodeVal==key.LR)
    response{pp} = 'LeftRight';
  elseif any(keyCodeVal==key.UD)
    response{pp} = 'UpDown';
  elseif any(keyCodeVal==key.FB)
    response{pp} = 'FrontBack';
  elseif keyCodeVal==key.none
    response{pp} = 'none';
  else
    response{pp} = 'NA';
  end
end

sca
end

