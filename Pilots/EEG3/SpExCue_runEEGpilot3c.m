% SpExCue_runEEGpilot

% Amp setting: -12 dB !!!

ID = 'RB'; % RS2
procedure = {...
  'familiarization';...
%   'behavioral';...
%   'eeg';...
  }; % type 1 or true to run with EEG setting

addpath(fullfile('..','..','MATLAB_general'))
addpath(fullfile('..','..','..','sofa','API_MO'))
addpath(fullfile('..','..','..','ltfat'))

ltfatstart
SOFAstart
sca

% rippleDepth = 0;

switch procedure{1}
  case 'familiarization'
    %% Familiarization 
    Nrep = 6;
    fnExtension = 'debug';
    SpExCue_EEGpilot3(ID,'M',[1,0],'azi',-90:90:90,'Nrep',Nrep,...
      'noFeedback','changeM','skipFamiliarization','noRoving',...
      'screenNumber',1,'fnExtension',fnExtension)
    
    cd analysis
    SpExCue_analyzeEEGpilot3_behavior([ID,fnExtension])
    cd ..
    
  case 'behavioral'
    %% Behavioral pilot (3 repetitions for 3 positions take 5 min) 
    Nrep = 2*6; % 252 trials -> 13 min presentation time -> 20 min with breaks
    fnExtension = 'behav';
    
    SpExCue_EEGpilot3(ID,'azi',-90:90:90,'Nrep',Nrep,...
      'D0detectionFeedback','repeateM','skipFamiliarization','noRoving',...
      'screenNumber',1,'fnExtension',fnExtension)
  
    cd analysis
    SpExCue_analyzeEEGpilot3_behavior([ID,fnExtension])
    cd ..
    
  case 'eeg'
    %% EEG monitored experiment 
    Nrep = 120; % 840 trials -> 42 min presentation time -> 60 min with breaks
    fnExtension = 'eeg';
    SpExCue_EEGpilot3(ID,'azi',-90,'Nrep',Nrep,...
      'D0detectionFeedback','repeateM','skipFamiliarization','noRoving',...
      'screenNumber',1,'fnExtension',fnExtension)
  
    cd analysis
    SpExCue_analyzeEEGpilot3_behavior([ID,fnExtension])
    cd ..
    
end
% Stimulus level calibration
% SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',0,'Nrep',12,'fnExtension','debug','debugMode','SPL',70,'dur',30)