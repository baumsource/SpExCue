% SpExCue_runEEGpilot

ID = 'RB'; % RS2
procedure = {...
  'familiarization';...
  'behavioral';...
  'eeg';...
  }; % type 1 or true to run with EEG setting

addpath(fullfile('..','..','MATLAB_general'))
addpath(fullfile('..','..','..','sofa','API_MO'))
addpath(fullfile('..','..','..','ltfat'))
addpath(fullfile('.','analysis'))

ltfatstart
SOFAstart
sca

switch procedure{1}
  case 'familiarization'
    %% Familiarization 
    Nrep = 6;
    SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',-90:90:90,'Nrep',Nrep,...
      'fnExtension','famili','noFeedback','changeM','skipFamiliarization')
    SpExCue_analyzeEEGpilot3_behavior([ID,'famili'])
    
  case 'behavioral'
    %% Behavioral pilot (3 repetitions for 3 positions take 5 min) 
    Nrep = 2*6; % 252 trials -> 13 min presentation time -> 20 min with breaks
    SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',-90:90:90,'Nrep',Nrep,...
      'fnExtension','behav','consistencyFeedback','changeM','skipFamiliarization')
    SpExCue_analyzeEEGpilot3_behavior([ID,'behav'])
    
  case 'eeg'
    %% EEG monitored experiment 
    Nrep = 120; % 840 trials -> 42 min presentation time -> 60 min with breaks
    SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',90,'Nrep',Nrep,...
      'fnExtension','eeg','consistencyFeedback','changeM','skipFamiliarization')
    SpExCue_analyzeEEGpilot3_behavior([ID,'eeg'])
    
end
% Stimulus level calibration
% SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',0,'Nrep',12,'fnExtension','debug','debugMode','SPL',70,'dur',30)