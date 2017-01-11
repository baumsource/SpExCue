% SpExCue_runExp1

% Check amp setting: -12 dB !!!

%% Listener-specific settings
ID = 'RB'; % RS2
procedure = {... 
%   'screening';...
  'behavioral';...
%   'eeg';...
  }; 
aziBehav = 30; % depends on screening result (only directions where distance is dominant movement cue)
aziEEG = 90; % depends on behavioral result (direction of max dprime)

%% General settings
M = [0,0.5,1];
HRTFs = 'HRCeq';
Feedback = 'noFeedback';
roving = 'noRoving';
screenNumber = 0;

%% Load dependencies
addpath(fullfile('..','MATLAB_general'))
addpath(fullfile('..','..','sofa','API_MO'))
addpath(fullfile('..','..','ltfat'))

ltfatstart
SOFAstart
sca

%% Start test procedure
switch procedure{1}
  case 'screening'
        
%     SpExCue_Exp1_screening(ID,'azi',[-90,0,90],HRTFs,'screenNumber',screenNumber)
    
  case 'behavioral'
    %% Behavioral pilot (3 repetitions for 3 positions take 5 min) 
    Nrep = 2*6; % 252 trials -> 13 min presentation time -> 20 min with breaks
    fnExtension = 'behav';
    
    SpExCue_Exp2(ID,'M',M,'azi',aziBehav,'Nrep',Nrep,...
      Feedback,roving,'fnExtension',fnExtension,HRTFs,'screenNumber',screenNumber,'changeM')
  
    cd analysis
    SpExCue_analyzeExp2BTL(ID,fnExtension)
%     SpExCue_analyzeExp2behav(ID,fnExtension);
    cd ..
    
  case 'eeg'
    %% EEG monitored experiment 
    Nrep = 120; % 840 trials -> 42 min presentation time -> 60 min with breaks
    fnExtension = 'eeg';
    SpExCue_Exp2(ID,'M',M,'azi',aziEEG,'Nrep',Nrep,...
      Feedback,roving,'fnExtension',fnExtension,HRTFs,'screenNumber',screenNumber)
  
    cd analysis
    SpExCue_analyzeExp1behav(ID,fnExtension);
    cd ..
    
end

% Stimulus level calibration
% SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',0,'Nrep',12,'fnExtension','debug','debugMode','SPL',70,'dur',30)