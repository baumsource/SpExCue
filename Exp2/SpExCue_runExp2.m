% SpExCue_runExp1

% Check amp setting: -12 dB !!!

%% check TDT
% myTDT = tdt('playback_2channel_16bit', 48, 7, 0.005 )

%% Listener-specific settings
ID = 'S07'; % RS2
procedure = {... 
%   'screening';...
%   'LR';...
  'distance';...
  }; 
azi = -90;

%% General settings
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
    
  case 'LR'
    %% L/R discrimination 
    Nrep = 2*6; % 252 trials -> 13 min presentation time -> 20 min with breaks
    M = [0,1];
    aziLR = [azi,0,-azi];
    
    SpExCue_Exp2(ID,'M',M,'azi',aziLR,'Nrep',Nrep,...
      Feedback,roving,procedure{1},HRTFs,'screenNumber',screenNumber,'changeM')
  
    cd analysis
    SpExCue_analyzeExp2BTL(ID,procedure{1})
%     SpExCue_analyzeExp2behav(ID,fnExtension);
    cd ..
    
  case 'distance'
    %% distance discrimination 
    Nrep = 10; % 840 trials -> 42 min presentation time -> 60 min with breaks
    M = [0,0.5,1];
    SpExCue_Exp2(ID,'M',M,'azi',azi,'Nrep',Nrep,'responseBox',...
      Feedback,roving,procedure{1},HRTFs,'screenNumber',screenNumber,'changeM')
  
    cd analysis
    SpExCue_analyzeExp2BTL(ID,procedure{1});
    cd ..
    
end

% Stimulus level calibration
% SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',0,'Nrep',12,'fnExtension','debug','debugMode','SPL',70,'dur',30)