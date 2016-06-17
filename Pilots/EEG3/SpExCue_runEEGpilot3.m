% SpExCue_runEEGpilot

ID = 'RB'; % RS2
eegflag = 0; % type 1 or true to run with EEG setting

addpath(fullfile('..','..','MATLAB_general'))
addpath(fullfile('..','..','..','sofa','API_MO'))
addpath(fullfile('..','..','..','ltfat'))

ltfatstart
SOFAstart
sca

if not(eegflag)
    %% behavioral pilot (3 repetitions for 3 positions take 5 min) 
    Nrep = 2*6; % 252 trials -> 13 min presentation time -> 20 min with breaks
    SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',-90:90:90,'Nrep',Nrep,...
      'fnExtension','behav','noFeedback','halfPosOrderPermutation')
    SpExCue_analyzeEEGpilot3_behavior([ID,'behav'])
else
    %% EEG monitored experiment 
    Nrep = 120; % 840 trials -> 42 min presentation time -> 60 min with breaks
    SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',90,'Nrep',120,...
      'fnExtension','eeg','blockedFeedback')
    SpExCue_analyzeEEGpilot3_behavior([ID,'eeg'])
end
% Stimulus level calibration
% SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',0,'Nrep',12,'fnExtension','debug','debugMode','SPL',70,'dur',30)