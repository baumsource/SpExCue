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
    %% behavioral pilot
    SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',-90:90:90,'Nrep',1*18,...
      'fnExtension','behav','noFeedback')
    SpExCue_analyzeEEGpilot3_behavior([ID,'behav'])
else
    %% EEG monitored experiment
    SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',90,'Nrep',...
      'fnExtension','eeg','blockedFeedback')
    SpExCue_analyzeEEGpilot3_behavior([ID,'eeg'])
end
% Stimulus level calibration
% SpExCue_EEGpilot3(ID,'screenNumber',1,'azi',0,'Nrep',12,'fnExtension','debug','debugMode','SPL',70,'dur',30)