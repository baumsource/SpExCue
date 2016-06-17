% SpExCue_runEEGpilot

ID = 'RB'; % RS2

% addpath(fullfile('..','..','MATLAB_general'))
% addpath(fullfile('..','..','MATLAB_general','sofa_new','API_MO'))
% addpath(fullfile('..','..','MATLAB_general','ltfat'))

% ltfatstart
% SOFAstart

% behavioral pilot
SpExCue_EEGpilot3(ID,'screenNumber',0,'azi',-90:90:90,'Nrep',12,...
  'fnExtension','checkDistortion','TDToff','debugMode','skipFamiliarization')

% EEG
% SpExCue_EEGpilot(ID,'screenNumber',1,'azi',-90)