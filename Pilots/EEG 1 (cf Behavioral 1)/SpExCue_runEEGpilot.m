% SpExCue_runEEGpilot

ID = 'RB'; % RS2

% addpath(fullfile('..','..','MATLAB_general'))
% addpath(fullfile('..','..','MATLAB_general','sofa_new','API_MO'))
% addpath(fullfile('..','..','MATLAB_general','ltfat'))

ltfatstart
SOFAstart

% behavioral pilot
SpExCue_EEGpilot(ID,'screenNumber',0,'azi',90,'Nrep',6,...
  'fnExtension','debug1','TDToff','TbTfeedback',...
  'rdepth',10,'flow',1e3,'fhigh',16e3,'debugMode')

% EEG
% SpExCue_EEGpilot(ID,'screenNumber',1,'azi',-90)