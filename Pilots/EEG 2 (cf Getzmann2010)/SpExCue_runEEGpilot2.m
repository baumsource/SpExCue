% SpExCue_runEEGpilot

if not(exist('fftreal','file'))
    addpath(fullfile('..','..','MATLAB_general'))
    addpath(fullfile('..','..','MATLAB_general','sofa_new','API_MO'))
    addpath(fullfile('..','..','MATLAB_general','ltfat'))
    ltfatstart
    SOFAstart
end

SpExCue_EEGpilot2('RB','M',[1,1/3,0],'responseTrialRatio',.6,'skipFamiliarization',...
  'screenNumber',0,'TDToff','debugMode')