% General processing pipeline for Exp1

epochThresh = 70; % �V

tmp = load('SpExCue_Exp1eeg_subjects.mat');
subjects = tmp.subjects;

%% ICA
for ii=1:length(subjects)
  SpExCue_analyzeExp1eeg_ICA(subjects{ii});
end

%% IC rejection to remove eye blinks and horizontal eye movements
for ii=1:length(subjects)
  SpExCue_analyzeExp1eeg_eyeBlinkICrejection(subjects{ii});
end

%% Data split according to physical change
% all ICs
SpExCue_analyzeExp1eeg_dataSplit('IDs',subjects,'noGrouping',...
  'filenameIDx','Exp1eeg_IDx_ICAraw.set');
% artificial ICs removed
SpExCue_analyzeExp1eeg_dataSplit('IDs',subjects,'noGrouping',...
  'filenameIDx','Exp1eeg_IDx_blinkICrej.set','epochThresh',epochThresh);

%% Data split according to behavioral response
% all ICs
SpExCue_analyzeExp1eeg_dataSplit('IDs',subjects,'groupE',...
  'filenameIDx','Exp1eeg_IDx_ICAraw.set');
% artificial ICs removed
SpExCue_analyzeExp1eeg_dataSplit('IDs',subjects,'groupE',...
  'filenameIDx','Exp1eeg_IDx_blinkICrej.set','epochThresh',epochThresh);

%% Statistical analaysis
% 