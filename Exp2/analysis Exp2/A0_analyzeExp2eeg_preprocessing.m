% General processing pipeline for SpExCue Exp. 2

% Primary path definitions required for this analysis
dirBDF = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp2/data/';%string with the filename
experimentString = 'Exp2eeg';  % String to label experiment number for the saved eeglab dataset files

cnfg.BDFfn = fullfile(dirBDF,[experimentString,'_Sxx']);

cnfg.dirLOCS = dirBDF;
cnfg.nChansLocsFilename = 'chanLocs_refA.locs';  % filename of the eeglab locations file

cnfg.refChanNum = [33 34]; % Indicies of reference channels in BDF file -> apply changes also to s
cnfg.refChanLabel = 'A1 A2';%'TP9 TP10'; % Labels of reference channels
cnfg.eyeChanNum = 35:37;
cnfg.eegChanNum = 1:32;

cnfg.FS = 500;  % In samps/sec
cnfg.locutoff = 0.5;  % In Hz
cnfg.hicutoff = 25;  % In Hz
cnfg.transitionBandwidth = 1;  % In Hz
cnfg.maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
cnfg.baselineRange = [-200,0]; % baseline time range in ms
cnfg.globalThreshold = [-200 800]; % µV
cnfg.finalThreshold = 70; % post IC rejection
cnfg.stimulusEpochRange = [-1,3];
cnfg.finalEpochRange = [-.200,.600];

cnfg.timeLockedEventNum = {11,22,33,211,222,233};

cnfg.cond =  {;...
   {110},'C0_center'; {120},'C1_center'; {130},'KEMAR_center'; ...
   {011},'C0_left'; {022},'C1_left'; {033},'KEMAR_left'; ...
   {211},'C0_right'; {222},'C1_right'; {233},'KEMAR_right'};

tmp = load('SpExCue_Exp2eeg_subjects.mat');
subjects = tmp.subject;

%% Add paths
if not(exist('eeglab','file'))
  addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab')
end
addpath('/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Matlab_general/eeglab_analysis')

%% pre-ICA processing
for ii=1:height(subjects)
  fn_bdf = strrep(cnfg.BDFfn,'Sxx',subjects(ii,:).name);
  fn_preICA{ii} = eeglab_preICA(cnfg,fn_bdf);
end

%% ICA
for ii=1:height(subjects)
  fn_ICA{ii} = eeglab_ICA(fn_preICA{ii});
end

%% IC rejection to remove eye blinks and horizontal eye movements
for ii=1:height(subjects)
  fn_clean{ii} = eeglab_postICA(cnfg,fn_ICA{ii});
end

%% Dipole fitting
% for ii=1:height(subjects)
%   eeglab_dipfit(subjects{ii});
% end

%% Data split according to physical change
for ii = 1:height(subjects)
  eeglab_dataSplit(cnfg,fn_clean{ii},subjects(ii,:).name)
end