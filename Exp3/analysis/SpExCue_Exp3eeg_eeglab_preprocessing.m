% General processing pipeline for SpExCue Exp. 3

% Primary path definitions required for this analysis
eeglab_analysis_scripts = fullfile('..','..','Matlab_general','eeglab_analysis');
dirBDF = fullfile('..','data'); % BDF data directory
experimentString = 'Exp3eeg';  % String to label experiment number for the saved eeglab dataset files
cnfg.BDFfn = fullfile(dirBDF,[experimentString,'_Sxx']);

cnfg.dirLOCS = dirBDF;
cnfg.nChansLocsFilename = 'chanLocs_refA.locs';  % filename of the eeglab locations file

cnfg.refChanNum = [33 34]; % Indicies of reference channels in BDF file -> apply changes also to s
cnfg.refChanLabel = 'A1 A2';%'TP9 TP10'; % Labels of reference channels
cnfg.eyeChanNum = 35:37;
cnfg.eegChanNum = 1:32;

cnfg.FS = 500;  % In samps/sec
cnfg.locutoff = 0.5;  % In Hz
cnfg.hicutoff = 20;  % In Hz
cnfg.transitionBandwidth = 1;  % In Hz
cnfg.maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
cnfg.baselineRange = [-100,0];%[-500,0]; % baseline time range in ms
cnfg.globalThreshold = [-200 800]; % µV
cnfg.finalThreshold = 50; % post IC rejection
cnfg.stimulusEpochRange = [-.5,4];
cnfg.finalEpochRange = [-0.6,1];%[-1.0,2.5];

cnfg.timeLockedEventNum = {250}; % cue trigger

% time-locked to first syllable
% cnfg.cond =  {cnfg.timeLockedEventNum(1:4),'ITD_right'; ...
%               cnfg.timeLockedEventNum(5:8),'ITD_left'; ...
%               cnfg.timeLockedEventNum(9:12),'ILD_right'; ...
%               cnfg.timeLockedEventNum(13:16),'ILD_left'; ...
%               cnfg.timeLockedEventNum(17:20),'HRTF_right'; ...
%               cnfg.timeLockedEventNum(21:24),'HRTF_left'};
            
% time-locked to second (embedded) syllable (NOTE: use/see
% SpExCue_Exp3_trigCode for help/info)
attend = num2cell(17:4:61);
unattend = num2cell(81:4:125);
cnfg.cond =  {attend(1:4),'ITD'; attend(5:8),'ILD'; attend(9:12),'HRTF';...
              unattend(1:4),'ITDu'; unattend(5:8),'ILDu'; unattend(9:12),'HRTFu'};
save([experimentString,'_cnfg'],'cnfg')
            
tmp = load('SpExCue_Exp3eeg_subjects.mat');
subjects = tmp.subject;

% Select subjects as required
subjects = subjects(5,:); disp(['only ',subjects.name])

%% Add paths
if not(exist('eeglab','file'))
  error('RB: eeglab not available.')
end
addpath(eeglab_analysis_scripts)

%% pre-ICA processing
for ii=1:height(subjects)
  fn_bdf = strrep(cnfg.BDFfn,'Sxx',subjects.name{ii});
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