function SpExCue_analyzeExp1eeg_dipfit(subjID)
% fit (bilateral) dipole locations to every IC

%% Input and paths
if not(exist('subjID','var'))
  subjID = input('Subject ID: ','s');
end

EEGfn = ['Exp1eeg_',subjID,'_ICAraw.set'];
% EEGfn = ['Exp1eeg_',subjID,'_blinkICrej.set'];
filepath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/data/';

%% Start EEGLAB
if not(exist('pop_loadset','file'))
  addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab/')
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%% Load data
EEG = pop_loadset('filename', EEGfn, 'filepath', filepath); 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw

%% Fit dipoles
nIC = size(EEG.icaact,1); % number of ICs

% EEG = pop_dipfit_settings( EEG, 'hdmfile','/Users/rbaumgartner/Documents/MATLAB/eeglab/plugins/dipfit2.3/standard_BESA/standard_BESA.mat','coordformat','Spherical','mrifile','/Users/rbaumgartner/Documents/MATLAB/eeglab/plugins/dipfit2.3/standard_BESA/avg152t1.mat','chanfile','/Users/rbaumgartner/Documents/MATLAB/eeglab/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp','chansel',[] );
EEG = pop_dipfit_settings( EEG, 'hdmfile','/Users/rbaumgartner/Documents/MATLAB/eeglab/plugins/dipfit2.3/standard_BEM/standard_vol.mat','coordformat','MNI','mrifile','/Users/rbaumgartner/Documents/MATLAB/eeglab/plugins/dipfit2.3/standard_BEM/standard_mri.mat','chanfile','/Users/rbaumgartner/Documents/MATLAB/eeglab/plugins/dipfit2.3/standard_BEM/elec/standard_1005.elc','coord_transform',[0.77975 -15.6368 2.9658 0.084126 0.0017284 -1.5735 99.9634 89.9092 97.0478] ,'chansel',[] );
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_dipfit_gridsearch(EEG, [1:nIC] ,[-85:17:85] ,[-85:17:85] ,[0:17:85] ,0.4);
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_multifit(EEG, 1:nIC,'threshold',100,'dipoles',1,'dipplot','off','plotopt',{'normlen' 'on'});
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = fitTwoDipoles(EEG, 'LRR', 35);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Save data
EEG = eeg_checkset(EEG);
EEG = pop_saveset(EEG,  'filename', fullfile(filepath,EEGfn));
end