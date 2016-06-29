function analyzeEEGpilot3_eeglab_ICbasedBlinkRej(subjID)

if not(exist('subjID','var'))
  subjID = input('Subject ID: ','s');
end

EEGfnRaw = ['SpExCue_EEGpilot3_',subjID,'_raw.set'];
EEGfnClean = ['SpExCue_EEGpilot3_',subjID,'_blinkICrej'];
ICAfn = ['paramsOut_analyzeEEGpilot3_eeglab_preprocessICA_',subjID];
filepath = '../data/';

if not(exist('pop_loadset','file'))
  addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab/')
  eeglab
end

%% Load data
EEG = pop_loadset('filename', EEGfnRaw, 'filepath', filepath); 
ICA = load(fullfile(filepath,ICAfn));

%% Transfer ICA data
icafields = {'icaweights','icasphere','icawinv','icachansind'};
for ii = 1:length(icafields)
  eval(['EEG.',icafields{ii},'=ICA.paramsOut.',icafields{ii},';'])
end

%% Filtering
fcutoffs = [.5,100];  % frequency cut-offs in Hz
transitionBandwidth = 1;  % Hz
maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
KaiserWindowBeta = pop_kaiserbeta(maxPassbandRipple);
filtOrder = pop_firwsord('kaiser', EEG.srate, transitionBandwidth, maxPassbandRipple);
EEG = pop_firws(EEG, 'fcutoff', fcutoffs, 'ftype', 'bandpass', 'wtype', 'kaiser', 'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);
EEGfnClean = [EEGfnClean,'_fmax',num2str(fcutoffs(2))];

%% Select components to remove
pop_eegplot( EEG, 0, 1, 1);
pop_selectcomps(EEG);

%% Remove IC representing eye blinks
EEG = pop_subcomp( EEG );
ICA.paramsOut.indICrej = input('Removed IC:');

%% Interpolate bad channels
if not(isempty(ICA.paramsOut.indelecReject))
  EEG = pop_interp(EEG,ICA.paramsOut.indelecReject);
end

%% Save 
EEG = eeg_checkset(EEG);
EEG = pop_saveset(EEG,  'filename', fullfile(filepath,EEGfnClean));
paramsOut = ICA.paramsOut;
save(fullfile(filepath,ICAfn),'paramsOut')

close all

end