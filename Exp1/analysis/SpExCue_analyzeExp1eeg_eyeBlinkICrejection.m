function SpExCue_analyzeExp1eeg_eyeBlinkICrejection(subjID)
% PROBLEM: detail plots in pop_selectcomps do not work

if not(exist('subjID','var'))
  subjID = input('Subject ID: ','s');
end

EEGfnRaw = ['Exp1eeg_',subjID,'_ICAraw.set'];
EEGfnClean = ['Exp1eeg_',subjID,'_blinkICrej'];
ICAfn = ['paramsOut_SpExCue_analyzeExp1eeg_ICA_',subjID];
filepath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/data/';

if not(exist('pop_loadset','file'))
  addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab/')
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%% Load data
EEG = pop_loadset('filename', EEGfnRaw, 'filepath', filepath); 
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw
ICA = load(fullfile(filepath,ICAfn));
% 
% %% Transfer ICA data
% icafields = {'icaweights','icasphere','icawinv','icachansind'};
% for ii = 1:length(icafields)
%   eval(['EEG.',icafields{ii},'=ICA.paramsOut.',icafields{ii},';'])
% end
% 
% %% Filtering
% fcutoffs = [.5,100];  % frequency cut-offs in Hz
% transitionBandwidth = 1;  % Hz
% maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
% KaiserWindowBeta = pop_kaiserbeta(maxPassbandRipple);
% filtOrder = pop_firwsord('kaiser', EEG.srate, transitionBandwidth, maxPassbandRipple);
% EEG = pop_firws(EEG, 'fcutoff', fcutoffs, 'ftype', 'bandpass', 'wtype', 'kaiser', 'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);
% EEGfnClean = [EEGfnClean,'_fmax',num2str(fcutoffs(2))];

%% Select components to remove
pop_eegplot( EEG, 0, 1, 1, [], 'winlength',10,'events','off');
pop_selectcomps(EEG,1:18);

%% Remove IC representing eye blinks
ICA.paramsOut.indICrej = 0;
while ICA.paramsOut.indICrej == 0
  EEG = pop_subcomp( EEG );
  ICA.paramsOut.indICrej = input('Removed IC:');
end

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