function fnsave = eeglab_postICA(cnfg,fnICA)
% eeglab_postICA - rejection of artifact ICs, final thresholding,
% interpolation of bad channels
% 
% PROBLEM: detail plots in pop_selectcomps do not work because 'eeglab
% redraw' malfunctions

% if not(exist('subjID','var'))
%   subjID = input('Subject ID: ','s');
% end
% 
% EEGfnRaw = ['Exp1eeg_',subjID,'_ICAraw.set'];
% EEGfnClean = ['Exp1eeg_',subjID,'_blinkICrej'];
% ICAfn = ['paramsOut_SpExCue_analyzeExp1eeg_ICA_',subjID];
% filepath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/data/';
% 
% if not(exist('pop_loadset','file'))
%   addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab/')
% end
% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%% Load data
[pathstr,name] = fileparts(fnICA);
EEG = pop_loadset('filename', [name,'.set'], 'filepath', pathstr); 
[ALLEEG EEG CURRENTSET] = eeg_store([], EEG);
eeglab redraw % DOES NOT WORK!
% ICA = load(fullfile(filepath,ICAfn));
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
pop_selectcomps(EEG,1:18); % does not work properly because EEG.icaact is empty (but I don't know why it is?!?)

%% Remove IC representing eye blinks
ICA.paramsOut.indICrej = 0;
while ICA.paramsOut.indICrej == 0
  EEG = pop_subcomp( EEG );
  ICA.paramsOut.indICrej = input('Removed IC:');
end

%% Final threshold
[EEG,indtrialRej] = pop_eegthresh(EEG,1,cnfg.eegChanNum,...
  -cnfg.finalThreshold,+cnfg.finalThreshold,...
  cnfg.stimulusEpochRange(1),cnfg.stimulusEpochRange(2),0,1);

%% Interpolate bad channels
fnBadChan = strrep(fnICA,'ICA','indelecReject.mat');
if exist(fnBadChan,'file') % not(isempty(ICA.paramsOut.indelecReject))
  badChan = load(fnBadChan);
  EEG = pop_interp(EEG,badChan.indelecReject);
end

%% Save 
EEG = eeg_checkset(EEG);
if strfind(name,'ICAraw')
  EEGfnClean = strrep(name,'ICAraw','clean');
else
  EEGfnClean = [name,'_clean'];
end
fnsave = fullfile(pathstr,[EEGfnClean,'.set']);
EEG = pop_saveset(EEG,  'filename', fnsave);
paramsOut = ICA.paramsOut;
save(fullfile(filepath,ICAfn),'paramsOut')

% close all

end