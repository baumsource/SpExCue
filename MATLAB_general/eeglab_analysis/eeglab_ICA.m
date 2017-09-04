function fnsave = eeglab_ICA(fn)
% eeglab_ICA - Perform ICA
%
% Usage: fn = eeglab_ICA(fnRaw)
%
% fnRaw ... filename (incl. path) to 

% AUTHOR: Robert Baumgartner, 2017/02/27

[pathstr,name] = fileparts(fn);
EEG = pop_loadset('filename', [name,'.set'], 'filepath', pathstr);

%% Run the ICA analysis
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'extended',1,'interupt','on','icatype','binica');
if strfind(fn,'preICA')
  fnsave = strrep(fn,'preICA','ICAraw');
else
  fnsave = [fn,'_ICAraw'];
end
EEG = pop_saveset(EEG,  'filename', fnsave);%, 'filepath', bdfFilepath);


%% Close EEGlab
% close all

%% Output and save the data parameters 
% paramsOut.timestamp = clock;
% paramsOut.subjID = subjID;
% paramsOut.FS = FS;
% paramsOut.locutoff = locutoff;
% paramsOut.transitionBandwidth = transitionBandwidth;
% paramsOut.maxPassbandRipple = maxPassbandRipple;
% paramsOut.filterOrder = filtOrder;
% paramsOut.KaiserWindowBeta = KaiserWindowBeta;
% paramsOut.indelecReject = indelecReject;
% paramsOut.icaweights = EEG.icaweights;
% paramsOut.icawinv = EEG.icawinv;
% paramsOut.icachansind = EEG.icachansind;
% paramsOut.icasphere = EEG.icasphere;
% paramsOut.indtrialRej = indtrialRej;
% 
% filename = [bdfFilepath, 'paramsOut_',mfilename, '_', subjID, '.mat'];
% save(filename, 'paramsOut');
end