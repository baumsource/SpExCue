% The only difference between this preprocessing script and the one without
% 'ICA' in the filename is that ICA is that ICA is also computed in this
% script. 

% [paramsOut] = dkrEEGLAB_pilot9_preprocessICA(subjID, expDate, condString)
% [paramsOut] = dkrEEGLAB_pilot9_preprocessICA(subjID, expDate, condString, FS)
% [paramsOut] = dkrEEGLAB_pilot9_preprocessICA(subjID, expDate, condString, FS, [locutoff, hicutoff])
% [paramsOut] = dkrEEGLAB_pilot9_preprocessICA(subjID, expDate, condString, FS, [locutoff, hicutoff], [transitionBW, maxRipple])
%
%  This function runs some basic preprocessing on the raw BDF file.
%  Namely, it: 
%   (1) re-references the data to the two external channels [33,34]
%   (2) removes the extra external channels [36:40], recall ch35 is ECOG
%   (3) defines the channel locations per an .locs file
%   (4) re-samples the data at defined sampling rate (DEFAULT: 512)
%   (5) filters the data given default or user input parameters
%
%  The eeglab data files are saved after Steps 2, 4, and 5.  Upone
%  completion of the preprocessing, this function also saves the parameters
%  of the preprocessing, i.e. 'paramsOut' as a .mat file with a timestamped
%  filename.
%
%  INPUTS:
%   subjID          STRING for the subject ID, e.g. '013'
%   expDate         STRING for the date of the experiment in the filename
%                   of the BDF file, e.g. '01112015'.
%   condString      STRING to identify if the data is for the attend
%                   background switch or attend figure swith condition. It
%                   MUST BE either 'Back' or 'Fig' and is CASE SENSITIVE!
%   FS              (optional) sampling frequency. Default is 512.
%   [locutoff, hicutoff]        ...
%                   (optional) filter cutoff frequencies in Hz. 
%                   Default is [0.5, 20].
%   [transitionBW, maxRipple]   ...
%                   (optional) other filter parameters.  
%                   Default is [1, 0.0002].
%
%  OUTPUT:
%   paramsOut       A structure with all of the important preprocessing
%                   settings.  NOTE: this structure is automatically saved
%                   upon completion of the function. The filename of the
%                   saved structure is timestamped to help identify the
%                   resulting eeglab data files with the .mat strucutre.
%
%  USAGE:
%       subjID = '013';  % STRING: for subject ID
%       expDate = '01112015';  % STRING: for date of experiment
%       condString = 'Fig';  % Either 'Back' or 'Fig'
%       FS = 512;
%       filtCutoffFreqs = [0.5 20];
%       filtParams = [1 0.0002];
%
%       [paramsOut] = dkrEEGLAB_pilot9_preprocessICA...
%                   (subjID, expDate, condString);
%    OR
%       [paramsOut] = dkrEEGLAB_pilot9_preprocessICA...
%                   (subjID, expDate, condString,...
%                   FS, filtCutoffFreqs, filtParams);
%
%
%  Created: Darrin K. Reed (Jan. 14, 2016)
%  Modified from 'dkrEEGLAB_pilot9_preprocess.m': 
%       Darrin K. Reed (Jan 31, 2016)... added the ICA analysis
                
%--------------------------------------------------------------------------

function [paramsOut] = analyzeEEGpilot3_eeglab_preprocessICA(subjID, expDate, condString, varargin)

if not(exist('subjID','var'))
  subjID = input('Subject ID: ','s');
end
  

% Primary path definitions required for this analysis
dirBDF = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Pilots/EEG3/data/';%string with the filename
dirLOCS = dirBDF;%'/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Tools/EEG/';

% dirBDF = 'C:\Users\dkreed\Documents\EEGdata\___TEST DATA\Exp2\';  % directory with the .bdf/.set files
% dirLOCS = '\\ad\eng\users\d\k\dkreed\Desktop\eeglabstuff\';       % directory with the eeglab locations file

% switch nargin
% %     case 3
%         FS = 512;  % In samps/sec
%         locutoff = 0.5;  % In Hz
%         hicutoff = 100;  % In Hz
%         transitionBandwidth = 1;  % In Hz
%         maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
%     
%     case 4
%         FS = varargin{1};
%         locutoff = 0.5;  % In Hz
%         hicutoff = 20;  % In Hz
%         transitionBandwidth = 1;  % In Hz
%         maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
% 
%     case 5
%         FS = varargin{1};
%         locutoff = varargin{2}(1);
%         hicutoff = varargin{2}(2);
%         transitionBandwidth = 1;  % In Hz
%         maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
%         
%     case 6
%         FS = varargin{1};
%         locutoff = varargin{2}(1);
%         hicutoff = varargin{2}(2);
%         transitionBandwidth = varargin{3}(1);
%         maxPassbandRipple = varargin{3}(2);
%         
%     otherwise
%         error('DKR: Invalid number of input parameters.')
% end

%--------------------------------------------------------------------------
% Input parameter checks
%--------------------------------------------------------------------------

% Check to ensure that condString is either 'Fig' or 'Back'
% if ~(strcmp(condString, 'Back') || strcmp(condString, 'Fig'))
%     error('DKR: The parameter ''condString'' MUST BE either ''Back'' or ''Fig'' (case sensitive).')
% end

% if locutoff >= hicutoff
%     error('DKR: Filter cutoff frequencies must be specified as [locutoff, hicutoff], e.g. [0.5 20]. Also, frequencies must not be equal.')
% end

experimentString = 'SpExCue_EEGpilot3';  % String to label experiment number for the saved eeglab dataset files
nChansLocsFilename = 'chanLocs_refA.locs';  % filename of the eeglab locations file

refChanNum = 32+[6,7];%[33 34]; % Indicies of reference channels in BDF file -> apply changes also to s
refChanLabel = 'A1 A2';%'TP9 TP10'; % Labels of reference channels
eyeChanNum = 35:37;
eegChanNum = 1:34;

FS = 200;  % In samps/sec
locutoff = 1;  % In Hz
% hicutoff = 100;  % In Hz
transitionBandwidth = 1;  % In Hz
maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002

globalThreshold = [-200 800];
stimulusEpochRange = [-1,1.9];

%% Open EEGlab
if not(exist('eeglab','file'))
  addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab')
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


%% Load BDF file
bdfFilepath = dirBDF;%[dirBDF, 'subj', subjID, '/'];  % NOTE: this path is VERY SPECIFIC TO THE USAGE
bdfFilename = [experimentString,'_',subjID];%['dkr_', experimentString, '_subj', subjID, '_', expDate, '_detect', condString, '.bdf'];  % NOTE: this filename is VERY SPECIFIC TO THE USAGE
fn = [bdfFilepath, bdfFilename];
EEG = pop_biosig([fn,'.bdf'], 'ref',refChanNum ,'refoptions',{'keepref' 'off'});
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
EEG = eeg_checkset(EEG);

%% Adjust event list (take only last 8 bits)
eventList = eeglabTrig2dkrTrigID([EEG.event.type]); 
for ii = 1:length(eventList)
	EEG.event(ii).type = eventList(ii);
end

%% Remove unused channels
EEG = pop_select(EEG,'nochannel',{'EXG8'});
EEG = pop_chanedit(EEG, 'load',{[dirLOCS, nChansLocsFilename] 'filetype' 'autodetect'},...  % Set channel locations
  'setref',{[num2str(refChanNum(1)),' ',num2str(refChanNum(2))] refChanLabel});             % and define reference channels
% [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset(EEG);
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'savenew',[fn,'_raw'],'gui','off'); 

%% Global artifact rejection for ICA

% Filter the data
KaiserWindowBeta = pop_kaiserbeta(maxPassbandRipple);
filtOrder = pop_firwsord('kaiser', FS, transitionBandwidth, maxPassbandRipple);
EEG = pop_firws(EEG, 'fcutoff', locutoff, 'ftype', 'highpass', 'wtype', 'kaiser', 'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);

% Epoching to stimulus onsets
eventSet = unique([EEG.event.type]);
tmp = eventSet(eventSet/10 == round(eventSet/10) & eventSet ~= 0);
stimulusOnsetEvents = cell(1,length(tmp));
for ii = 1:length(tmp)
  stimulusOnsetEvents{ii} = num2str(tmp(ii));
end
EEG = pop_epoch(EEG, stimulusOnsetEvents, stimulusEpochRange);
EEG = pop_rmbase(EEG,[-200,0]);
 
% Resample
EEG = pop_resample(EEG, FS);

% Threshold EEG channels
[EEG,indtrialRej1] = pop_eegthresh(EEG,1,eegChanNum,globalThreshold(1),globalThreshold(2),...
  stimulusEpochRange(1),stimulusEpochRange(2),0,1);
% Threshold eye channels
[EEG,indtrialRej2] = pop_eegthresh(EEG,1,eyeChanNum,-1*globalThreshold(2),-1*globalThreshold(1),...
  stimulusEpochRange(1),stimulusEpochRange(2),0,1);
indtrialRej = [indtrialRej1,indtrialRej2];

% Remove bad channels
[EEG, indelecReject] = pop_rejchan(EEG,'elec',eegChanNum);
% Update channel numbers
% indelecRemain = true(length(eegChanNum),1);
% indelecRemain(indelecReject) = false;
% eegChanNum = eegChanNum(indelecRemain);
% eyeChanNum = eyeChanNum - length(indelecReject);

% Save
EEG = eeg_checkset(EEG);
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'savenew',[fn,'_preICA'],'gui','off'); 

%% Resample the data to a specified sampling rate (DEFAULT: 512 samps/sec)
% EEG = pop_resample(EEG, FS);
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'savenew',[bdfFilepath, experimentString 'subj', subjID, '_detect', condString, '_resamp' num2str(FS)],'gui','off'); 

%% Bandpass filtering

% % Calculate the filter parameters  
% KaiserWindowBeta = pop_kaiserbeta(maxPassbandRipple);
% filtOrder = pop_firwsord('kaiser', FS, transitionBandwidth, maxPassbandRipple);
% 
% % Filter the data
% EEG = pop_firws(EEG, 'fcutoff', [locutoff hicutoff], 'ftype', 'bandpass', 'wtype', 'kaiser', 'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);
% fn = [fn,'_filt', num2str(hicutoff)];
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'savenew',fn,'gui','off'); 


%% Run the ICA analysis
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'extended',1,'interupt','on');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
fn = [fn,'_ICAraw'];
EEG = pop_saveset(EEG,  'filename', fn);%, 'filepath', bdfFilepath);


%% Reject components by map
% pop_selectcomps(EEG, [1:35] );
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% EEG = eeg_checkset( EEG );
% EEG = pop_subcomp( EEG );
% fn = strrep(fn,'raw','clean');
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'savenew',fn,'gui','off'); 


%% Close EEGlab
close all

%% Output and save the data parameters 
paramsOut.timestamp = clock;
paramsOut.subjID = subjID;
% paramsOut.expDate = expDate;
% paramsOut.condString = condString;
paramsOut.FS = FS;
paramsOut.locutoff = locutoff;
% paramsOut.hicutoff = hicutoff;
paramsOut.transitionBandwidth = transitionBandwidth;
paramsOut.maxPassbandRipple = maxPassbandRipple;
paramsOut.filterOrder = filtOrder;
paramsOut.KaiserWindowBeta = KaiserWindowBeta;
paramsOut.indelecReject = indelecReject;
paramsOut.icaweights = EEG.icaweights;
paramsOut.icawinv = EEG.icawinv;
paramsOut.icachansind = EEG.icachansind;
paramsOut.icasphere = EEG.icasphere;
paramsOut.indtrialRej = indtrialRej;

filename = [bdfFilepath, 'paramsOut_',mfilename, '_', subjID, '.mat'];
save(filename, 'paramsOut');
end