% The only difference between this preprocessing script and the one without
% 'ICA' in the filename is that ICA is that ICA is also computed in this
% script. 

% [paramsOut] = SpExCue_EEGLABpreprocessing(subjID, expDate)
% [paramsOut] = SpExCue_EEGLABpreprocessing(subjID, expDate, FS)
% [paramsOut] = SpExCue_EEGLABpreprocessing(subjID, expDate, FS, [locutoff, hicutoff])
% [paramsOut] = SpExCue_EEGLABpreprocessing(subjID, expDate, FS, [locutoff, hicutoff], transitionBW, maxRipple )
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
%
% Optional key value pairs:
%   FS = 100;
%   cutoffs = [0.5,25]; % filter cutoff frequencies in Hz.
%   transitionBandwidth = 1;  % In Hz
%   maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
%   experimentString = 'EEGpilot';  % String to label experiment number for the saved eeglab dataset files
%   dirBDF = '/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Experiments/Pilot/Results/';  % directory with the .bdf/.set files
%   dirLOCS = '/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Tools/EEG/';       % directory with the eeglab locations file
%   referenceChannels = [33,34]; % Indicies of reference channels in BDF file
%   nChansLocsFilename = 'biosemi_eloc.locs';  % filename of the eeglab locations file
%
% Optional flags:
%   ICA = {'','ICA'}; % run ICA or not
%
%  OUTPUT:
%   paramsOut       A structure with all of the important preprocessing
%                   settings.  NOTE: this structure is automatically saved
%                   upon completion of the function. The filename of the
%                   saved structure is timestamped to help identify the
%                   resulting eeglab data files with the .mat strucutre.
%
%
%  Created: Darrin K. Reed (Jan. 14, 2016)
%  Modified from 'dkrEEGLAB_exp2_preprocessICA.m': 
%       Darrin K. Reed (Jan 31, 2016)... added the ICA analysis
% Modified by Robert Baumgartner for SpExCue (Mar 8, 2016)
                
%--------------------------------------------------------------------------

function [paramsOut] = SpExCue_EEGLABpreprocessing(subjID, expDate, varargin)

definput.keyvals.FS = 100;
definput.keyvals.cutoff = [0.5,25];
definput.keyvals.transitionBandwidth = 1;  % In Hz
definput.keyvals.maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
definput.keyvals.experimentString = 'EEGpilot';  % String to label experiment number for the saved eeglab dataset files
definput.keyvals.dirBDF = '/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Experiments/Pilot/Results/';  % directory with the .bdf/.set files
definput.keyvals.dirLOCS = '/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Tools/EEG/';       % directory with the eeglab locations file
definput.keyvals.referenceChannels = [33,34]; % Indicies of reference channels in BDF file
definput.keyvals.nChansLocsFilename = 'biosemi_eloc.locs';  % filename of the eeglab locations file
definput.flags.ICA = {'','ICA'};

posdepnames = {'FS','cutoff','transitionBandwidth','maxPassbandRipple'};
[flags,kv]=ltfatarghelper(posdepnames,definput,varargin);

% switch nargin
%     case 3
%         FS = 512;  % In samps/sec
%         locutoff = 0.5;  % In Hz
%         hicutoff = 20;  % In Hz
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

%**************************************************************************
% Open EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%--------------------------------------------------------------------------
% Load BDF file
bdfFilepath = fullfile(kv.dirBDF,subjID,'');  % NOTE: this path is VERY SPECIFIC TO THE USAGE
bdfFilename = ['SpExCue_', kv.experimentString, '_', subjID, '_', expDate, '.bdf'];  % NOTE: this filename is VERY SPECIFIC TO THE USAGE
EEG = pop_biosig(fullfile(bdfFilepath, bdfFilename), 'ref',kv.referenceChannels ,'refoptions',{'keepref' 'off'});
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
EEG = eeg_checkset(EEG);

%--------------------------------------------------------------------------
% Remove unused channels
EEG = pop_select(EEG,'nochannel',{'EXG4', 'EXG5' 'EXG6' 'EXG7' 'EXG8'});
EEG = pop_chanedit(EEG, 'load',{[dirLOCS, kv.nChansLocsFilename] 'filetype' 'autodetect'});  % Set channel locations
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset(EEG);
fnprefix = [kv.experimentString '_subj', subjID];
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'savenew',fullfile(bdfFilepath,[fnprefix,'_raw']),'gui','off'); 

%--------------------------------------------------------------------------
% Resample the data to a specified sampling rate (DEFAULT: 512 samps/sec)
EEG = pop_resample(EEG, kv.FS);
fn = [fnprefix,'_resamp' num2str(kv.FS)];
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'savenew',fullfile(bdfFilepath,fn),'gui','off'); 

%--------------------------------------------------------------------------
% Filter the data
  
KaiserWindowBeta = pop_kaiserbeta(kv.maxPassbandRipple);% Calculate the filter parameters
filtOrder = pop_firwsord('kaiser', kv.FS, kv.transitionBandwidth, kv.maxPassbandRipple);

EEG = pop_firws(EEG, 'fcutoff', kv.cutoff, 'ftype', 'bandpass', 'wtype', 'kaiser', 'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);
fn = [fn,'_filt', num2str(kv.cutoff(2))];
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'savenew',fullfile(bdfFilepath,fn),'gui','off'); 

%--------------------------------------------------------------------------
% Run the ICA analysis
if flags.do_ICA
  [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
  EEG = eeg_checkset( EEG );
  EEG = pop_runica(EEG, 'extended',1,'interupt','on');
  [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
  fn = [fn,'_ICAraw'];
  EEG = pop_saveset(EEG,  'filename', fn, 'filepath', bdfFilepath);
end

%--------------------------------------------------------------------------
% Close EEGlab
close all

%--------------------------------------------------------------------------
% Output and save the data parameters 
paramsOut.timestamp = clock;
paramsOut.subjID = subjID;
paramsOut.expDate = expDate;
paramsOut.kv = kv;
paramsOut.flags = flags;
% paramsOut.FS = FS;
% paramsOut.locutoff = locutoff;
% paramsOut.hicutoff = hicutoff;
% paramsOut.transitionBandwidth = transitionBandwidth;
% paramsOut.maxPassbandRipple = maxPassbandRipple;
paramsOut.filterOrder = filtOrder;
paramsOut.KaiserWindowBeta = KaiserWindowBeta;

filename = fullfile(bdfFilepath, [kv.experimentString, '_subj', subjID, '_paramsOut', num2str(paramsOut.timestamp(1)), '_', num2str(paramsOut.timestamp(2)), '_', num2str(paramsOut.timestamp(3)), '_', num2str(paramsOut.timestamp(4)), '_', num2str(paramsOut.timestamp(5)), '.mat']);
save(filename, 'paramsOut');
