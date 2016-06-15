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

% Primary path definitions required for this analysis
dirBDF = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Pilots/EEG3/data/';%string with the filename
dirLOCS = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Tools/EEG/';

% dirBDF = 'C:\Users\dkreed\Documents\EEGdata\___TEST DATA\Exp2\';  % directory with the .bdf/.set files
% dirLOCS = '\\ad\eng\users\d\k\dkreed\Desktop\eeglabstuff\';       % directory with the eeglab locations file

% switch nargin
%     case 3
        FS = 512;  % In samps/sec
        locutoff = 0.5;  % In Hz
        hicutoff = 100;  % In Hz
        transitionBandwidth = 1;  % In Hz
        maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
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

if locutoff >= hicutoff
    error('DKR: Filter cutoff frequencies must be specified as [locutoff, hicutoff], e.g. [0.5 20]. Also, frequencies must not be equal.')
end


experimentString = 'SpExCue_EEGpilot3';  % String to label experiment number for the saved eeglab dataset files
referenceChannels = [33 34]; % Indicies of reference channels in BDF file
nChansLocsFilename = 'biosemi_eloc.locs';  % filename of the eeglab locations file

%**************************************************************************
% Open EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%--------------------------------------------------------------------------
% Load BDF file
bdfFilepath = dirBDF;%[dirBDF, 'subj', subjID, '/'];  % NOTE: this path is VERY SPECIFIC TO THE USAGE
bdfFilename = [experimentString,'_',subjID];%['dkr_', experimentString, '_subj', subjID, '_', expDate, '_detect', condString, '.bdf'];  % NOTE: this filename is VERY SPECIFIC TO THE USAGE
fn = [bdfFilepath, bdfFilename];
EEG = pop_biosig([fn,'.bdf'], 'ref',referenceChannels ,'refoptions',{'keepref' 'off'});
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
EEG = eeg_checkset(EEG);

%--------------------------------------------------------------------------
% Remove unused channels
EEG = pop_select(EEG,'nochannel',{'EXG4', 'EXG5' 'EXG6' 'EXG7' 'EXG8'}); % adjust for RS !!!!
EEG = pop_chanedit(EEG, 'load',{[dirLOCS, nChansLocsFilename] 'filetype' 'autodetect'});  % Set channel locations
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset(EEG);
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'savenew',[fn,'_raw'],'gui','off'); 

%--------------------------------------------------------------------------
% Resample the data to a specified sampling rate (DEFAULT: 512 samps/sec)
% EEG = pop_resample(EEG, FS);
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'savenew',[bdfFilepath, experimentString 'subj', subjID, '_detect', condString, '_resamp' num2str(FS)],'gui','off'); 

%--------------------------------------------------------------------------
% Calculate the filter parameters  
KaiserWindowBeta = pop_kaiserbeta(maxPassbandRipple);
filtOrder = pop_firwsord('kaiser', FS, transitionBandwidth, maxPassbandRipple);

%--------------------------------------------------------------------------
% Filter the data
EEG = pop_firws(EEG, 'fcutoff', [locutoff hicutoff], 'ftype', 'bandpass', 'wtype', 'kaiser', 'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);
fn = [fn,'_filt', num2str(hicutoff)];
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'savenew',fn,'gui','off'); 


%--------------------------------------------------------------------------
% Run the ICA analysis
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'extended',1,'interupt','on');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
fn = [fn,'_ICAraw'];
EEG = pop_saveset(EEG,  'filename', fn);%, 'filepath', bdfFilepath);


% EEG = eeg_checkset( EEG );
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
% EEG = eeg_checkset( EEG );
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'savenew',[bdfFilepath, experimentString 'subj', subjID, '_detect', condString, '_resamp', num2str(FS), 'filt', num2str(hicutoff),'_ICAraw'],'gui','off'); 
% EEG = pop_saveset( EEG, 'filename','exp1subj021_detectFig_resamp512filt20_ICAclean.set','filepath','C:\\Users\\dkreed\\Documents\\EEGdata\\___TEST DATA\\Exp1\\subj021\\');


%--------------------------------------------------------------------------
% Close EEGlab
close all

%--------------------------------------------------------------------------
% Output and save the data parameters 
paramsOut.timestamp = clock;
paramsOut.subjID = subjID;
% paramsOut.expDate = expDate;
% paramsOut.condString = condString;
paramsOut.FS = FS;
paramsOut.locutoff = locutoff;
paramsOut.hicutoff = hicutoff;
paramsOut.transitionBandwidth = transitionBandwidth;
paramsOut.maxPassbandRipple = maxPassbandRipple;
paramsOut.filterOrder = filtOrder;
paramsOut.KaiserWindowBeta = KaiserWindowBeta;

filename = [bdfFilepath, experimentString, 'Subj', subjID, '_paramsOut', num2str(paramsOut.timestamp(1)), '_', num2str(paramsOut.timestamp(2)), '_', num2str(paramsOut.timestamp(3)), '_', num2str(paramsOut.timestamp(4)), '_', num2str(paramsOut.timestamp(5)), '.mat'];
save(filename, 'paramsOut');
