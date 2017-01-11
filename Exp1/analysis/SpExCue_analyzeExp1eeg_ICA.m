function [paramsOut] = SpExCue_analyzeExp1eeg_ICA(subjID)

if not(exist('subjID','var'))
  subjID = input('Subject ID: ','s');
end
  
% Primary path definitions required for this analysis
dirBDF = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/data/';%string with the filename
dirLOCS = dirBDF;%'/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Tools/EEG/';

experimentString = 'Exp1eeg';  % String to label experiment number for the saved eeglab dataset files
nChansLocsFilename = 'chanLocs_refA.locs';  % filename of the eeglab locations file

refChanNum = 32+[6,7];%[33 34]; % Indicies of reference channels in BDF file -> apply changes also to s
refChanLabel = 'A1 A2';%'TP9 TP10'; % Labels of reference channels
eyeChanNum = 35:37;
eegChanNum = 1:34;

FS = 100;  % In samps/sec
locutoff = 0.5;  % In Hz
hicutoff = 20;  % In Hz
transitionBandwidth = 1;  % In Hz
maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002

globalThreshold = [-200 800];
stimulusEpochRange = [-1,3];

timeLockedEventNum = {31,32,21,11,22,33,12,23,13};

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
% check if bad channels marked during recording
if exist([fn,'_badchan.mat'],'file')
  tmp = load([fn,'_badchan.mat']);
  EEG = eeg_interp(EEG,tmp.badchan,'spherical');
%   for ii=1:length(tmp.badchan);
%     EEG.chanlocs(tmp.badchan(ii)).badchan = true;
%   end
%   EEG = pop_select( EEG,'nochannel',{EEG.chanlocs(find([EEG.chanlocs.badchan])).labels});
%   EEG = eeg_interp(EEG,find([EEG.chanlocs.badchan]),'spherical');
end
EEG = eeg_checkset(EEG);
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'savenew',[fn,'_raw'],'gui','off'); 

%% Global artifact rejection for ICA

% Filter the data
KaiserWindowBeta = pop_kaiserbeta(maxPassbandRipple);
filtOrder = pop_firwsord('kaiser', FS, transitionBandwidth, maxPassbandRipple);
EEG = pop_firws(EEG, 'fcutoff', [locutoff,hicutoff], 'wtype', 'kaiser',...
  'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);

% Remove first digit coding stimulus direction
for ii = 1:length(EEG.event)
  EEG.event(ii).type = mod(EEG.event(ii).type,100);
end

% Epoching to stimulus change
EEG = pop_epoch(EEG, timeLockedEventNum, stimulusEpochRange);
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
if indelecReject
  disp(['Automatically detected bad channels for ',subjID,':'])
  disp(indelecReject)
  pause
end

% Save
EEG = eeg_checkset(EEG);
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'savenew',[fn,'_preICA'],'gui','off'); 


%% Run the ICA analysis
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'extended',1,'interupt','on','icatype','binica');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
fn = [fn,'_ICAraw'];
EEG = pop_saveset(EEG,  'filename', fn);%, 'filepath', bdfFilepath);


%% Close EEGlab
close all

%% Output and save the data parameters 
paramsOut.timestamp = clock;
paramsOut.subjID = subjID;
paramsOut.FS = FS;
paramsOut.locutoff = locutoff;
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