function [fnsave,ALLEEG, EEG] = eeglab_preICA(cnfg,fn_bdf)
% eeglab_preICA - preprocessing of EEG data including referencing,
% adjustment of event list, removal of unused channels, bandpass filtering,
% epoching, resampling, thresholding, removal of bad channels
%
% Usage: eeglab_preICA(cnfg,fn_bdf)

% AUTHOR: Robert Baumgartner, 2017/02/27

ALLEEG = eeglab;

%% Load BDF file
fn = fn_bdf;
if not(exist([fn,'_raw.set'],'file'))
  EEG = pop_biosig([fn,'.bdf'], 'ref',cnfg.refChanNum ,'refoptions',{'keepref' 'off'});
  [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
  EEG = eeg_checkset(EEG);

  %% Adjust event list (take only last 8 bits)
  eventList = eeglabTrig2dkrTrigID([EEG.event.type]); 
  for ii = 1:length(eventList)
    EEG.event(ii).type = eventList(ii);
  end
  if strcmp(fn(end-2:end),'S27')
    EEG = fixEventsS27(EEG);
  end

  %% Remove unused channels
  EEG = pop_select(EEG,'nochannel',{'EXG8'});
  EEG = pop_chanedit(EEG, 'load',{fullfile(cnfg.dirLOCS, cnfg.nChansLocsFilename) 'filetype' 'autodetect'},...  % Set channel locations
    'setref',{[num2str(cnfg.refChanNum(1)),' ',num2str(cnfg.refChanNum(2))] cnfg.refChanLabel});             % and define reference channels
  % check if bad channels marked during recording
  if exist([fn,'_badchan.mat'],'file')
    tmp = load([fn,'_badchan.mat']);
    EEG = eeg_interp(EEG,tmp.badchan,'spherical');
  end
  EEG = eeg_checkset(EEG);
  [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'savenew',[fn,'_raw'],'gui','off'); 
else
  EEG = pop_loadset('filename', [fn,'_raw.set']);
end


%% Global artifact rejection for ICA

% Filter the data
KaiserWindowBeta = pop_kaiserbeta(cnfg.maxPassbandRipple);
filtOrder = pop_firwsord('kaiser', cnfg.FS, cnfg.transitionBandwidth, cnfg.maxPassbandRipple);
EEG = pop_firws(EEG, 'fcutoff', [cnfg.locutoff,cnfg.hicutoff], 'wtype', 'kaiser',...
  'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);

% Remove first digit coding stimulus direction
% for ii = 1:length(EEG.event)
%   EEG.event(ii).type = mod(EEG.event(ii).type,100);
% end

% Epoching to specified events
EEG = pop_epoch(EEG, cnfg.timeLockedEventNum, cnfg.stimulusEpochRange);
EEG = pop_rmbase(EEG,cnfg.baselineRange);
 
% Resample
EEG = pop_resample(EEG, cnfg.FS);

% Threshold EEG channels
[EEG,indtrialRej1] = pop_eegthresh(EEG,1,cnfg.eegChanNum,...
  cnfg.globalThreshold(1),cnfg.globalThreshold(2),...
  cnfg.stimulusEpochRange(1),cnfg.stimulusEpochRange(2),0,1);
% Threshold eye channels
[EEG,indtrialRej2] = pop_eegthresh(EEG,1,cnfg.eyeChanNum,...
  -1*cnfg.globalThreshold(2),-1*cnfg.globalThreshold(1),...
  cnfg.stimulusEpochRange(1),cnfg.stimulusEpochRange(2),0,1);
indtrialRej = [indtrialRej1,indtrialRej2];

% Remove bad channels
[EEG, indelecReject] = pop_rejchan(EEG,'elec',cnfg.eegChanNum);
if indelecReject
  disp(['Automatically detected bad channels for ',fn_bdf,':'])
  disp(indelecReject)
  pause
  save([fn,'_indelecReject'],indelecReject)
end

% Save
EEG = eeg_checkset(EEG);
fnsave = [fn,'_preICA'];
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'savenew',fnsave,'gui','off');

end