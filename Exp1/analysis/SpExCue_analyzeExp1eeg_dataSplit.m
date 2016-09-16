function SpExCue_analyzeExp1eeg_dataSplit(varargin)
% Analyze SpExCue_EEGpilot

definput.keyvals.filenameIDx = 'Exp1eeg_IDx_ICAraw.set';
definput.keyvals.filepath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/data/';
definput.keyvals.IDs = []; % subject IDs (as cell array)
definput.keyvals.epochInterval = [-1,2]; % for export
definput.keyvals.erpInterval = [-.2,.8]; % for plotting and thresholding
definput.keyvals.baselineInterval = [-.200,0]; % for baseline removal and thresholding
definput.keyvals.chanNum = 1:32; % channel numbers
% definput.keyvals.ERPfmax = 30; % cut-off frequency of low-pass filter for ERPs
% Settings for automatic epoch rejection
definput.keyvals.epochThresh = Inf; % threshold in µV within entire epoch
definput.keyvals.baselineThresh = Inf; % threshold in µV within baseline interval
% definput.keyvals.slopeThresh = 80; % max slope in µV/epoch
% definput.keyvals.stepThresh = 50; % threshold in µV for step function based artifact rejection
% definput.keyvals.latencyRanges = [[.075,.150];[.175,.250];[.300,.600]]; % N1, P2, P3 latency ranges in seconds
% definput.keyvals.latencyLabels = {'N1 (75-150ms)','P2 (175-250ms)','P3 (300-600ms)'};

% definput.flags.site = {'central','frontal','parietal','topo'}; 
definput.flags.timeLock = {'Change','Onset'}; 
definput.flags.position = {'all','left','right','top'}; 
% definput.flags.contralaterality = {...
%   'off';...
% %   'ipsiVsContra';...
%   'rightHemiDominance';... % Getzmann & Lewald (2010)
% %   'intrahemiContralaterality';... % Palomäki et al. (2005)
%   };
% definput.flags.TFanalysis = {'','TFanalysis'};
% definput.flags.plotersp = {'on','off'};
% definput.flags.plotitc = {'on','off'};
% definput.flags.eyeChan = {'bipolarEyeChan','monopolarEyeChan'};
definput.flags.grouping = {'noGrouping','groupM','groupD','groupE','groupE_Mp3'};
% definput.flags.print = {'','print'};
definput.flags.export = {'export','noExport'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

%% Check missing variables and functions

if isempty(kv.IDs)
  IDs = {input('Subject ID: ','s')};
elseif ischar(kv.IDs)
  IDs = {kv.IDs};
else
  IDs = kv.IDs;
end

if not(exist('pop_loadset','file'))
  addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab/')
  eeglab
end

if not(exist('pop_artstep_EEGlab','file'))
  addpath('/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/MATLAB_general/');
end

%% Settings

switch flags.timeLock

  case 'Onset'
    MlegendLabel = {'M_o = 1  ','M_o = 1/3','M_o = 0  '};
    timeLockedEventNum = {30,20,10};
    M = [1,1/3,0];
    
  case 'Change'
    
%     MlegendLabel = {...
%       '1   \rightarrow 0  ','1   \rightarrow 1/3','1/3 \rightarrow 0',...
%       '0   \rightarrow 0  ','1/3 \rightarrow 1/3','1   \rightarrow 1',...
%       '0   \rightarrow 1/3','1/3 \rightarrow 1  ','0   \rightarrow 1'};
    timeLockedEventNum = {31,32,21,11,22,33,12,23,13};
    MlegendLabel = {'M1-0','M1-i','Mi-0',...
                    'M0-0','Mi-i','M1-1',...
                    'M0-i','Mi-1','M0-1'};
    
end

if flags.do_groupD
  grouping = {31,32,21,[11,22,33],12,23,13};
  MlegendLabel = {'D =-1   ','D =-2/3','D =-1/3','D = 0   ','D = 1/3','D = 2/3','D = 1   '};
  
elseif flags.do_groupM
  grouping = {[33,23,13],[32,22,12],[31,21,11]};
  MlegendLabel = {'M_c = 1  ','M_c = 1/3','M_c = 0  '};
  
elseif flags.do_groupE || flags.do_groupE_Mp3
  grouping = {5,6};
  MlegendLabel = {'Closer','Farther'};
  if flags.do_groupE_Mp3
    timeLockedEventNum = {32,12};
    flags.do_groupE = true;
  else
    timeLockedEventNum = {31,32,21,12,23,13}; % exclude D = 0
  end
  
else
  grouping = timeLockedEventNum;
  
end

NumTrials = nan(length(IDs),length(grouping));
for ll = 1:length(IDs)
  
  ID = IDs{ll};

  %% Load the actual data
  
  filename = strrep(kv.filenameIDx,'IDx',ID);
  EEG = pop_loadset('filename', filename, 'filepath', kv.filepath); 
  EEG = eeg_checkset(EEG);  % verify consistency of the fields of an EEG dataset

  %% Event selection
  
  % Remove first digit coding stimulus direction
  if flags.do_all 
    for ii = 1:length(EEG.event)
      EEG.event(ii).type = mod(EEG.event(ii).type,100);
    end
  end

  % Select Epochs and correct for baseline offset
  EEGall = pop_epoch(EEG,timeLockedEventNum,kv.epochInterval);
  event_indices = cell(length(grouping),1);
  for ii = 1:length(grouping)
    [EEG(ii),event_indices{ii}] = pop_selectevent(EEGall,'type',grouping{ii});
    EEG(ii) = pop_rmbase(EEG(ii), 1e3*kv.baselineInterval);
    EEG(ii) = eeg_checkset(EEG(ii));
  end
  if length(grouping) == 2 && sum(ismember(event_indices{1},event_indices{2})) > 0
    error('RB: There are overlaping epochs between groups. Reduce epoch interval!')
  end
  
  %% Artifactual epoch rejection
  for ii = 1:length(EEG)

%     Npre = EEG(ii).trials; % initial trial count

    % Threshold rejection
    EEG(ii) = pop_eegthresh(EEG(ii),1,kv.chanNum,-kv.epochThresh,kv.epochThresh,...
      kv.erpInterval(1),kv.erpInterval(2),0,1);
%     Nrej.epochTresh(ll) = Nrej.epochTresh(ll) + (Npre-EEG(ii).trials);
%     Ntmp = EEG(ii).trials;
    if kv.epochThresh > kv.baselineThresh
      EEG(ii) = pop_eegthresh(EEG(ii),1,kv.chanNum,-kv.baselineThresh,kv.baselineThresh,...
        kv.baselineInterval(1),kv.baselineInterval(2),0,1);
%       Nrej.baselineThresh(ll) = Nrej.baselineThresh(ll) + (Ntmp-EEG(ii).trials);
%       Ntmp = EEG(ii).trials;
    end
%     Npost = EEG(ii).trials;
 
%     % Slope-based rejection
%     if not(isnan(kv.slopeThresh))
%       EEG(ii) = pop_rejtrend(EEG(ii),1,chNum,512,kv.slopeThresh,0.3,0,1,0);
%       Nrej.slopeThresh(ll) = Nrej.slopeThresh(ll) + (Ntmp-EEG(ii).trials);
%       Ntmp = EEG(ii).trials;
%     end
% 
%     % Step-function-based rejection
%     if flags.do_bipolarEyeChan
%       EEG(ii) = pop_artstep_EEGlab(EEG(ii),1e3*kv.erpInterval,kv.stepThresh,200,50,eyeChan,1);
%       Nrej.stepThresh(ll) = Nrej.stepThresh(ll) + (Ntmp-EEG(ii).trials);
%     end
%     Npost = EEG(ii).trials;
%     
%     if Npost/Npre < 0.75
%       disp('*****')
%       warning('  RB: More than 25% of trials removed for ',ID,'!')
%       disp('*****')
%     end
    
    EEG(ii) = eeg_checkset(EEG(ii));
  end
    
    %% Dipole fitting
    
    %% Check # trials
    
    
    %% Export dataset
  for ii = 1:length(EEG)
    % Check # trials
    NumTrials(ll,ii) = EEG(ii).trials;
    
    if flags.do_export
      fnExport = strrep(filename,'.set',['_',MlegendLabel{ii},'.set']);
  %     fnExport = regexprep(fnExport,'^[\d\w~!@#$%^&()_\-{}]',''); % problem: removes first character
      EEG(ii).subject = ID;
      EEG(ii).condition = MlegendLabel{ii};
      pop_saveset(EEG(ii),'filename', fullfile(kv.filepath,fnExport));
    end
  end
end

NumTrialsTab = array2table(NumTrials,'RowNames',IDs,'VariableNames',strrep(MlegendLabel,'-','to'));
disp(NumTrialsTab)

end