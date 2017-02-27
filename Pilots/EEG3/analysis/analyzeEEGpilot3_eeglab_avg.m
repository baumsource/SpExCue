function analyzeEEGpilot3_eeglab_avg(varargin)
% Analyze SpExCue_EEGpilot

definput.keyvals.filenameIDx = 'SpExCue_EEGpilot3_IDx_blinkICrej_fmax100.set';
definput.keyvals.filepath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Pilots/EEG3/data/';
definput.keyvals.IDs = []; % subject IDs (as cell array)
definput.keyvals.epochInterval = [-1,3]; % for exposrt
definput.keyvals.erpInterval = [-.2,.8]; % for plotting and thresholding
definput.keyvals.baselineInterval = [-.200,0]; % for baseline removal and thresholding
definput.keyvals.ERPfmax = 30; % cut-off frequency of low-pass filter for ERPs
% Settings for automatic epoch rejection
definput.keyvals.epochThresh = 80; % threshold in µV within entire epoch
definput.keyvals.baselineThresh = 50; % threshold in µV within baseline interval
definput.keyvals.slopeThresh = 80; % max slope in µV/epoch
definput.keyvals.stepThresh = 50; % threshold in µV for step function based artifact rejection
definput.keyvals.latencyRanges = [[.075,.150];[.175,.250];[.300,.600]]; % N1, P2, P3 latency ranges in seconds
definput.keyvals.latencyLabels = {'N1 (75-150ms)','P2 (175-250ms)','P3 (300-600ms)'};

definput.flags.site = {'central','frontal','parietal','topo'}; 
definput.flags.timeLock = {'Change','Onset'}; 
definput.flags.position = {'all','left','right','top'}; 
definput.flags.contralaterality = {...
  'off';...
%   'ipsiVsContra';...
  'rightHemiDominance';... % Getzmann & Lewald (2010)
%   'intrahemiContralaterality';... % Palomäki et al. (2005)
  };
definput.flags.TFanalysis = {'','TFanalysis'};
definput.flags.plotersp = {'on','off'};
definput.flags.plotitc = {'on','off'};
definput.flags.eyeChan = {'bipolarEyeChan','monopolarEyeChan'};
definput.flags.grouping = {'noGrouping','groupM','groupD','groupE','groupE_Mp3'};
definput.flags.print = {'','print'};
definput.flags.export = {'','export'};

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

Color = 0.75*colormap(hsv(3));

switch flags.timeLock

  case 'Onset'
    MlegendLabel = {'M_o = 1  ','M_o = 1/3','M_o = 0  '};
    timeLockedEventNum = {30,20,10};
    LineStyle = {'-','-.',':'};
    M = [1,1/3,0];
    
  case 'Change'
    
    MlegendLabel = {...
      '1   \rightarrow 0  ','1   \rightarrow 1/3','1/3 \rightarrow 0',...
      '0   \rightarrow 0  ','1/3 \rightarrow 1/3','1   \rightarrow 1',...
      '0   \rightarrow 1/3','1/3 \rightarrow 1  ','0   \rightarrow 1'};
    LineStyle = {'-','-','-',':',':',':','-.','-.','-.'};
    Color = [Color;Color;flipud(Color)];
    timeLockedEventNum = {31,32,21,11,22,33,12,23,13};
    MlegendLabel = {'M1-0','M1-i','Mi-0',...
                    'M0-0','Mi-i','M1-1',...
                    'M0-i','Mi-1','M0-1'};
    
end

if flags.do_groupD
  grouping = {31,32,21,[11,22,33],12,23,13};
  MlegendLabel = {'D =-1   ','D =-2/3','D =-1/3','D = 0   ','D = 1/3','D = 2/3','D = 1   '};
  LineStyle = {'-','-','-',':','-.','-.','-.'};
  Color = Color([1:4,7:9],:);
  
elseif flags.do_groupM
  grouping = {[33,23,13],[32,22,12],[31,21,11]};
  MlegendLabel = {'M_c = 1  ','M_c = 1/3','M_c = 0  '};
  LineStyle = {'-','-.',':'};
  
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

if flags.do_TFanalysis
  
  kv.epochInterval = [-3,3];

  % Frequency bands for T-F-Analysis
  freqBand(1).name = 'delta';
  freqBand(1).range = [2,4]; % Hz
  freqBand(2).name = 'theta';
  freqBand(2).range = [4,8]; % Hz
  freqBand(3).name = 'alpha';
  freqBand(3).range = [8,12]; % Hz
  freqBand(4).name = 'beta';
  freqBand(4).range = [15,30]; % Hz
  
else

  populationMean = cell(length(grouping),1);
  
end

% Evaluation site
if strcmp(flags.contralaterality,'off')
  switch flags.site
    case 'frontal'
      chNum = 31; % Fz
    case 'central'
      chNum = 32; % Cz
    case 'parietal'
      chNum = 13; % Pz
    case 'topo'
      chNum = 1:32;
  end
else
  chNum = 1:32;
end

% Initialize epoch rejection counter
Nrej.epochTresh = zeros(1,length(IDs));
Nrej.baselineThresh = zeros(1,length(IDs));
Nrej.slopeThresh = zeros(1,length(IDs));
Nrej.stepThresh = zeros(1,length(IDs));

for ll = 1:length(IDs)
  
  ID = IDs{ll};

  %% Load the actual data
  
  filename = strrep(kv.filenameIDx,'IDx',ID);
  if flags.do_TFanalysis % -> load data with low-pass at relatively high freq
    EEG = pop_loadset('filename', filename, 'filepath', kv.filepath); 
  else % ERP -> use low-pass filtered data
    fnERP = strrep(filename,'fmax100',['fmax',num2str(kv.ERPfmax)]);
    try
      EEG = pop_loadset('filename', fnERP, 'filepath', kv.filepath);
      
    catch % -> perform low-pass filterting and store result for next run
      EEG = pop_loadset('filename', filename, 'filepath', kv.filepath); 
      % Low-pass filter for ERP analysis
      transitionBandwidth = 1;  % In Hz
      maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
      KaiserWindowBeta = pop_kaiserbeta(maxPassbandRipple);
      filtOrder = pop_firwsord('kaiser', EEG.srate, transitionBandwidth, maxPassbandRipple);
      EEG = pop_firws(EEG, 'fcutoff', kv.ERPfmax, 'ftype', 'lowpass', 'wtype', 'kaiser', 'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);
      EEG = pop_saveset(EEG,  'filename', fullfile(kv.filepath,fnERP));
    end
    
  end
  EEG = eeg_checkset(EEG);  % verify consistency of the fields of an EEG dataset

  %% Create bipolar eye channels
  if flags.do_bipolarEyeChan
    EEG = pop_eegchanoperator( EEG,... 
      {'ch35=ch30-ch35 label vEOG',...
       'ch36=ch36-ch37 label RHEOG',...
       'ch37=ch37-ch36 label LRVEOG'});
    eyeChan = 35:37;
  end

  %% Event selection
  
  % Remove first digit coding stimulus direction
  if flags.do_all 
    for ii = 1:length(EEG.event)
      EEG.event(ii).type = mod(EEG.event(ii).type,100);
    end
  end
  
  % Include externalization response trigger if not existent (for RB, S01, and S04)
  % and correct stimulus offset trigger (0 instead of 9)
  NcorrectedEventTypes = 0;
  if flags.do_groupE && not(ismember(5,[EEG.event.type]))
    for ii = 3:length(EEG.event)
      if ismember(EEG.event(ii-2).type,[timeLockedEventNum{:}]) % previous event was M change
        tmp = num2str(EEG.event(ii-2).type); % change to string for digit separation
        D = str2double(tmp(2)) - str2double(tmp(1)); % subtract digits to get D
        if (D > 0 && EEG.event(ii).type == 8)... % D>0 and correct
            || (D < 0 && EEG.event(ii).type == 9) % D<0 and wrong
          EEG.event(ii).type = 6; % farther
        else
          EEG.event(ii).type = 5; % closer
        end
        % correct offset trigger
        if EEG.event(ii-1).type == 9
          EEG.event(ii-1).type = 0;
        end
        NcorrectedEventTypes = NcorrectedEventTypes+1;
      end
    end
  end
  disp([ID,': ',num2str(NcorrectedEventTypes),' corrected event types.'])

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

    Npre = EEG(ii).trials; % initial trial count

    % Threshold rejection
    EEG(ii) = pop_eegthresh(EEG(ii),1,chNum,-kv.epochThresh,kv.epochThresh,...
      kv.erpInterval(1),kv.erpInterval(2),0,1);
    Nrej.epochTresh(ll) = Nrej.epochTresh(ll) + (Npre-EEG(ii).trials);
    Ntmp = EEG(ii).trials;
    if kv.epochThresh > kv.baselineThresh
      EEG(ii) = pop_eegthresh(EEG(ii),1,chNum,-kv.baselineThresh,kv.baselineThresh,...
        kv.baselineInterval(1),kv.baselineInterval(2),0,1);
      Nrej.baselineThresh(ll) = Nrej.baselineThresh(ll) + (Ntmp-EEG(ii).trials);
      Ntmp = EEG(ii).trials;
    end

    % Slope-based rejection
    if not(isnan(kv.slopeThresh))
      EEG(ii) = pop_rejtrend(EEG(ii),1,chNum,512,kv.slopeThresh,0.3,0,1,0);
      Nrej.slopeThresh(ll) = Nrej.slopeThresh(ll) + (Ntmp-EEG(ii).trials);
      Ntmp = EEG(ii).trials;
    end

    % Step-function-based rejection
    if flags.do_bipolarEyeChan
      EEG(ii) = pop_artstep_EEGlab(EEG(ii),1e3*kv.erpInterval,kv.stepThresh,200,50,eyeChan,1);
      Nrej.stepThresh(ll) = Nrej.stepThresh(ll) + (Ntmp-EEG(ii).trials);
    end
    Npost = EEG(ii).trials;
    
    if Npost/Npre < 0.75
      disp('*****')
      warning('  RB: More than 25% of trials removed for ',ID,'!')
      disp('*****')
    end

    EEG(ii) = eeg_checkset(EEG(ii));
    
    %% Dipole fitting
    
    
    %% Export dataset
    if flags.do_export
      fnExport = strrep(fnERP,'.set',['_',MlegendLabel{ii},'.set']);
  %     fnExport = regexprep(fnExport,'^[\d\w~!@#$%^&()_\-{}]',''); % problem: removes first character
      EEG(ii).subject = ID;
      EEG(ii).condition = MlegendLabel{ii};
      pop_saveset(EEG(ii),'filename', fullfile(kv.filepath,fnExport));
    end

  end
  
%   disp(Nrej); pause;

  %% Average epochs and evaluate component measures
  switch flags.contralaterality

    case 'rightHemiDominance'

      idLH = ismember({EEGall.chanlocs.labels},... % left hemisphere
        {'AF3','AF7','C1','C3','C5','CP1','CP3','CP5','F3','F7','FC1','FC3',...
        'FC5','Fp1','FT9','O1','P1','P3','P7','PO3','PO7','PO9','T7','TP7'});
      idRH = ismember({EEGall.chanlocs.labels},... % right hemisphere
        {'AF4','AF8','C2','C4','C6','CP2','CP4','CP6','F4','F8','FC2','FC4',...
        'FC6','Fp2','FT10','O2','P2','P4','P8','PO10','PO4','PO8','T8','TP8'});
      time = kv.epochInterval(1):1/EEGall.srate:kv.epochInterval(2)-1/EEGall.srate;
      idTimeN1 = time>=kv.latencyRanges(1,1) & time<=kv.latencyRanges(1,2);
      idTimeP2 = time>=kv.latencyRanges(2,1) & time<=kv.latencyRanges(2,2);
      for ii = 1:length(EEG)
        N1amp = min(mean(EEG(ii).data(:,idTimeN1,:),3),[],2);
        P2amp = max(mean(EEG(ii).data(:,idTimeP2,:),3),[],2);
        N1RH(ii) = rms(N1amp(idRH));
        N1LH(ii) = rms(N1amp(idLH));
        P2RH(ii) = rms(P2amp(idRH));
        P2LH(ii) = rms(P2amp(idLH));
      end
      N1RHdom = N1RH - N1LH; % right-hemispheric dominance
      P2RHdom = P2RH - P2LH; % right-hemispheric dominance

    otherwise

      if flags.do_TFanalysis
        ERSP = [];
        for ii = 1:length(EEG)
          if strcmp(kv.plotersp,'on') || strcmp(kv.plotitc,'on')
            figure
          end
          [ERSP(:,:,ii,ll),itc,powbase,ERSPtimes,ERSPfreqs] = pop_newtimef(...
            EEG(ii),1,chNum,[tmpEEG.xmin,tmpEEG.xmax]*1e3,0,'freqs',[0,20],...
            'alpha',.05,'plotitc',flags.plotitc,'plotersp',flags.plotersp);
        end
      else

        for ii = 1:length(EEG)

          % Across-trial average of ERPs
          Ntrials(ii,ll) = size(EEG(ii).data,3);
          if length(chNum) == 1
            populationMean{ii}(:,:,ll) = mean(EEG(ii).data(chNum,:,:),3);

            % ERP component amplitudes
            time = kv.epochInterval(1):1/EEGall.srate:kv.epochInterval(2)-1/EEGall.srate;
            for tt = 1:size(kv.latencyRanges,1)
              idTime{tt} = time>=kv.latencyRanges(tt,1) & time<=kv.latencyRanges(tt,2);
%             idTimeP2 = time>=kv.latencyRanges(2,1) & time<=kv.latencyRanges(2,2);
%             idTimeP3 = time>=kv.latencyRanges(3,1) & time<=kv.latencyRanges(3,2);
              Camp{tt}(ii,ll) = min(populationMean{ii}(:,idTime{tt},ll),[],2);
%             P2amp(ii,ll) = max(populationMean{ii}(:,idTimeP2,ll),[],2);
%             P3amp(ii,ll) = max(populationMean{ii}(:,idTimeP3,ll),[],2);
            end
            
            % Peak-to-peak measure
            Camp{4} = Camp{2}-Camp{1};
            kv.latencyLabels{4} = 'P2-N1';
            
          else
            populationMean{ii} = cat(3,populationMean{ii},EEG(ii).data(chNum,:,:));
          end
        end
      end
  end
end

% Display artifact rejection counters
disp(Nrej)

%% Average across listeners
switch flags.contralaterality
  case 'rightHemiDominance'
  otherwise
    
    if flags.do_TFanalysis
        ERSP = mean(ERSP,4);
    else
      for ii = 1:length(EEG)%MLabel)
        populationMean{ii} = mean(populationMean{ii},3);
      end
      if length(chNum) == 1
        for ii = 1:length(EEG)%MLabel)
          populationMean{ii} = squeeze(populationMean{ii});
        end
        if length(IDs) > 1
          N = length(IDs);
          Camp_SE = cell(length(Camp),1);
          for tt=1:length(Camp)
            Camp_SE{tt} = std(Camp{tt},0,2)/sqrt(N);
            Camp{tt} = mean(Camp{tt},2);
          end
%           N1amp_SE = std(N1amp,0,2)/sqrt(N);
%           N1amp = mean(N1amp,2);
%           P2amp_SE = std(P2amp,0,2)/sqrt(N);
%           P2amp = mean(P2amp,2);
%           P3amp_SE = std(P3amp,0,2)/sqrt(N);
%           P3amp = mean(P3amp,2);
        end
      end
      Ntrials = sum(Ntrials,2);
    end
    
end

%% Plots

switch flags.contralaterality
  
  case 'off' % mean epoch
    
    if flags.do_TFanalysis
      
      fig(1) = figure;
      for ff = 1:length(freqBand)
        subplot(1,length(freqBand)+1,ff)
        idf = ERSPfreqs >= freqBand(ff).range(1) & ERSPfreqs <= freqBand(ff).range(2);
        plot([0,0],[-10,10],'k:')
        hold on
        for iD = 1:length(EEG)
          hERSP(iD) = plot(ERSPtimes,squeeze(mean(ERSP(idf,:,iD))));
          set(hERSP(iD),'LineStyle',LineStyle{iD},'Color',Color(iD,:));
        end
        xlabel('Time (ms)')
        ylabel('ERSP power (dB)')
        axis([1000*[-.200,.600],2.5*[-1,1]])
        title([freqBand(ff).name,' (',num2str(freqBand(ff).range(1)),'-',...
          num2str(freqBand(ff).range(2)),' Hz)'])
      end
      legend(hERSP,MlegendLabel,'Location','eastoutside')
      
    elseif flags.do_topo
      
      if length(populationMean) == 2 % include difference
        populationMean{3} = populationMean{2} - populationMean{1};
        MlegendLabel{3} = [MlegendLabel{2},' - ',MlegendLabel{1}];
      end
      
      t = kv.epochInterval(1):1/EEGall.srate:kv.epochInterval(2)-1/EEGall.srate;
%       latencies = [125 210 370];
      
      fig = figure; 
      Nlat = size(kv.latencyRanges,1);
      for ii = 1:length(populationMean)
%         maplimit = max(abs(populationMean{ii}(:)));
        for tt = 1:Nlat
          subplot(length(populationMean),Nlat,tt+(ii-1)*Nlat)
          idt = t >= kv.latencyRanges(tt,1) & t <= kv.latencyRanges(tt,2);
          topoplot(mean(populationMean{ii}(:,idt),2),... 
            EEG(1).chanlocs(chNum),'maplimits',[-4,4]);
          title([MlegendLabel{ii},', ',kv.latencyLabels{tt}])
          if tt==Nlat && ii == length(populationMean)
            colorbar('eastoutside')
          end
        end
      end
    else
      
      fig(1) = figure;
      t = 1000*(kv.epochInterval(1):1/EEGall.srate:kv.epochInterval(2)-1/EEGall.srate);
      plot([0,0],[-10,10],'k:')
      hold on
      for ii = 1:length(EEG)
        h(ii) = plot(t,populationMean{ii});
        set(h(ii),'LineStyle',LineStyle{ii},'Color',Color(ii,:));
      end
      set(gca,'XMinorTick','on','XLim',1e3*kv.erpInterval,'YLim',5.9*[-1,1])
      if length(IDs) == 1
        title(['Listener: ',ID,'; Site: ',flags.site])
      else
        title(['Site: ',flags.site])
      end
      xlabel('Time (ms)')
      ylabel('Amplitude ({\mu}V)')
  %     legend(h,MlegendLabel)
      % legend also with # of averaged trials
      MlegendLabelNtrials = MlegendLabel;
      for ii = 1:length(EEG)
        MlegendLabelNtrials{ii} = [MlegendLabel{ii},' (',num2str(Ntrials(ii)),' epochs)'];
      end
      legend(h,MlegendLabelNtrials,'Location','eastoutside')

      fig(2) = figure;
      for tt=1:length(Camp)
        if length(IDs) == 1
          plot(Camp{tt})
        else
          errorbar(Camp{tt},Camp_SE{tt})
        end
        hold on
      end
      legend(kv.latencyLabels,'Location','eastoutside')
      if length(IDs) == 1
        title(['Listener: ',ID,'; Site: ',flags.site])
      else
        title(['Site: ',flags.site])
      end
      set(gca,'XTick',1:length(EEG),'XTickLabel',MlegendLabel,'YLim',[-6,8])
      ylabel('Amplitude ({\mu}V)')
      if strcmp(flags.timeLock,'Combinations')
        xlabel('D = M_{change} - M_{onset}')
        set(gca,'XTickLabel',{'-1','-2/3','-1/3','0','1/3','2/3','1'})
      end
    end
    
  case 'rightHemiDominance'
    fig = figure;
    switch flags.timeLock
      case 'Onset'
        plot(M,[N1RHdom(:);P2RHdom(:)]')
        xlabel('Spectral contrast, M')
        set(gca,'XTick',sort(M))
        axis([-.2,1.2,-.4,1.6])
      otherwise % 'Combinations
        plot([N1RHdom(:);P2RHdom(:)]')
        XTickLabel = strrep(MlegendLabel,'M:','');
        XTickLabel = strrep(XTickLabel,' ','');
        set(gca,'XTick',1:length(EEG),'XTickLabel',XTickLabel)
        xlabel('Spectral contrast change, M_{onset} -> M_{change}')
        axis([0.5,length(EEG)+.5,-.4,1.6])
    end
    legend('N1','P2','Location','northwest')
    title('Right-hemispheric dominance')
    ylabel('Amplitude difference (\muV)')
end

%% Print
if flags.do_print && exist('fig','var')
  FontSize = 8;
  Resolution = '-r600';
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  set(fig,'PaperUnits','centimeters','PaperPosition',[100,100,12,8])
  
  % common filename
  fn = fullfile('.',mfilename);
  if length(IDs) == 1
    fn = fullfile(fn,IDs{1});
  end
  if not(exist(fn,'dir'))
    mkdir(fn)
  end
  fn = fullfile(fn,[mfilename '_' flags.timeLock '_' flags.site]);
  if not(strcmp(flags.contralaterality,'off'))
    fn = [fn,'_' flags.contralaterality];
  end
  if not(flags.do_all)
    fn = [fn,'_',flags.position];
  end
  if not(flags.do_noGrouping)
    fn = [fn,'_',flags.grouping];
  end
  
  % individual filename 
  if flags.do_TFanalysis
    fn1 = [fn,'_TF'];
    set(fig(1),'PaperUnits','centimeters','PaperPosition',[100,100,20,10])
    print(fig(1),Resolution,'-dpng',fn1)
  elseif flags.do_topo
    set(fig(1),'PaperUnits','centimeters','PaperPosition',[100,100,15,15])
    print(fig(1),Resolution,'-dpng',fn)
  else
    fn1 = [fn,'_ERP'];
    print(fig(1),Resolution,'-dpng',fn1)
  end
  if length(fig) > 1
    fn2 = [fn,'_ERPcompMeasures'];
    print(fig(2),Resolution,'-dpng',fn2)
  end
end