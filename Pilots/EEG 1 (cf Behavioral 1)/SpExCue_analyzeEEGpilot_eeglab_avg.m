function SpExCue_analyzeEEGpilot_eeglab_avg(IDs)
% Analyze SpExCue_EEGpilot

flags.do_print = false;

timeflag = {...
  'Combinations';...
%   'Onset';...
%   'Change';...
  };
flags.position = {...
  'left';...
%   'right';...
%   'top';...
%   'all';...
  };
flags.contralaterality = {...
  'off';...
%   'ipsiVsContra';...
%   'rightHemiDominance';... % Getzmann & Lewald (2010)
%   'intrahemiContralaterality';... % Palomäki et al. (2005)
  };
flags.do_TFanalysis = true;

%% 

if not(exist('IDs','var'))
  IDs = {input('Subject ID: ','s')};
elseif ischar(IDs)
  IDs = {IDs};
end

if not(exist('pop_loadset','file'))
  addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab/')
  eeglab
end

%% Settings

epochStart = -.2;  % seconds before the trigger event
epochEnd = .8; % seconds after the trigger event
baselineCorrectInterval = [-100 0];
N1latency = [.120,.220]; % seconds
P2latency = [.220,.300]; % seconds
% P2latency = [.250,.600]; % seconds
P3latency = [.300,.600]; % seconds

M = [1,0.5,0];
MLabel = {'M1','Mp5','M0'};
LineStyle = {'-','-.',':'};

chNum = 13; % Cz: 32; Pz: 13

thresh = 70;

dirSupp = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/MATLAB_general/';  % directory with required analysis functions, e.g. 'eeglabTrig2dkrTrigID.m' 
addpath(dirSupp);

for ll = 1:length(IDs)
  
  ID = IDs{ll};

for pp = 1:length(flags.position)

%% Load the actual data
filename = ['SpExCue_EEGpilot_',ID,'_filt20_ICAclean.set'];
filepath = './Results/';
fnLOCS = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Tools/EEG/biosemi_eloc.locs';
EEG = pop_loadset('filename', filename, 'filepath', filepath); 
EEG = eeg_checkset(EEG);  % verify consistency of the fields of an EEG dataset
EEG = pop_chanedit(EEG, 'load',{fnLOCS 'filetype' 'autodetect'});

%% Convert and adjust event list
% Take only last 8 bits
eventList = eeglabTrig2dkrTrigID([EEG.event.type]); 
% Adjust to code M combination in event list (bug in first run of EEGpilot)
idx = floor(mod(eventList,100)/10) == 0 & mod(eventList,10) < 8 & mod(eventList,10) > 0;
eventList(idx) = mod(eventList(circshift(idx,-1)),100) + eventList(idx);
% Write into EEG object
for i = 1:length(eventList)
    EEG.event(i).type = eventList(i);
end

for tt = 1:length(timeflag)
  
%% Event selection

switch flags.position{pp} % set undesired events to 0
  case 'left' % 0
    eventList(eventList>=100) = 0;
  case 'top' % 100
    eventList(eventList>=200 | eventList<100) = 0;
  case 'right' % 200
    eventList(eventList<200) = 0;
%   otherwise % all
end

switch timeflag{tt}
  
  case 'Onset'
    MlegendLabel = {'M_o = 1','M_o = 0.5','M_o = 0'};
    MTrig = [30,20,10];
    for ii = 1:length(M)
      idx = mod(eventList,100) == MTrig(ii);
      eval(['event',MLabel{ii},' = num2cell(unique(eventList(idx)));'])
    end
    
  case 'Change'
    MlegendLabel = {'M_c = 1','M_c = 0.5','M_c = 0'};
    MTrig = [3,2,1];
    for ii = 1:length(M)
      idx = mod(eventList,10) == MTrig(ii);
      eval(['event',MLabel{ii},' = num2cell(unique(eventList(idx)));'])
    end
    
  case 'Combinations'
    MLabel = {'M1_0','M1_p5','Mp5_0','Mp5_1','M0_p5','M0_1'};
    MlegendLabel = {'M: 1  -> 0 ','M: 1  -> .5','M: .5 -> 0 ','M: .5 -> 1 ','M: 0  -> .5','M: 0  -> 1 '};
    LineStyle = {'-','-','-.','-.',':',':'};
    MTrig = [31,32,21,23,12,13];
    for ii = 1:length(MTrig)
      idx = mod(eventList,100) == MTrig(ii);
      eval(['event',MLabel{ii},' = num2cell(unique(eventList(idx)));'])
    end
    
end

% Select Epochs and correct for baseline offset
for ii = 1:length(MLabel)
  eval(['EEG_',MLabel{ii},'= pop_epoch(EEG, event',MLabel{ii},', [epochStart epochEnd]);'])
  eval(['EEG_',MLabel{ii},'= pop_rmbase(EEG_',MLabel{ii},', baselineCorrectInterval);'])
  eval(['EEG_',MLabel{ii},'= eeg_checkset(EEG_',MLabel{ii},');'])
end

%% Thresholding
for ii = 1:length(MLabel)
  eval(['MaxP = squeeze(max(max(abs(EEG_',MLabel{ii},'.data),[],1),[],2));'])
  idsel{ii} = find(MaxP < thresh);
  eval(['EEG_',MLabel{ii},'.data = EEG_',MLabel{ii},'.data(:,:,idsel{ii});'])
end

%% Average epochs and evaluate component measures
switch flags.contralaterality{end}
  case 'ipsiVsContra'
    for ii = 1:length(MLabel)
      eval(['tmp.type = cat(1,EEG_',MLabel{ii},'.event.type);'])
      eval(['tmp.epoch = cat(1,EEG_',MLabel{ii},'.event.epoch);'])
      eval(['idleft = tmp.epoch(tmp.type == event',MLabel{ii},'{1});'])
      eval(['idright = tmp.epoch(tmp.type == event',MLabel{ii},'{3});'])
      idleft = ismember(idsel{ii},idleft);
      idright = ismember(idsel{ii},idright);
      eval(['popMeanLeft = mean(EEG_',MLabel{ii},'.data(:,:,idleft),3);'])
      eval(['popMeanRight = mean(EEG_',MLabel{ii},'.data(:,:,idright),3);'])
      eval(['populationMean',MLabel{ii},' = popMeanLeft - popMeanRight;']);
    end
  case {'rightHemiDominance','intrahemiContralaterality'}
    idLH = ismember({EEG.chanlocs.labels},... % left hemisphere
      {'AF3','AF7','C1','C3','C5','CP1','CP3','CP5','F3','F7','FC1','FC3',...
      'FC5','Fp1','FT9','O1','P1','P3','P7','PO3','PO7','PO9','T7','TP7'});
    idRH = ismember({EEG.chanlocs.labels},... % right hemisphere
      {'AF4','AF8','C2','C4','C6','CP2','CP4','CP6','F4','F8','FC2','FC4',...
      'FC6','Fp2','FT10','O2','P2','P4','P8','PO10','PO4','PO8','T8','TP8'});
    time = epochStart:1/EEG.srate:epochEnd-1/EEG.srate;
    idTimeN1 = time>=N1latency(1) & time<=N1latency(2);
    idTimeP2 = time>=P2latency(1) & time<=P2latency(2);
    for ii = 1:length(MLabel)
      eval(['N1amp = min(mean(EEG_',MLabel{ii},'.data(:,idTimeN1,:),3),[],2);'])
      eval(['P2amp = max(mean(EEG_',MLabel{ii},'.data(:,idTimeP2,:),3),[],2);'])
      N1RH(pp,ii) = rms(N1amp(idRH));
      N1LH(pp,ii) = rms(N1amp(idLH));
      P2RH(pp,ii) = rms(P2amp(idRH));
      P2LH(pp,ii) = rms(P2amp(idLH));
    end
    N1RHdom = N1RH - N1LH; % right-hemispheric dominance
    P2RHdom = P2RH - P2LH; % right-hemispheric dominance
  otherwise
    
    if flags.do_TFanalysis
      % Downsampling parameters
      fsLow = 40;
      fsDiv = gcd(fsLow,EEG.srate);
      P = fsLow/fsDiv;
      Q = EEG.srate/fsDiv;
      % DGT parameters
      WinFunct={'gauss',1,'2'}; % Gaussian window with time-freq ratio of 1 and normalized L2 norm
      DownSampFact=1; % Downsampling factor in time
      NumChan=40; % Total number of channels, only 101 will be computed
      ersp = nan(12,200,length(MLabel),length(IDs));
    end
    
    for ii = 1:length(MLabel)
      
      % T-F characteristics
      if flags.do_TFanalysis
        eval(['ERPs = squeeze(EEG_',MLabel{ii},'.data(chNum,:,:));'])
        [ersp(:,:,ii,ll),itc,powbase,times,freqs,erspboot,itcboot] = newtimef(...
          ERPs,EEG_M0_1.pnts,[EEG_M0_1.xmin,EEG_M0_1.xmax]*1000, EEG.srate); % 'freqs',[1,20]
        
%         % Downsampling
%         ERPs_fsLow = resample(double(ERPs),P,Q);
%         % DGT
%         PSD = db(abs(dgtreal(ERPs_fsLow,WinFunct,DownSampFact,NumChan)));
%         eval(['PSD_',MLabel{ii},'(:,:,ll) = mean(PSD,3);'])
      end
      
      % Across-trial average of ERPs
      eval(['populationMean',MLabel{ii},'(:,ll) = squeeze(mean(EEG_',MLabel{ii},'.data(chNum,:,:),3));']);
      % ERP component amplitudes
      time = epochStart:1/EEG.srate:epochEnd-1/EEG.srate;
      idTimeN1 = time>=N1latency(1) & time<=N1latency(2);
      idTimeP2 = time>=P2latency(1) & time<=P2latency(2);
      idTimeP3 = time>=P3latency(1) & time<=P3latency(2);
      eval(['N1amp(pp,ii,ll) = min(populationMean',MLabel{ii},'(idTimeN1,ll),[],1);'])
      eval(['P2amp(pp,ii,ll) = max(populationMean',MLabel{ii},'(idTimeP2,ll),[],1);'])
      eval(['P3amp(pp,ii,ll) = max(populationMean',MLabel{ii},'(idTimeP3,ll),[],1);'])
    end
end

end
end
end

%% Average across listeners
switch flags.contralaterality{end}
  case 'ipsiVsContra'
  case {'rightHemiDominance','intrahemiContralaterality'}
  otherwise
    for ii = 1:length(MLabel)
      eval(['populationMean',MLabel{ii},' = mean(populationMean',MLabel{ii},',2);']);
      if flags.do_TFanalysis
        ersp = mean(ersp,4);
      end
    end
    N1amp = mean(N1amp,3);
    P2amp = mean(P2amp,3);
    P3amp = mean(P3amp,3);
end

%% Plots

switch flags.contralaterality{end}
  case 'off' % mean epoch
    
    fig(1) = figure;
    t = 1000*(epochStart:1/EEG.srate:epochEnd-1/EEG.srate);
    plot([0,0],[-10,10],'k:')
    hold on
    for ii = 1:length(MLabel)
      eval(['h(ii) = plot(t,populationMean',MLabel{ii},')']);
      set(h(ii),'LineStyle',LineStyle{ii});
    end
    set(gca,'XMinorTick','on','XLim',1000*[epochStart,epochEnd],'YLim',4.9*[-1,1])
    title([timeflag{tt},', listener: ',ID])
    xlabel('Time (ms)')
    ylabel('Amplitude ({\mu}V)')
    legend(h,MlegendLabel)
    
    fig(2) = figure;
    plot([N1amp(pp,:);P2amp(pp,:);P2amp(pp,:)-N1amp(pp,:);P3amp(pp,:)]')
    legend('N1','P2','P2-N1','P3')
    title([timeflag{tt},', listener: ',ID])
    set(gca,'XTick',1:length(MLabel),'XTickLabel',MlegendLabel)
    ylabel('Amplitude ({\mu}V)')
    
    if flags.do_TFanalysis
      fig(3) = figure;
      for ff = 1:4
        subplot(2,2,ff)
        plot(times,squeeze(ersp(ff,:,:)))
        title([num2str(freqs(ff)),' Hz'])
        legend(MlegendLabel)
      end
%       dynrange = 20;
%       for ii = 1:length(MLabel)
%         subplot(2,3,ii)
%         eval(['plotdgtreal(PSD_',MLabel{ii},',DownSampFact,NumChan,fsLow,dynrange);']);
%         set(gca,'XTickLabel',get(gca,'XTick')+baselineCorrectInterval(1)/1e3)
%         title(MlegendLabel{ii})
%       end
    end
    
  case 'rightHemiDominance'
    fig = figure;
    switch timeflag{tt}
      case 'Onset'
        plot(M,[N1RHdom(pp,:);P2RHdom(pp,:)]')
        xlabel('Spectral contrast, M')
        set(gca,'XTick',sort(M))
        axis([-.2,1.2,-.4,1.6])
      otherwise % 'Combinations
        plot([N1RHdom(pp,:);P2RHdom(pp,:)]')
        XTickLabel = strrep(MlegendLabel,'M:','');
        XTickLabel = strrep(XTickLabel,' ','');
        set(gca,'XTick',1:length(MLabel),'XTickLabel',XTickLabel)
        xlabel('Spectral contrast change, M_{onset} -> M_{change}')
        axis([0.5,length(MLabel)+.5,-.4,1.6])
    end
    legend('N1','P2','Location','northwest')
    title(['Right-hemispheric dominance (',flags.position{pp},')'])
    ylabel('Amplitude difference (\muV)')
  case 'intrahemiContralaterality'
    if pp == length(flags.position)
      fig = figure;
      idL = ismember(flags.position,'left');
      idR = ismember(flags.position,'right');
      N1RHdiff = N1RH(idL,:) - N1RH(idR,:); % contra-ipsilateral stimulus 
      N1LHdiff = N1LH(idR,:) - N1LH(idL,:); % contra-ipsilateral stimulus 
      plot([-10,10],zeros(2,1),'k-')
      hold on
      switch timeflag{tt}
        case 'Onset'
          h = plot(M,[N1RHdiff;N1LHdiff]');
          xlabel('Spectral contrast, M')
          set(gca,'XTick',sort(M),'XLim',[-.2,1.2])
        otherwise % 'Combinations
          h = plot([N1RHdiff;N1LHdiff]');
          xlabel('Spectral contrast change, M_{onset} -> M_{change}')
          XTickLabel = strrep(MlegendLabel,'M:','');
          XTickLabel = strrep(XTickLabel,' ','');
          set(gca,'XTick',1:length(MLabel),'XTickLabel',XTickLabel,'XLim',[.8,length(MLabel)+.2])
      end
      legend(h,'Right hemisphere','Left hemisphere','Location','northwest')
      ylabel('N1 amplitude difference (\muV)')
      title('Contra- vs. ipsilateral stimulation')
    end
end

%% Print
if flags.do_print && exist('fig','var')
  FontSize = 8;
  Resolution = '-r600';
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  set(fig,'PaperUnits','centimeters','PaperPosition',[100,100,12,8])
  fn = fullfile('.',mfilename);
  if not(exist(fn,'dir'))
    mkdir(fn)
  end
  fn = fullfile(fn,[mfilename '_' timeflag{tt}]);
  if not(strcmp(flags.contralaterality{end},'off'))
    fn = [fn,'_' flags.contralaterality{end}];
  end
  fn = [fn,'_',flags.position{pp}];
  print(fig(1),Resolution,'-dpng',fn)
  if length(fig) > 1
    fn2 = [fn,'_compMeasures'];
    print(fig(2),Resolution,'-dpng',fn2)
  end
  if length(fig) > 2
    fn3 = [fn,'_sgram'];
    set(fig(3),'PaperUnits','centimeters','PaperPosition',[100,100,16,10])
    print(fig(3),Resolution,'-dpng',fn3)
  end
end