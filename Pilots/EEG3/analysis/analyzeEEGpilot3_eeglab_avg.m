function analyzeEEGpilot3_eeglab_avg(IDs)
% Analyze SpExCue_EEGpilot

flags.do_print = true;

site = {...
%   'frontal';...
%   'central';...
  'parietal';...
  };
timeflag = {...
  'Combinations';...
%   'Onset';...
%   'Change';...
  };
flags.position = {...
%   'left';...
%   'right';...
%   'top';...
  'all';...
  };
flags.contralaterality = {...
  'off';...
%   'ipsiVsContra';...
%   'rightHemiDominance';... % Getzmann & Lewald (2010)
%   'intrahemiContralaterality';... % Palomäki et al. (2005)
  };
flags.do_TFanalysis = false;

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
if flags.do_TFanalysis
  epochStart = -.6;  % seconds before the trigger event
  epochEnd = 1.6; % seconds after the trigger event
else
  epochStart = -.2;  % seconds before the trigger event
  epochEnd = .6; % seconds after the trigger event
end

baselineCorrectInterval = [-200 0];
ERPfhigh = 50; % cut-off frequency of low-pass filter for ERPs

% Settings for automatic epoch rejection
epochThresh = 75; % threshold in µV within entire epoch
baselineThresh = 20; % threshold in µV within baseline interval
maxslope = 5; % max slope in µV/epoch

N1latency = [.120,.220]; % seconds
P2latency = [.220,.300]; % seconds
P3latency = [.300,.600]; % seconds

% Frequency bands for T-F-Analysis
freqBand(1).name = 'delta';
freqBand(1).range = [2,4]; % Hz
freqBand(2).name = 'theta';
freqBand(2).range = [4,8]; % Hz
freqBand(3).name = 'alpha';
freqBand(3).range = [8,12]; % Hz
freqBand(4).name = 'beta';
freqBand(4).range = [15,30]; % Hz

M = [1,1/3,0];
MLabel = {'M1','Mp3','M0'};
LineStyle = {'-','-.',':'};
Color = 0.7*colormap(hsv(4));

dirSupp = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/MATLAB_general/';  % directory with required analysis functions, e.g. 'eeglabTrig2dkrTrigID.m' 
addpath(dirSupp);

% Evaluation site
if strcmp(flags.contralaterality,'off')
  switch site{1}
    case 'frontal'
      chNum = 31; % Fz
    case 'central'
      chNum = 32; % Cz
    case 'parietal'
      chNum = 13; % Pz
  end
else
  chNum = 1:32;
end

for ll = 1:length(IDs)
  
  ID = IDs{ll};

  for pp = 1:length(flags.position)

    %% Load the actual data
    filename = ['SpExCue_EEGpilot3_',ID,'_filt100_ICAclean.set'];
    filepath = '../data/';
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

    for tt = 1%:length(timeflag)

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
          MlegendLabel = {'M_o = 1','M_o = 1/3','M_o = 0'};
          MTrig = [30,20,10];
          for ii = 1:length(M)
            idx = mod(eventList,100) == MTrig(ii);
            eval(['event',MLabel{ii},' = num2cell(unique(eventList(idx)));'])
          end

        case 'Change'
          MlegendLabel = {'M_c = 1','M_c = 1/3','M_c = 0'};
          MTrig = [3,2,1];
          for ii = 1:length(M)
            idx = mod(eventList,10) == MTrig(ii);
            eval(['event',MLabel{ii},' = num2cell(unique(eventList(idx)));'])
          end

        case 'Combinations'
          MLabel = {'Dm1','Dmp7','Dmp3','D0','Dp3','Dp7','D1'};
          MlegendLabel = {'D =-1   ','D =-2/3','D =-1/3','D = 0   ','D = 1/3','D = 2/3','D = 1   '};
          LineStyle = {'-','-','-',':','-.','-.','-.'};
           
          Color = cat(1,Color(1:4,:),flipud(Color(1:3,:)));
          MTrig = {31,32,21,[11,22,33],12,23,13};
          for ii = 1:length(MTrig)
            idx = mod(eventList,100) == MTrig{ii}(1);
            for jj = 2:length(MTrig{ii})
              idx = idx | mod(eventList,100) == MTrig{ii}(jj);
            end
            eval(['event',MLabel{ii},' = num2cell(unique(eventList(idx)));'])
          end

      end

      % Select Epochs and correct for baseline offset
      for ii = 1:length(MLabel)
        eval(['EEG_',MLabel{ii},'= pop_epoch(EEG, event',MLabel{ii},', [epochStart epochEnd]);'])
        eval(['EEG_',MLabel{ii},'= pop_rmbase(EEG_',MLabel{ii},', baselineCorrectInterval);'])
        eval(['EEG_',MLabel{ii},'= eeg_checkset(EEG_',MLabel{ii},');'])
      end

      %% Artifactual epoch rejection
      for ii = 1:length(MLabel)
        % by thresholding
        eval(['EEG_',MLabel{ii},'= pop_eegthresh(EEG_',MLabel{ii},...
          ',1,chNum,-epochThresh,epochThresh,epochStart,epochEnd,0,1);'])
        if epochThresh > baselineThresh
          eval(['EEG_',MLabel{ii},'= pop_eegthresh(EEG_',MLabel{ii},...
            ',1,chNum,-baselineThresh,baselineThresh,',...
            'baselineCorrectInterval(1),baselineCorrectInterval(2),0,1);'])
        end
        % by detecting linear drifts
        eval(['EEG_',MLabel{ii},'= pop_rejtrend(EEG_',MLabel{ii},...
          ',1,chNum,512, maxslope, 0.3,0,1,0);'])
        eval(['EEG_',MLabel{ii},'= eeg_checkset(EEG_',MLabel{ii},');'])
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
            ERSP = [];%nan(12,200,length(MLabel),length(IDs));
            eval(['tmpEEG = EEG_',MLabel{1},';'])
            for ii = 1:length(MLabel)
              eval(['tmpEEG = EEG_',MLabel{ii},';'])
%               figure
              [ERSP(:,:,ii,ll),itc,powbase,ERSPtimes,ERSPfreqs] = pop_newtimef(...
                tmpEEG, 1, chNum,[tmpEEG.xmin,tmpEEG.xmax]*1000,0,...
                'freqs',[0,20],'alpha',.05,'plotitc','on','plotersp','on'); % 
              figure
%               pop_newtimef(EEG_Dm1, 1, chNum,[tmpEEG.xmin,tmpEEG.xmax]*1000,1);
%               eval(['ERPs = squeeze(EEG_',MLabel{ii},'.data(chNum,:,:));'])
%               figure
%               [ERSP(:,:,ii,ll),itc,powbase,ERSPtimes,ERSPfreqs,erspboot,itcboot] = newtimef(...
%                 ERPs,tmpEEG.pnts,[tmpEEG.xmin,tmpEEG.xmax]*1000,tmpEEG.srate,0,...
%                 'freqs',[0,20],'plotitc','off');%,'plotersp','off'
            end
          else

            for ii = 1:length(MLabel)
              % Low-pass filter and baseline-correct the data again
              transitionBandwidth = 1;  % In Hz
              maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
              KaiserWindowBeta = pop_kaiserbeta(maxPassbandRipple);
              filtOrder = pop_firwsord('kaiser', EEG.srate, transitionBandwidth, maxPassbandRipple);
              eval(['tmpEEG = EEG_',MLabel{ii},';'])
              tmpEEG = pop_firws(tmpEEG, 'fcutoff', ERPfhigh, 'ftype', 'lowpass', 'wtype', 'kaiser', 'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);
              tmpEEG = pop_rmbase(tmpEEG, baselineCorrectInterval);
              eval(['EEG_',MLabel{ii},' = tmpEEG;'])

              % Across-trial average of ERPs
              Ntrials(ii,ll) = size(tmpEEG.data,3);
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
end

%% Average across listeners
switch flags.contralaterality{end}
  case 'ipsiVsContra'
  case {'rightHemiDominance','intrahemiContralaterality'}
  otherwise
    
    if flags.do_TFanalysis
        ERSP = mean(ERSP,4);
    else
      for ii = 1:length(MLabel)
        eval(['populationMean',MLabel{ii},' = mean(populationMean',MLabel{ii},',2);']);
      end
      N1amp = mean(N1amp,3);
      P2amp = mean(P2amp,3);
      P3amp = mean(P3amp,3);
      Ntrials = sum(Ntrials,2);
    end
    
end

%% Plots

switch flags.contralaterality{end}
  
  case 'off' % mean epoch
    
    if flags.do_TFanalysis
      
      fig(1) = figure;
      for ff = 1:length(freqBand)
        subplot(1,length(freqBand)+1,ff)
        idf = ERSPfreqs >= freqBand(ff).range(1) & ERSPfreqs <= freqBand(ff).range(2);
        plot([0,0],[-10,10],'k:')
        hold on
        for iD = 1:length(MLabel)
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
      
    else
      
      fig(1) = figure;
      t = 1000*(epochStart:1/EEG.srate:epochEnd-1/EEG.srate);
      plot([0,0],[-10,10],'k:')
      hold on
      for ii = 1:length(MLabel)
        eval(['h(ii) = plot(t,populationMean',MLabel{ii},');']);
        set(h(ii),'LineStyle',LineStyle{ii},'Color',Color(ii,:));
      end
      set(gca,'XMinorTick','on','XLim',1000*[epochStart,epochEnd],'YLim',5.9*[-1,1])
      if length(IDs) == 1
        title(['Listener: ',ID,'; Site: ',site{1}])
      else
        title(['Site: ',site{1}])
      end
      xlabel('Time (ms)')
      ylabel('Amplitude ({\mu}V)')
  %     legend(h,MlegendLabel)
      % legend also with # of averaged trials
      MlegendLabelNtrials = MlegendLabel;
      for ii = 1:length(MLabel)
        MlegendLabelNtrials{ii} = [MlegendLabel{ii},' (',num2str(Ntrials(ii)),' epochs)'];
      end
      legend(h,MlegendLabelNtrials,'Location','eastoutside')

      fig(2) = figure;
      plot([N1amp(pp,:);P2amp(pp,:);P2amp(pp,:)-N1amp(pp,:);P3amp(pp,:)]')
      legend('N1','P2','P2-N1','P3','Location','eastoutside')
      if length(IDs) == 1
        title(['Listener: ',ID,'; Site: ',site{1}])
      else
        title(['Site: ',site{1}])
      end
      set(gca,'XTick',1:length(MLabel),'XTickLabel',MlegendLabel,'YLim',[-6,8])
      ylabel('Amplitude ({\mu}V)')
      if strcmp(timeflag{1},'Combinations')
        xlabel('D = M_{change} - M_{onset}')
        set(gca,'XTickLabel',{'-1','-2/3','-1/3','0','1/3','2/3','1'})
      end
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
      legend(h,'Right hemisphere','Left hemisphere','Location','westoutside')
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
  
  % common filename
  fn = fullfile('.',mfilename);
  if length(IDs) == 1
    fn = fullfile(fn,IDs{1});
  end
  if not(exist(fn,'dir'))
    mkdir(fn)
  end
  fn = fullfile(fn,[mfilename '_' timeflag{tt} '_' site{1}]);
  if not(strcmp(flags.contralaterality{end},'off'))
    fn = [fn,'_' flags.contralaterality{end}];
  end
  fn = [fn,'_',flags.position{pp}];
  
  % individual filename 
  if flags.do_TFanalysis
    fn1 = [fn,'_TF'];
    set(fig(1),'PaperUnits','centimeters','PaperPosition',[100,100,20,10])
    print(fig(1),Resolution,'-dpng',fn1)
  else
    fn1 = [fn,'_ERP'];
    print(fig(1),Resolution,'-dpng',fn1)
  end
  if length(fig) > 1
    fn2 = [fn,'_ERPcompMeasures'];
    print(fig(2),Resolution,'-dpng',fn2)
  end
end