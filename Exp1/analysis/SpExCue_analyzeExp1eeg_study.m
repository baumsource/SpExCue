clear all

%% Settings

flags.do_blinkICrej = true;
conditionSet = 'Mchange'; % 'CloserVsFarther' or 'Mchange'
% designs for 'Mchange': 
% 1: all, 
% 2: M1-0vsM0-1, 3: M1-ivsMi-1, 4: Mi-0vsM0-i, 5: M1-ivsM0-i, 6: D<0vsD>0
% 7: M1-0vsMi-0, 8: M0-1vsMi-1, 9: M0-ivsM1-i
design = 6;

% Precomp settings
precomp.redo = 'off';
precomp.data = 'channels'; % 'components' or 'channels'
precomp.erp = 'on';
precomp.scalp = 'on';
precomp.spec = 'off';
precomp.ersp = 'off'; % 'on'
erpparams = { 'recompute',precomp.redo };
specparams = { 'freqrange', [2,20], 'timerange', [-200 800],'recompute',precomp.redo };
erspparams = { 'cycles', [2,10], 'freqs', [2 20],'nfreqs',10,...
  'alpha', 0.01, 'padratio' 1,'timelimits', [-200 800],'recompute',precomp.redo};

% Clustering
flags.do_cluster = strcmp(precomp.data,'components'); % always cluster in case of component-based analysis
flags.do_channels = not(flags.do_cluster);
preclustsetting = {...
  {'erp' 'npca' 5 'norm' 1 'weight' 1 'timewindow' [-200 800]},...
  {'scalp' 'npca' 5 'norm' 1 'weight' 1},...
  {'spec' 'npca' 5 'norm' 1 'weight' 1 'freqrange' [2 20]}};
  % 	{'ersp' 'npca' 10 'norm' 1 'weight' 1}
  % 	{'dipoles' 'norm' 1 'weight' 10} , ... 
clus_num = 12; % number of clusters

% Statistical analysis
stats.mcorrect = 'none'; % 'none' or 'MonteCarloCluster'

% Plots
plotflag.topo = false;
plotflag.erp = true;
plotflag.corr = true;
flags.do_print = true;

% Data
filepath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/data';

%% Subjects and data files

if flags.do_blinkICrej
  fnX = 'Exp1eeg_Sxx_blinkICrej_Cxx.set';
else
  fnX = 'Exp1eeg_Sxx_ICAraw_Cxx.set';
end
tmp = load('SpExCue_Exp1eeg_subjects');
subject = tmp.subjects; %{'S01','S02','S04','S07','S08','S09','S11'}; % 'S01','S02','S04','S07','S08',  ;RS: M = 0.5; S03: reversed externalization judgements


%% Start eeglab

if not(exist('pop_loadset','file'))
  addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab/')
  eeglab
end

%% Create/load STUDY

fnStudy = strrep(fnX,'Sxx','ALL');
fnStudy = strrep(fnStudy,'Cxx',conditionSet);
fnStudy = strrep(fnStudy,'set','study');
  
switch conditionSet
  case 'Mchange'
    condition = {'M1-0','M1-i','Mi-0','M0-0','Mi-i','M1-1','M0-i','Mi-1','M0-1'};
  case 'CloserVsFarther'
    condition = {'Closer','Farther'};
  otherwise
    error('RB: Check conditionSet!')
end

if exist(fullfile(filepath,fnStudy),'file')
  
  [STUDY ALLEEG] = pop_loadstudy('filename', fnStudy,'filepath',filepath);
  STUDY.currentdesign = design;
  
else

  commands = cell(1,length(subject)*length(condition));
  fn = commands;
  ii = 0;
  for ss = 1:length(subject)
    for cc = 1:length(condition)
      ii = ii+1;
      fn{ii} = strrep(fnX,'Sxx',subject{ss});
      fn{ii} = strrep(fn{ii},'Cxx',condition{cc});
      commands{ii} = { 'index' ii 'load' fullfile(filepath,fn{ii}),...
        'subject' subject{ss} 'condition' condition{cc} };
    end
  end

  % Save
  [STUDY,ALLEEG] = std_editset([],[],'commands',commands,'name','SpExCue',...
    'task', 'closer vs. farther for second of two consecutive noise bursts'); 
  STUDY.design(1).name = 'all';

  if strcmp(conditionSet,'Mchange')
    % Create additional designs
    STUDY = std_makedesign(STUDY, ALLEEG, 2, 'name', 'M1-0vsM0-1',...
      'variable1', 'condition', 'values1',{'M1-0','M0-1'},'subjselect',subject);
    STUDY = std_makedesign(STUDY, ALLEEG, 3, 'name', 'M1-ivsMi-1',...
      'variable1', 'condition', 'values1',{'M1-i','Mi-1'},'subjselect',subject);
    STUDY = std_makedesign(STUDY, ALLEEG, 4, 'name', 'Mi-0vsM0-i',...
      'variable1', 'condition', 'values1',{'Mi-0','M0-i'},'subjselect',subject);
    STUDY = std_makedesign(STUDY, ALLEEG, 5, 'name', 'M1-ivsM0-i',...
      'variable1', 'condition', 'values1',{'M1-i','M0-i'},'subjselect',subject);
    STUDY = std_makedesign(STUDY, ALLEEG, 6, 'name', 'D<0vsD>0','subjselect',subject,...
      'variable1', 'condition', 'values1',{{'M1-i','M1-0','Mi-0'},{'Mi-1','M0-1','M0-i'}});
    STUDY = std_makedesign(STUDY, ALLEEG, 7, 'name', 'M1-0vsMi-0',...
      'variable1', 'condition', 'values1',{'M1-0','Mi-0'},'subjselect',subject);
    STUDY = std_makedesign(STUDY, ALLEEG, 8, 'name', 'M0-1vsMi-1',...
      'variable1', 'condition', 'values1',{'M0-1','Mi-1'},'subjselect',subject);
    STUDY = std_makedesign(STUDY, ALLEEG, 9, 'name', 'M0-ivsM1-i',...
      'variable1', 'condition', 'values1',{'M0-i','M1-i'},'subjselect',subject); 
    STUDY.currentdesign = design;
  end

  % Save
  [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG,'filename', fnStudy,'filepath',filepath);
  
end

% update workspace variables and redraw EEGLAB window
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
eeglab redraw

%% Compute measures

% Number of trials per subject and condition
Ntrials = zeros(length(condition),length(subject)+1);
for ii = 1:length(STUDY.design(design).cell)
  ss = ismember(subject,STUDY.design(design).cell(ii).case);
  cc = ismember(condition,STUDY.design(design).cell(ii).value{1});
  Ntrials(cc,ss) = length(STUDY.design(design).cell(ii).trials{1});
end
Ntrials(:,end) = sum(Ntrials(:,1:length(subject)),2);
NtrialsTab = array2table(Ntrials,'RowNames',condition,'VariableNames',[subject,{'All'}]);

% Precomputation
precompsettings = {'interp', 'on', 'recompute', precomp.redo,...
  'scalp', precomp.scalp, ...
  'erp', precomp.erp, 'erpparams', erpparams,...
  'spec', precomp.spec, 'specparams', specparams,...
  'ersp', precomp.ersp, 'itc', precomp.ersp, 'erspparams', erspparams};
[STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, precomp.data,precompsettings{:});

if flags.do_cluster
  % Clusstering based on scalp map correlation to a specified template
  % (here: EEG set 1, IC 6)
%   [correlationMap,STUDY,ALLEEG] = corrmap(STUDY,ALLEEG,1,6,'th','.75');
  
  % Preclustering
  [STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, [], preclustsetting{:});

  % Clustering
  [STUDY ALLEEG] = pop_clust(STUDY, ALLEEG,'algorithm','kmeans',...
    'clus_num',clus_num,'outliers',3);
end

% Channel labels
ChanLbls = {ALLEEG(1).chanlocs.labels};
              
%% Set statistics
switch stats.mcorrect
  case 'none'
    STUDY = pop_statparams(STUDY,'condstats','on','mode','eeglab','mcorrect','none','alpha',nan);
  case 'MonteCarloCluster'
    STUDY = pop_statparams(STUDY,'condstats','on','mode','fieldtrip',...
      'fieldtripmethod','montecarlo','fieldtripmcorrect','cluster');
end

%% Explore clusters
if flags.do_cluster
  pop_clustedit(STUDY,ALLEEG);
  pause
end

%% ERP analysis
if flags.do_channels
  
  % grand average ERP: evaluate site with maximum deflection -> Cz 
%     for jj = 1:32; 
%       [STUDY, erpdata, erptimes, ~, pcond] = std_erpplot(...
%         STUDY,ALLEEG,'channels',ChanLbls(jj),'noplot','on'); 
%       ERPmax(jj) = max(max(abs(mean(erpdata{1},2)),abs(mean(erpdata{2},2)))); 
%     end
%     [~,iERPmax] = max(ERPmax);
  
  iERPmax = 32; % Cz (evaluated by code above)
  chanLabel = ChanLbls(iERPmax);
  peakTimes = [120,220]; % extracted from Cz grand average ERP
  timeRange = 50;
  compLabel = {'N1','P2'};
  
  figERP = [];
  if plotflag.erp
    STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
      'timerange',[-200,800],'averagechan','on','topotime',[]);
    
    [STUDY,erpdata,erptimes,~,pcond] = std_erpplot(STUDY,ALLEEG,'channels',chanLabel);
    figERP = gcf;

    % BUG: EEGlab does not plot wave forms for fixed alpha in statparams ->
    % workaround with alpha = nan in statparams:
    alpha = .05;
    tsigdiff = erptimes(pcond{1} < alpha); % sig. diff. time instances 
    delete(subplot(122));
    pos = get(subplot(121),'Position');
    set(subplot(121),'Position',pos.*[1,1,2,1])
    YLim = get(gca,'YLim');
    hold on
    plot(tsigdiff,(YLim(1)+.1*diff(YLim))*ones(length(tsigdiff),1),'r.')
  end
  
  figTopo = [];
  if plotflag.topo
    STUDY = pop_erpparams(STUDY, 'plotconditions', 'apart','averagechan','off');
    for ii = 1:length(peakTimes)
      topotime = peakTimes(ii) + timeRange*[-1,1];
      [STUDY, C(ii).erpdata, ~, ~, C(ii).pcond] = std_erpplot(STUDY,ALLEEG,...
        'channels',ChanLbls(1:32), 'topotime',topotime,'ylim',[-2,2]);
      figTopo(ii) = gcf;
    end
%     topotimeLabel = [num2str(tsigdiff(1),'%u'),'-',num2str(tsigdiff(end),'%u'),'ms'];
%     topotime = [tsigdiff(1), tsigdiff(end)];
%     [STUDY, erpdata, ~, ~, pcond] = std_erpplot(STUDY,ALLEEG,...
%         'channels',ChanLbls(1:32), 'topotime',topotime,'ylim',[-2,2]);
%     figTopo = gcf;
  end
  
  figCorr = [];
  if plotflag.corr
    
    [STUDY,erpdata,erptimes,~,pcond] = std_erpplot(STUDY,ALLEEG,...
      'channels',ChanLbls(iERPmax),'noplot','on');
    
    behav = load('SpExCue_analyzeExp1behav_avg_eeg.mat');
      
    compAmp = nan(length(erpdata),length(STUDY.subject),length(peakTimes));
    for ii = 1:length(peakTimes) % N1 and P2
      
      timeBounds = peakTimes(ii)+timeRange*[-1,1]; % ms
      iComp = erptimes >= timeBounds(1) & erptimes <= timeBounds(2);
      for jj=1:length(erpdata)
        compAmp(jj,:,ii) = mean(erpdata{jj}(iComp,:));
      end
      
%       x = -diff(compAmp(:,:,ii))./mean(compAmp(:,:,ii)); % amplitude difference normalized by mean amplitude
% %       y = nanmean(behav.dprime);
%       y = nanmean(behav.bias);
%       [B,BINT,R,RINT,STATS] = regress(y(:),[ones(length(x),1),x(:)]);
%       r = sqrt(STATS(1));
%       p = STATS(3);
%       
%       figCorr(ii) = figure;
%       symb = '<s>';
%       legLbl = {'left','front','right'};
%       for pp = 1:size(behav.dprime,1) % directions
%         idir = not(isnan(behav.dprime(pp,:)));
%         h(pp) = plot(x(idir),y(idir),['k',symb(pp)]);
%         hold on
%       end
%       set(h,'MarkerFaceColor','k')
%       yRegress = B(1) + B(2)*x;
%       plot(x,yRegress,'k')
%       title([STUDY.design(design).name,', r = ',num2str(r,'%1.2f'),', p = ',num2str(p,'%1.2f')])
%       ylabel('bias towards closer (c)')
%       xlabel('Normalized P2 amplitude difference (µV)')
%       legend(legLbl,'Location','northwest')
      
    end
    
%     N1range = 120+50*[-1,1]; % ms
%     P2range = 220+50*[-1,1]; % ms
%     Xrange = [140,190]; % ms
%     iN1 = erptimes >= N1range(1) & erptimes <= N1range(2);
%     iX = erptimes >= Xrange(1) & erptimes <= Xrange(2);
%     iP2 = erptimes >= P2range(1) & erptimes <= P2range(2);
%     for ii=1:length(erpdata)
%       N1(ii,:) = mean(erpdata{ii}(iN1,:));
%       X(ii,:) = mean(erpdata{ii}(iX,:));
%       P2(ii,:) = mean(erpdata{ii}(iP2,:));
%     end
%     R=0:10:100; for jj=1:length(R); P2range = 220+R(jj)*[-1,1]; iP2 = erptimes >= P2range(1) & erptimes <= P2range(2); for ii=1:length(erpdata); P2(ii,:) = mean(erpdata{ii}(iP2,:)); end; [h(jj),p(jj)]=ttest(diff(P2)); end
%     [h,p] = ttest(diff(P2))

    
%     [r,p] = corrcoef(-diff(P2)./mean(P2),nanmean(behav.dprime));
%     figure; scatter(-diff(P2)./mean(P2),nanmean(behav.dprime)); 
%     title([STUDY.design(design).name,', r = ',num2str(r(2),'%1.2f'),', p = ',num2str(p(2),'%1.2f')])
%     ylabel('dprime')
%     xlabel('Normalized P2 amplitude difference (µV)')
  end
  
end

%% IC cluster analysis

if flags.do_cluster

  % plot topographies
  std_topoplot(STUDY,ALLEEG,'clusters',3:length(STUDY.cluster));
  figTopo = gcf;

  % plot ERPs
  if plotflag.erp
    STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
      'timerange',[-200,800],'averagechan','on','topotime',[]);
    alpha = .05;
    for jj=3:length(STUDY.cluster)
      [STUDY, erpdata, erptimes, ~, pcond] = std_erpplot(...
            STUDY,ALLEEG,'clusters',jj);
    end
  end

end

%% Topographic analysis 
% STUDY = pop_erpparams(STUDY, 'plotconditions', 'apart',...
%     'timerange',[-200,800],'averagechan','off','topotime',[]);
% figTopo = [];
% if plotflag.topo
%   if not(flags.do_cluster)
%     topotime = {[100 150];[180,200];[200 250];[300 600]};
%     topotimeLabel = {'N1','P2a','P2','P3'};
%     for ii = 1:length(topotime)
%       [STUDY, C(ii).erpdata, ~, ~, C(ii).pcond] = std_erpplot(STUDY,ALLEEG,...
%         'channels',ChanLbls(1:32), 'topotime',topotime{ii},'ylim',[-2,2]);
%       figTopo(ii) = gcf;
%     end
%   else % Cluster topographies
%     for ii=1:length(STUDY.cluster)
%       std_topoplot(STUDY,ALLEEG,'clusters',ii); 
%       figTopo(ii) = gcf;
%     end
%   end
% end

%% ERP analysis

% figERP = [];
% if plotflag.erp
%   STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
%     'timerange',[-200,800],'averagechan','on','topotime',[]);
%   alpha = .05;
%   if not(flags.do_cluster)
%     % evaluate site with maximum deflection of grand average ERP
% %     for jj = 1:32; 
% %       [STUDY, erpdata, erptimes, ~, pcond] = std_erpplot(...
% %         STUDY,ALLEEG,'channels',ChanLbls(jj),'noplot','on'); 
% %       MaxERP(jj) = max(max(abs(mean(erpdata{1},2)),abs(mean(erpdata{2},2)))); 
% %     end
% %     [~,N] = max(MaxERP);
%     N = 32; % Cz (evaluated by code above)
%     lbl = ChanLbls(N);
%     % significant topography cluster average:
% %     N=1;
% %     jj=C(ii).pcond{1}<.05;
%   else
%     N = 3:length(STUDY.cluster);
%     lbl = {STUDY.cluster.name};
%   end
%   ChLAB = [];    
%   kk = 1; % counter
%   for jj = N
%     if not(flags.do_cluster)
%       [STUDY, erpdata, erptimes, ~, pcond] = std_erpplot(...
%         STUDY,ALLEEG,'channels',ChanLbls(jj));
%     else
%       if length(STUDY.cluster(jj).comps) > length(subject)
%         [STUDY, erpdata, erptimes, ~, pcond] = std_erpplot(...
%           STUDY,ALLEEG,'clusters',jj);
%       else
%         continue
%       end
%     end
% 
%     % BUG: EEGlab does not plot wave forms for fixed alpha in statparams ->
%     % workaround with alpha = nan in statparams:
%     tsigdiff = erptimes(pcond{1} < alpha); % sig. diff. time instances 
%     delete(subplot(122));
%     pos = get(subplot(121),'Position');
%     set(subplot(121),'Position',pos.*[1,1,2,1])
%     YLim = get(gca,'YLim');
%     hold on
%     plot(tsigdiff,(YLim(1)+.1*diff(YLim))*ones(length(tsigdiff),1),'r.')
%     if any(tsigdiff)
%       figERP(kk) = gcf;
%       ChLAB{kk} = lbl{jj};
%       kk = kk+1;
%     end
%   end
% end

% STUDY = pop_chanplot(STUDY,ALLEEG);

disp(NtrialsTab)

%% Print
if flags.do_print 
  FontSize = 8;
  Resolution = '-r600';
  set(findall([figTopo,figERP,figCorr],'-property','FontSize'),'FontSize',FontSize)

  if not(exist(mfilename,'dir'))
    mkdir(mfilename)
  end

  set(figTopo,'PaperUnits','centimeters','PaperPosition',[100,100,20,5])
  for ii = 1:length(figTopo)
    fn = fullfile(mfilename,['topo_',compLabel{ii},'_',conditionSet,'_',STUDY.design(design).name]);
    print(figTopo(ii),Resolution,'-dpng',fn)
    saveas(figTopo(ii),fn)
  end

  set(figERP,'PaperUnits','centimeters','PaperPosition',[100,100,12,8])
  for ii = 1:length(figERP)
    fn = fullfile(mfilename,['erp_',chanLabel{ii},'_',conditionSet,'_',STUDY.design(design).name]);
    print(figERP(ii),Resolution,'-dpng',fn)
    saveas(figERP(ii),fn)
  end

%   set(figCorr,'PaperUnits','centimeters','PaperPosition',[100,100,10,10])
%   for ii = 1:length(figCorr)
%     fn = fullfile(mfilename,['corr_dprime-',compLabel{ii},'AmpDiff_',conditionSet,'_',STUDY.design(design).name]);
%     print(figCorr(ii),Resolution,'-dpng',fn)
%     saveas(figCorr(ii),fn)
%   end
  if plotflag.corr
    fn = fullfile(mfilename,['compAmp_',conditionSet,'_',STUDY.design(design).name]);
    save(fn,'compAmp')
  end
  
end
