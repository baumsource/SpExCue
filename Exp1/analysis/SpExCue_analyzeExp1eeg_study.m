clear

%% Settings

flags.do_blinkICrej = true;
conditionSet = 'Onset'; % 'CloserVsFarther' or 'Mchange' or 'Onset','OnsetToCloserVsFarther'
% all: 1
% grand average: 2
% additional designs for 'Mchange': 
% 3: M1-0vsM0-1, 4: M1-ivsMi-1, 5: Mi-0vsM0-i, 6: M1-ivsM0-i, 7: D<0vsD>0
% 8: M1-0vsMi-0, 9: M0-1vsMi-1, 10: M0-ivsM1-i
% 11: all with D=0 collapsed, 12: Mi-1vsMi-0 (no prior expectation)
% 13: all except D=0
design = 1;


  
iERPmax = 32; % Cz (evaluated by code above), use 27 for F4, 13 for Pz

% Statistical analysis
stats.mcorrect = 'none'; % 'none' or 'MonteCarloCluster'

% Plots
plotflag.topo = true;
plotflag.erp = true;
plotflag.compAmp = true;
flags.do_print = true;

% ERP evaluation option
flags.do_ihd = false; % evaluate inter-hemispheric differences in ERPs (right - left

% Precomp settings
precomp.redo = 'on';
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
  {'dipoles' 'norm' 0 'weight' 1} ; ... 
%   {'erp' 'npca' 5 'norm' 1 'weight' 1 'timewindow' [0 400]};...
%   {'scalp' 'npca' 5 'norm' 1 'weight' 0};...
%   {'spec' 'npca' 5 'norm' 1 'weight' 1 'freqrange' [2 20]};...
% 	{'ersp' 'npca' 10 'norm' 1 'weight' 1};...
  };
clus_num = 8; % number of clusters

% Data
filepath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/data';

% component amplitudes
compLabel = {'N1','P2'};
switch conditionSet% extracted from Cz grand average ERP
  case {'Onset' , 'OnsetToCloserVsFarther'}
    compWin = {[80,140],[140,280]}; 
    ylimTopo = 4*[-1,+1];
  case {'Mchange' , 'CloserVsFarther'}
    compWin = {[90,150],[150,290]}; 
    ylimTopo = 2*[-1,+1];
end

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
fnStudy = strrep(fnStudy,'.set','.study');
  
switch conditionSet
  case 'Onset'
    condition = {'M1','Mi','M0'};
  case 'OnsetToCloserVsFarther'
    condition = {'OnsetToCloser','OnsetToFarther'};
  case 'Mchange'
    condition = {'M1-0','M1-i','Mi-0','Mx-x','M0-i','Mi-1','M0-1'};
  case 'CloserVsFarther'
    condition = {'Closer','Farther'};
  otherwise
    error('RB: Check conditionSet!')
end

if exist(fullfile(filepath,fnStudy),'file') && strcmp(precomp.redo,'off')
  
  [STUDY,ALLEEG] = pop_loadstudy('filename', fnStudy,'filepath',filepath);
%   STUDY.currentdesign = design;
  
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
        'subject' subject{ss} 'condition' condition{cc}};
    end
  end

  % Save
  [STUDY,ALLEEG] = std_editset([],[],'commands',commands,'name','SpExCue',...
    'task', 'closer vs. farther for second of two consecutive noise bursts'); 
%   [STUDY,ALLEEG] = std_editset(STUDY, ALLEEG,... % apply dipfit residual var threshold of 15%
%     'commands',{{'inbrain' 'on' 'dipselect' 0.15}},'updatedat','on','rmclust','off' );
  STUDY.design(1).name = 'all';

  % Create additional designs
  STUDY = std_makedesign(STUDY, ALLEEG, 2, 'name', 'grandAvg',...
    'variable1', 'condition','subjselect',subject,'values1',{condition});
  if strcmp(conditionSet,'Mchange')
    MchangeDesigns = {'M1-0vsM0-1',{'M1-0','M0-1'};...
      'M1-ivsMi-1',{'M1-i','Mi-1'};...
      'Mi-0vsM0-i',{'Mi-0','M0-i'};...
      'M1-ivsM0-i',{'M1-i','M0-i'};...
      'D<0vsD>0',{{'M1-i','M1-0','Mi-0'},{'Mi-1','M0-1','M0-i'}};...
      'M1-0vsMi-0',{'M1-0','Mi-0'};...
      'M0-1vsMi-1',{'M0-1','Mi-1'};...
      'M0-ivsM1-i',{'M0-i','M1-i'};...
      'all_D0pooled',{'M1-i','M1-0','Mi-0','Mi-1','M0-1','M0-i','Mx-x'};...
      'Mi-1vsMi-0',{'Mi-1','Mi-0'};...
      'noD0',{'M1-i','M1-0','Mi-0','Mi-1','M0-1','M0-i'}};
    for dd = 1:length(MchangeDesigns)
      STUDY = std_makedesign(STUDY, ALLEEG, dd+2, 'name',MchangeDesigns{dd,1},...
        'variable1', 'condition', 'values1',MchangeDesigns{dd,2},'subjselect',subject);
    end
      
  end

  % Save
  [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG,'filename', fnStudy,'filepath',filepath);
  
end
[STUDY] = std_selectdesign(STUDY, ALLEEG, design);
% STUDY.currentdesign = design;

% update workspace variables and redraw EEGLAB window
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
% eeglab redraw

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
[STUDY,ALLEEG] = std_precomp(STUDY, ALLEEG, precomp.data,precompsettings{:});

if flags.do_cluster
  % Clusstering based on scalp map correlation to a specified template
  % (here: EEG set 1, IC 6)
%   [correlationMap,STUDY,ALLEEG] = corrmap(STUDY,ALLEEG,1,6,'th','.75');
  
  % Preclustering
  [STUDY,ALLEEG] = std_preclust(STUDY, ALLEEG, [], preclustsetting{:});

  % Clustering
  [STUDY,ALLEEG] = pop_clust(STUDY, ALLEEG,'algorithm','kmeans',...
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

%% ERP analysis
if flags.do_channels
  
  % grand average ERP: evaluate site with maximum deflection -> Cz 
%     for jj = 1:32; 
%       [STUDY, erpdata, erptimes, ~, pcond] = std_erpplot(...
%         STUDY,ALLEEG,'channels',ChanLbls(jj),'noplot','on'); 
%       ERPmax(jj) = max(max(abs(mean(erpdata{1},2)),abs(mean(erpdata{2},2)))); 
%     end
%     [~,iERPmax] = max(ERPmax);

  chanLabel = ChanLbls(iERPmax);
  
  figERP = [];
  if plotflag.erp
    STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
      'timerange',[-200,600],'averagechan','on','topotime',[]);
    
    [STUDY,erpdata,erptimes,~,pcond] = std_erpplot(STUDY,ALLEEG,'channels',chanLabel);
    figERP(1) = gcf;

    if length(pcond) >= 1
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
    
    std_plotcurve(erptimes, erpdata, 'plotconditions', 'together',...
      'plotstderr', 'on', 'figure', 'on');
    figERP(2) = gcf;
    
  end
  
  figTopo = [];
  if plotflag.topo
    STUDY = pop_erpparams(STUDY, 'plotconditions', 'apart','averagechan','off');
    for ii = 1:length(compWin)
      [STUDY, C(ii).erpdata, ~, ~, C(ii).pcond] = std_erpplot(STUDY,ALLEEG,...
        'channels',ChanLbls(1:32), 'topotime',compWin{ii},'ylim',ylimTopo);
      figTopo(ii) = gcf;
    end
%     topotimeLabel = [num2str(tsigdiff(1),'%u'),'-',num2str(tsigdiff(end),'%u'),'ms'];
%     topotime = [tsigdiff(1), tsigdiff(end)];
%     [STUDY, erpdata, ~, ~, pcond] = std_erpplot(STUDY,ALLEEG,...
%         'channels',ChanLbls(1:32), 'topotime',topotime,'ylim',[-2,2]);
%     figTopo = gcf;
  end
  
  figComp = [];
  if plotflag.compAmp
    
    if flags.do_ihd % Inter-hemispheric difference
      chL = [1:12,14,15];
      chR = 17:30;
      [STUDY,erpdataL,erptimes,~,pcond] = std_erpplot(STUDY,ALLEEG,...
        'channels',ChanLbls(chL),'noplot','on');
      [~,erpdataR] = std_erpplot(STUDY,ALLEEG,...
        'channels',ChanLbls(chR),'noplot','on');
      erpdata = cell(length(erpdataR),1);
      for jj=1:length(erpdata)
        erpdata{jj} = squeeze(mean(erpdataR{jj} - erpdataL{jj},2));
      end
    else
      [STUDY,erpdata,erptimes,~,pcond] = std_erpplot(STUDY,ALLEEG,...
        'channels',ChanLbls(iERPmax),'noplot','on');
    end
    
    behav = load('SpExCue_analyzeExp1behav_avg_eeg.mat');
      
    condLbl = STUDY.design(design).variable(1).value(:); % condition labels
    compAmp = nan(length(erpdata),length(STUDY.subject),length(compWin));
    for ii = 1:length(compWin) % N1 and P2
      
      iComp = erptimes >= compWin{ii}(1) & erptimes <= compWin{ii}(2);
      for jj=1:length(erpdata)
        compAmp(jj,:,ii) = mean(erpdata{jj}(iComp,:));
      end
      
      if length(erpdata) == 2
        [h,p,ci,stats] = ttest(compAmp(1,:,ii),compAmp(2,:,ii));
        disp([compLabel{ii},', ',STUDY.design(design).name,': ',...
          't = ',num2str(stats.tstat,'%1.2f'),', p = ',num2str(p,'%1.5f')])
        
      elseif length(erpdata) > 2
      
        t = array2table(compAmp(:,:,ii)');
        
        if strcmp(conditionSet,'Mchange') && ismember(STUDY.design(design).name,{'noD0'})
          direction(ismember(condLbl,{'M1-0','M1-i','Mi-0'})) = {'decrease'};
          direction(ismember(condLbl,{'M0-1','Mi-1','M0-i'})) = {'increase'};
          direction(ismember(condLbl,{'M0-0','Mi-i','M1-1'})) = {'constant'};
          
          combination(ismember(condLbl,{'M1-0','M0-1'})) = {'0,1'};
          combination(ismember(condLbl,{'M0-i','Mi-0'})) = {'0,i'};
          combination(ismember(condLbl,{'Mi-1','M1-i'})) = {'i,1'};
          
          figComp(ii) = figure; 
          idsort = [2,1,3,5,4,6];
          compAmpWithin = compAmp-repmat(mean(compAmp),[size(compAmp,1),1,1]);
          boxplot(compAmpWithin(idsort,:,ii)',{direction(idsort),combination(idsort)})
          ylabel(['Within-subject ',compLabel{ii},' amplitude difference (µV)'])
          
          within = table(direction(:),combination(:),'VariableNames',{'direction','combination'});
          rm = fitrm(t,['Var1-Var' num2str(size(compAmp,1)) ' ~ 1'],'WithinDesign',within);
          [tbl.ranova,~,C,~] = ranova(rm,'WithinModel','direction*combination');
        else
          within = table(condLbl,'VariableNames',{conditionSet});
          rm = fitrm(t,['Var1-Var' num2str(size(compAmp,1)) ' ~ 1'],'WithinDesign',within);
          [tbl.ranova,~,C,~] = ranova(rm);
        end
        
        tbl.ranova.Properties.RowNames = strrep(tbl.ranova.Properties.RowNames,'(Intercept):','');
    
        % Mauchly's test for sphericity
%         tbl.mauchly = mauchly(rm,C);
%         tbl.mauchly.Properties.RowNames = tbl.ranova.Properties.RowNames(1:2:end);
    
        % Sphericity corrections
        tbl.eps = epsilon(rm,C);
%         tbl.eps.Properties.RowNames = tbl.ranova.Properties.RowNames(1:2:end);
    
        % Add corrected DFs to ranova table
        idrep = round(0.5:0.5:length(tbl.eps.GreenhouseGeisser)); % repeat iteratively
        tbl.ranova.DFGG = tbl.ranova.DF .* ...
          reshape(tbl.eps.GreenhouseGeisser(idrep),size(tbl.ranova.DF));
    
        % Add effect sizes to ranova table
        SSeffect = tbl.ranova.SumSq(1:2:end);
        SSerror = tbl.ranova.SumSq(2:2:end);
        eta_pSq = nan(2*length(SSerror),1);
        eta_pSq(1:2:end) = SSeffect./(SSeffect+SSerror); % effect size per (eta_partial)^2
        tbl.ranova.eta_pSq = eta_pSq;
        
        % Display results
        disp(['Repeated-measures ANOVA for ' compLabel{ii},' (',STUDY.design(design).name,')'])
        disp(tbl.ranova)
%         disp('Mauchly test and sphericity corrections')
%         disp([tbl.mauchly,tbl.eps])
    
        % Post-hoc analysis
        if any(tbl.ranova.pValueGG(1:2:end) < .05)
          disp('Posthoc analysis')
          if exist('direction','var')
            if tbl.ranova.pValueGG(3) < .05
              tbl.posthoc.direction = multcompare(rm,'direction');
              disp(tbl.posthoc.direction)
            end
            if tbl.ranova.pValueGG(5) < .05
              tbl.posthoc.combination = multcompare(rm,'combination');
              disp(tbl.posthoc.combination)
            end
          else
            tbl.posthoc = multcompare(rm,{conditionSet});
            disp(tbl.posthoc)
          end
        end
        
      end
      
    end
  end
  
end

%% IC cluster analysis

% Explore clusters
% if flags.do_cluster
%   pop_clustedit(STUDY,ALLEEG);
%   pause
% end

if flags.do_cluster

  % plot topographies
  std_topoplot(STUDY,ALLEEG,'clusters',3:length(STUDY.cluster));
  figTopo = gcf;
  
  % plot dipole locations
  std_dipplot(STUDY,ALLEEG,'clusters',2:length(STUDY.cluster)); % including outlier cluster
  figDip = gcf;

  %% plot ERPs
  if plotflag.erp
    STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
      'timerange',[-200,600],'averagechan','on','filter',20);
    alpha = .05;
    for jj=3:length(STUDY.cluster)
      [STUDY, erpdata, erptimes, ~, pcond] = std_erpplot(...
            STUDY,ALLEEG,'clusters',jj);
      figERP(jj-2) = gcf;
      title(STUDY.cluster(jj).name)
    end
  end
  
  %% evaluate time-windowed activity
  if plotflag.compAmp
    
    P2 = cell(length(STUDY.cluster),1);
    for jj=1:length(STUDY.cluster)
      idt = STUDY.cluster(jj).erptimes >= compWin{2}(1) & ...
            STUDY.cluster(jj).erptimes <= compWin{2}(2);
      for ee=1:length(STUDY.cluster(jj).erpdata)
        P2{jj}(ee,:) = rms(STUDY.cluster(jj).erpdata{ee}(idt,:));
      end
    end
      
    % Statistics
    % within-model definition for rm ANOVA
    if strcmp(conditionSet,'Mchange') && ismember(STUDY.design(design).name,{'noD0'})
      condLbl = STUDY.design(design).variable(1).value(:); % condition labels
      direction(ismember(condLbl,{'M1-0','M1-i','Mi-0'})) = {'decrease'};
      direction(ismember(condLbl,{'M0-1','Mi-1','M0-i'})) = {'increase'};
      direction(ismember(condLbl,{'M0-0','Mi-i','M1-1'})) = {'constant'};
      combination(ismember(condLbl,{'M1-0','M0-1'})) = {'0,1'};
      combination(ismember(condLbl,{'M0-i','Mi-0'})) = {'0,i'};
      combination(ismember(condLbl,{'Mi-1','M1-i'})) = {'i,1'};
      combination(ismember(condLbl,{'M0-0'})) = {'0,0'};
      combination(ismember(condLbl,{'Mi-i'})) = {'i,i'};
      combination(ismember(condLbl,{'M1-1'})) = {'1,1'};
      within = table(direction(:),combination(:),'VariableNames',{'direction','combination'});
    end
      
    for jj=1:length(STUDY.cluster)
      if length(STUDY.cluster(jj).erpdata) == 2
        [h(jj),p(jj)] = ttest(P2{jj}(1,:),P2{jj}(2,:));
        
      elseif exist('within','var') && size(P2{jj},1) == length(within.direction)
        
        t = array2table(P2{jj}');
        
        rm = fitrm(t,['Var1-Var' num2str(size(P2{jj},1)) ' ~ 1'],'WithinDesign',within);
        [tbl.ranova,~,C,~] = ranova(rm,'WithinModel','direction*combination');
        
        tbl.ranova.Properties.RowNames = strrep(tbl.ranova.Properties.RowNames,'(Intercept):','');
    
        % Sphericity corrections
        tbl.eps = epsilon(rm,C);
    
        % Add corrected DFs to ranova table
        idrep = round(0.5:0.5:length(tbl.eps.GreenhouseGeisser)); % repeat iteratively
        tbl.ranova.DFGG = tbl.ranova.DF .* ...
          reshape(tbl.eps.GreenhouseGeisser(idrep),size(tbl.ranova.DF));
    
        % Add effect sizes to ranova table
        SSeffect = tbl.ranova.SumSq(1:2:end);
        SSerror = tbl.ranova.SumSq(2:2:end);
        eta_pSq = nan(2*length(SSerror),1);
        eta_pSq(1:2:end) = SSeffect./(SSeffect+SSerror); % effect size per (eta_partial)^2
        tbl.ranova.eta_pSq = eta_pSq;
        
        % Display results
%         disp(['Repeated-measures ANOVA for ',STUDY.cluster(jj).name,' (',STUDY.design(design).name,')'])
%         disp(tbl.ranova)
        p(jj) = tbl.ranova.pValue(3);
      end
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
  set(findall([figTopo,figERP,figComp],'-property','FontSize'),'FontSize',FontSize)

  if not(exist(mfilename,'dir'))
    mkdir(mfilename)
  end

  set(figTopo,'PaperUnits','centimeters')
  for ii = 1:length(figTopo)
    set(figTopo(ii),'PaperPosition',[100,100,7*(length(get(figTopo(ii),'Children'))-1),5])
    fn = fullfile(mfilename,['topo_',compLabel{ii},'_',conditionSet,'_',STUDY.design(design).name]);
    print(figTopo(ii),Resolution,'-dpng',fn)
    saveas(figTopo(ii),fn)
  end

  set(figERP,'PaperUnits','centimeters','PaperPosition',[100,100,12,8])
  for ii = 1:length(figERP)
    fn = fullfile(mfilename,['erp_',chanLabel{1},'_',conditionSet,'_',STUDY.design(design).name]);
    if ii == 2
      fn = [fn,'_stderr'];
    end
    print(figERP(ii),Resolution,'-dpng',fn)
    saveas(figERP(ii),fn)
  end

  set(figComp,'PaperUnits','centimeters','PaperPosition',[100,100,10,10])
  for ii = 1:length(figComp)
    fn = fullfile(mfilename,[compLabel{ii},conditionSet,'_',STUDY.design(design).name]);
%     print(figComp(ii),Resolution,'-dpng',fn)
    saveas(figComp(ii),fn)
    saveas(figComp(ii),[fn,'.eps'])
  end
  if plotflag.compAmp
    fn = fullfile(mfilename,['compAmp_',chanLabel{1},'_',conditionSet,'_',STUDY.design(design).name]);
    save(fn,'compAmp','condLbl')
  end
  
end
