% cluster-based analysis for SpExCue Exp.1

addpath('/Users/rbaumgartner/Documents/MATLAB/fieldtrip-20160526')


%% Condtion settings

conditionSet = 'Mchange';%'Onset';

switch conditionSet
  case 'Onset'
    condition = {'M1','Mi','M0'};
  case 'OnsetToCloserVsFarther'
    condition = {'OnsetToCloser','OnsetToFarther'};
  case 'Mchange'
    condition = {'M1-0','M1-i','Mi-0','M0-1','Mi-1','M0-i'};
  case 'CloserVsFarther'
    condition = {'Closer','Farther'};
  otherwise
    error('RB: Check conditionSet!')
end


%% Load data
datapath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/data';
analysispath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/analysis';
load('SpExCue_Exp1eeg_subjects'); % subjects
fnData = fullfile(datapath,[mfilename,'_',conditionSet,'.mat']);
if exist(fnData,'file')
  tmp = load(fnData); % data cell array
  data = tmp.data;
else
  if not(exist('pop_loadset','file'))
    addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab/')
    eeglab
  end
  fnSxCx = 'Exp1eeg_Sxx_blinkICrej_Cxx.set';
  data = cell(length(subjects),length(condition));
  for isub = 1:length(subjects)
    fnCx = strrep(fnSxCx,'Sxx',subjects{isub});
    for icon = 1:length(condition)
      fn = strrep(fnCx,'Cxx',condition{icon});
      EEG = pop_loadset(fullfile(datapath,fn));
      data_allTrials = eeglab2fieldtrip( EEG, 'preprocessing', 'none' );
      tl.cfg.channel = 1:32; % only scalp EEG channels
      data{isub,icon} = ft_timelockanalysis(tl.cfg, data_allTrials);
    end
  end
  save(fnData,'data','condition')
end

Nsub = length(subjects);                       
Ncond = size(data,2);


%% General statistics settings
stat_cfg = [];
stat_cfg.method = 'montecarlo';
stat_cfg.correctm = 'cluster';
stat_cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that 
                                 % will be used for thresholding
stat_cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the 
                                 % permutation distribution. 
stat_cfg.minnbchan = 1;               % minimum number of neighborhood channels that is 
                                 % required for a selected sample to be included 
                                 % in the clustering algorithm (default=0).
stat_cfg.numrandomization = 500;      % number of draws from the permutation distribution
stat_cfg.latency       = [0 0.6];       % time interval over which the experimental 
                                 % conditions must be compared (in seconds)
stat_cfg.uvar  = 1;                    % row in design matrix specifying independent variable(s) -> conditions
stat_cfg.ivar  = 2;                 % row in design matrix specifying unit variable
% Define neighbouring channels
nb_cfg.method        = 'triangulation'; %, 'triangulation' or 'template'
stat_cfg.neighbours = ft_prepare_neighbours(nb_cfg, data{1});

switch conditionSet
  case 'Mchange'
    
    %% Statistics
    stat_cfg.statistic = 'depsamplesT'; % use the dependent samples T-statistic as a measure to 
                                     % evaluate the effect at the sample level
    stat_cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
    stat_cfg.clustertail = 0;
    stat_cfg.alpha = 0.025;               % alpha level of the permutation test

    % for bias
    stat_cfg.design = repmat(1:Nsub,[1,Ncond]);             % design row for subjects
    stat_cfg.design(2,:) = ceil((1:Nsub*Ncond)/(Ncond/2*Nsub)); % design row for bias

    fnStat = fullfile(analysispath,[mfilename,'_',conditionSet,'_stat.mat']);
    if exist(fnStat,'file')
      load(fnStat)
    else
      stat = ft_timelockstatistics(stat_cfg, data{:});
      save(fnStat,'stat')
    end

    %% Prepare data for plotting
    cfg = [];
    cfg.channel   = 'all';
    cfg.latency   = 'all';
    cfg.parameter = 'avg';
    Cdecrease = data(:,1:3);
    Cincrease = data(:,4:6);
    GA_Cd         = ft_timelockgrandaverage(cfg,Cdecrease{:});  
    GA_Ci         = ft_timelockgrandaverage(cfg,Cincrease{:});

    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    GA_Bias = ft_math(cfg,GA_Cd,GA_Ci);

    % get relevant (significant) values
    pos_cluster_pvals = [stat.posclusters(:).prob];
    pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
    pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
    % For each condition and subject
    ClustPot.mean = nan(15,6);
    ClustPot.rms = nan(15,6);
    ClustPot.gfp = nan(15,6);
    itime = data{1}.time>=stat.time(1) & data{1}.time<=stat.time(end);
    for icond = 1:Ncond
      for isub = 1:Nsub
        tmp = data{isub,icond}.avg(:,itime);
        ClustPot.mean(isub,icond) = mean(tmp(pos));
        ClustPot.max(isub,icond) = max(tmp(pos));
        ClustPot.rms(isub,icond) = rms(tmp(pos)); 
        ClustPot.gfp(isub,icond) = mean(std(tmp));
      end
    end
    save(strrep(fnStat,'stat','ClustPot'),'ClustPot')

    % get relevant time range
    timePoints = stat.time(any(pos,1)); % time range in seconds
    n = find(any(pos,1)); % time range in EEG samples
    timestep = diff(timePoints(1:2));
    % sigTimeRange(1) = min(stat.time(any(pos,1)));
    % sigTimeRange(2) = max(stat.time(any(pos,1)));
    % 
    % timestep = 1/data_allTrials.fsample;      %(in seconds)
    % timePoints = sigTimeRange(1):range(sigTimeRange)/length(stat.:sigTimeRange(end);   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    % n = round(timePoints*data_allTrials.fsample);  % temporal endpoints in EEG samples

    % get relevant bias range
    ZLim = max(abs(GA_Bias.avg(:)));

    %% Plot
    figure; 
    cfg = [];   
    % cfg_lay.elec = data_allTrials.elec;
    cfg_lay.layout = 'biosemi32';
    cfg.layout = ft_prepare_layout(cfg_lay,[]);% data_allTrials);
    cfg.zlim = ZLim*[-1,1];
    cfg.highlight = 'on';
    cfg.comment = 'no';  % xlim 
    cfg.commentpos = 'lefttop'; 
    cfg.colorbar = 'no'; 
    colormap jet
    for k = 1:2:length(timePoints);
       subplot(ceil(length(timePoints)/5/2),5,ceil(k/2)); 
       cfg.xlim= timePoints(k)+[-.5,.5]*timestep;  
       pos_int = all(pos(:, n(k)), 2);
       cfg.highlightchannel = find(pos_int);       
       ft_topoplotER(cfg, GA_Bias);
       axis tight
       title([num2str(timePoints(k)*1e3),' ms'])
    end  
    c = colorbar;
    set(c,'Position',[0.92 0.3 0.01 0.4])
    set(get(c,'YLabel'),'String','µV','Rotation',0)
    % RB_print(gcf,[16,5],'Exp1_ft_analysis_BiasCluster')
    saveas(gcf,'Exp1_ft_analysis_BiasCluster.svg')

    %% for contrast pair
%     stat_cfg.statistic = 'depsamplesFmultivariate';
%     stat_cfg.tail = 1;
%     stat_cfg.clustertail = 1;
%     stat_cfg.alpha = 0.05; 
%     % cfg.design = repmat(1:Nsub,[1,Ncond]);             % design row for subjects
%     % cfg.design(2,:) = ceil((1:Nsub*Ncond)/Nsub);
%     % contrastcoefs = repmat([1,-1,0 ; 1,0,-1],[1,2]); % compare C-pairs 
%     % cfg.contrastcoefs = repmat([1,-1,0 ; 1,0,-1 ; 0,1,-1],[1,2]); % compare C-pairs 
%     % [stat] = ft_timelockstatistics(cfg, data{:});%C01{:}, Ci1{:}, C0i{:});
% 
%     % Group across directions
%     stat_cfg.design = repmat(1:Nsub,[1,Ncond/2]); 
%     stat_cfg.design(2,:) = ceil((1:Nsub*Ncond/2)/Nsub);
%     dataC = cell(Nsub,Ncond/2);
%     mathcfg.operation = 'add';
%     mathcfg.parameter = 'avg';
%     for isub = 1:Nsub
%       for ipair = 1:Ncond/2
%         dataC{isub,ipair} = ft_math(mathcfg,data{isub,ipair},data{isub,Ncond/2+ipair});  
%       end
%     end
% 
%     fnStat = fullfile(analysispath,[mfilename,'_',conditionSet,'_CPair_stat.mat']);
%     if exist(fnStat,'file')
%       load(fnStat)
%     else
%       stat = ft_timelockstatistics(stat_cfg, dataC{:});
%       save(fnStat,'stat')
%     end
% 
%     pos_cluster_pvals = [stat.posclusters(:).prob];
%     pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
%     pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
%     sigtime = stat.time(any(pos,1))*1e3;
%     disp([sigtime(1),sigtime(end)]) % 90-320 ms

  case 'Onset'
    
    stat_cfg.statistic = 'depsamplesFmultivariate'; % use the dependent samples T-statistic as a measure to 
    stat_cfg.tail = 1;
    stat_cfg.clustertail = 1;
    stat_cfg.alpha = 0.05; 
    stat_cfg.design = repmat(1:Nsub,[1,Ncond]);             % design row for subjects
    stat_cfg.design(2,:) = ceil((1:Nsub*Ncond)/(Nsub)); % design row for bias

    fnStat = fullfile(analysispath,[mfilename,'_',conditionSet,'_stat.mat']);
    if exist(fnStat,'file')
      load(fnStat)
    else
      stat = ft_timelockstatistics(stat_cfg, data{:});
      save(fnStat,'stat')
    end
    
end