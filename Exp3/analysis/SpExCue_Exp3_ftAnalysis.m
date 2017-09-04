% addpath('/Users/rbaumgartner/Documents/MATLAB/fieldtrip-20160526')
% ft_defaults

dirBDF = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp3/data/';%string with the filename
experimentString = 'Exp3eeg';  % String to label experiment number for the saved eeglab dataset files
nChansLocsFilename = 'biosemi37.lay';  % filename of the eeglab locations file

trg.offset = 65280; % trigger offset because first 8 bits were always set high

trg.attend = trg.offset+(17:4:61);
trg.unattend = trg.offset+(81:4:125);
% trg.HRTF = 
% trg.ILD = 
% trg.ITD = 

layoutFn = fullfile(dirBDF,nChansLocsFilename);

cfg = [];
% cfg.layout = 'biosemi32.lay';
cfg.layout = layoutFn;
layout = ft_prepare_layout(cfg);
ft_layoutplot(cfg)
%%
cfg = [];
% cfg.layout = layout;
cfg.dataset = fullfile(dirBDF,[experimentString,'_S15.bdf']);
% cfg.dataset = filename;
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = [trg.attend,trg.unattend];
cfg.trialdef.prestim    = 1;
cfg.trialdef.poststim   = 2;
cfg = ft_definetrial(cfg);
%
cfg.reref = 'yes';
cfg.refchannel = {'EXG1','EXG2'};
cfg.demean = 'yes';
cfg.baselinewindow = [-0.2 0];
data = ft_preprocessing(cfg);

%%
cfg = [];
cfg.channel = {'all','-EXG6','-EXG7','-EXG8'};
data = ft_selectdata(cfg, data);
data.label(1:37) = layout.label(1:37);
%%
cfg = [];
cfg.keeptrials = 'yes';
timelock = ft_timelockanalysis(cfg, data);

%% stats
cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to 
                                 % evaluate the effect at the sample level
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that 
                                 % will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the 
                                 % permutation distribution. 
cfg.minnbchan = 2;               % minimum number of neighborhood channels that is 
                                 % required for a selected sample to be included 
                                 % in the clustering algorithm (default=0).
cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.025;               % alpha level of the permutation test
cfg.numrandomization = 100;      % number of draws from the permutation distribution

design = zeros(1,size(timelock.trial,1));
design(ismember(timelock.trialinfo,trg.attend)) = 1;
design(ismember(timelock.trialinfo,trg.unattend)) = 2;

cfg.design = design;             % design matrix
cfg.ivar  = 1;                   % number or list with indices indicating the independent variable(s)

cfg_neighb        = [];
cfg_neighb.method = 'triangulation';
cfg_neighb.layout = layout;
neighbours        = ft_prepare_neighbours(cfg_neighb, data);

cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with 
                                 % which other sensors it can form clusters
                                 
ExclChans = {'-Tp9','-Tp10','-vEOG','-hEOGl','-hEOGr','-Status'};
cfg.channel       = {'all',ExclChans{:}};     % cell-array with selected channel labels
cfg.latency       = [-0.5 1];       % time interval over which the experimental 
                                 % conditions must be compared (in seconds)
                                 
[stat] = ft_timelockstatistics(cfg, timelock);

%% plot sig cluster
alpha = .05;
neg_signif_clust = find([stat.negclusters.prob] < alpha);
neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
neg_time = stat.time(any(neg,1));

cfg.highlight = 'on';
% Get the index of each significant channel
cfg.highlightchannel = find(any(neg,2));
cfg.comment = 'xlim';   
cfg.commentpos = 'title';   
cfg.layout = layout;
cfg.interactive = 'no';
cfg.xlim = [neg_time(1),neg_time(end)];
ft_topoplotER(cfg, timelock);  

%% ERP waveform
cfg = [];
cfg.showlabels = 'yes'; 
cfg.fontsize = 6; 
cfg.layout = layout;
% cfg.ylim = [-3e-13 3e-13];
% only one 
% ft_multiplotER(cfg, timelock); 
% or multiple ones overlayed
% ft_multiplotER(cfg, avgFC, avgIC, avgFIC);
% or only for the cannel of interest
cfg.channel = 'Cz';
clf;
ft_singleplotER(cfg,timelock);

%% T-F analysis
cfg              = [];
cfg.output       = 'pow';
% cfg.channel      = 'EEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:30;                         % analysis 2 to 30 Hz in steps of 2 Hz 
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.t_ftimwin    = 5./cfg.foi;                     % or 5 cycles per freq
cfg.toi          = -0.5:0.05:1.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
% cfg.tapsmofrq  = 0.4*cfg.foi; % freq-dependent multitapering (uncomment cfg.taper)
TFR = ft_freqanalysis(cfg, data);

%%
cfg = [];
cfg.baseline     = [-0.5 -0.1];
% cfg.baselinetype = 'absolute';  
% cfg.zlim         = [-3e-27 3e-27];	
cfg.baselinetype = 'relchange';  
% cfg.zlim         = [-2.0 2.0];
% cfg.maskstyle    = 'saturation';	
cfg.channel      = 'Cz';
figure 
ft_singleplotTFR(cfg, TFR);
%%
cfg = [];
cfg.baseline     = [-0.5 -0.1]; 
cfg.baselinetype = 'relchange'; 	        
% cfg.zlim         = [-3e-25 3e-25];
cfg.showlabels   = 'yes';	        
cfg.layout       = layout;
figure
ft_multiplotTFR(cfg, TFR)