%% Plot GFPs

addpath('/Users/rbaumgartner/Documents/MATLAB/fieldtrip-20160526')

conditionSet = 'MChange';
timeRange = [-199,400];

datapath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/data';
fnData = fullfile(datapath,['Exp1_ft_analysis_',conditionSet,'.mat']);
tmp = load(fnData);
data = tmp.data;
condition = tmp.condition;
Ncond = length(condition);

cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';

%% GFP for all conditions
GA = cell(Ncond,1);
GFP = nan(length(data{1}.time),Ncond);
GFP_individ = nan(size(data,1),length(data{1}.time),Ncond);
for ii = 1:Ncond
	 GA{ii} = ft_timelockgrandaverage(cfg,data{:,ii}); 
   GFP(:,ii) = std(GA{ii}.avg);
   for isub = 1:size(data,1)
      GFP_individ(isub,:,ii) = std(data{isub,ii}.avg);
   end
end

%% Get sig. time range
analysispath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp1/analysis';
fnStat = fullfile(analysispath,['Exp1_ft_analysis_',conditionSet,'_stat.mat']);
load(fnStat)

pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
sigtime = stat.time(any(pos,1))*1e3; 

%% all conditions
figure
time = GA{1}.time*1e3;
YLim = [0,1];
if not(isempty(sigtime))
  x = [min(sigtime),max(sigtime)];
  p = patch([x,fliplr(x)],YLim([1,1,2,2]),0.8*ones(1,3),'FaceAlpha',.1,'EdgeAlpha',0.1);
  hold on
end
h = plot(time,GFP);
xlabel('Time (ms)')
ylabel('Global field power (µV)')
axis([timeRange,YLim])
box off
switch conditionSet
  case 'Mchange'
    Color = [170,0,150;180,180,0;0,150,150]/255;
    for ii = 1:length(h)/2
      set(h(ii+[0,3]),'Color',Color(ii,:))
    end
    set(h(1:3),'LineStyle','-')
    set(h(4:6),'LineStyle',':')
  case 'Onset'
    Color = 0.7*eye(3);
    for ii = 1:length(h)
      set(h(ii),'Color',Color(ii,:))
    end
end

RB_print(gcf,[8.5,6.5],['Exp1_ft_plotGA_',conditionSet])

%% only bias
if strcmp(conditionSet,'Mchange')
  figure
  time = GA{1}.time*1e3;
  YLim = [0,1];
  if not(isempty(sigtime))
    x = [min(sigtime),max(sigtime)];
    p = patch([x,fliplr(x)],YLim([1,1,2,2]),0.8*ones(1,3),'FaceAlpha',.1,'EdgeAlpha',0.1);
    hold on
  end
  GFBbias = [mean(GFP(:,1:3),2),mean(GFP(:,4:6),2)];
  h = plot(time,GFBbias,'k-');
  xlabel('Time (ms)')
  ylabel('Global field power (µV)')
  axis([timeRange,YLim])
  box off
  set(h(1),'LineStyle','-')
  set(h(2),'LineStyle',':')

  RB_print(gcf,[8.5,5.5],['Exp1_ft_plotGA_',conditionSet,'_onlyBias'])
end