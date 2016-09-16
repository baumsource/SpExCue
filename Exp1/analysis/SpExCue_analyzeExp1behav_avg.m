% evaluate behavioral results of Exp1

saveflag = true;

postfix = 'behav'; % famili, behav, or eeg
tmp = load('SpExCue_Exp1eeg_subjects.mat');
subject = tmp.subjects;% {'S01','S02','S03','S04'};

% postfix = 'eeg';
% subject = {'RB','S01','S02','S03','S04'};

%% Evaluate listener-specific results
[~,~,meta] = SpExCue_analyzeExp1behav(subject{1},postfix,0,0);
pFarther = nan(length(meta.Dset),length(meta.position),length(subject));
dprime = nan(length(meta.position),length(subject));
dprime_M = nan(length(meta.Mcomb),length(subject));
bias = dprime;
for ii = 1:length(subject)
  [pFarther(:,:,ii),stats] = SpExCue_analyzeExp1behav(subject{ii},postfix,0,0);
  dprime(:,ii) = stats.Pos.dprime;
  bias(:,ii) = stats.Pos.bias;
  dprime_M(:,ii) = stats.M.dprime;
end

%% Plot
colors = 0.75*fliplr(eye(3));
symb = '<^>';
dx = [-.05,0,.05];
boxplotSettings = {'boxstyle','outline','plotstyle','traditional',...
  'colors',colors,'labels',meta.positionLabel};

% psychometric Function
fig(1) = figure;
legendLabel = meta.positionLabel;
for ii = 1:length(meta.position)
  h(ii) = errorbar(meta.Dset+dx(ii),nanmean(pFarther(:,ii,:),3),nanstd(pFarther(:,ii,:),0,3),['-',symb(ii)]);
  set(h(ii),'Color',colors(ii,:),'MarkerFaceColor',colors(ii,:))
  hold on
  legendLabel{ii} = [meta.positionLabel{ii},' (N=',num2str(sum(not(isnan(pFarther(1,ii,:))))),')'];
end
set(gca,'XLim',[-1.2,1.2],'XTick',meta.Dset,'XTickLabel',meta.Dlabels)
xlabel('Spectral contrast change (M_{onset}\rightarrowM_{change})')
ylabel('% ´farther´ judgments')
set(gca,'YLim',[0,100])
if length(meta.position) > 1
  legend(legendLabel,'Location','northwest')
end

% dprime and bias
fig(2) = figure;
subplot(121)
boxplot(dprime',boxplotSettings{:});
set(gca,'YLim',[-2.4,5.1])
ylabel('dprime')
subplot(122)
boxplot(bias',boxplotSettings{:});
set(gca,'YLim',[-0.9,1.1])
ylabel('bias towards "closer"')

% paired t-test for difference between decreasing and increasing spectral
% contrasts
pCorCloser = 100-squeeze(nanmean(pFarther(1:3,:,:),2)); % percent correct closer judgments
pCorFarther = squeeze(nanmean(pFarther(4:6,:,:),2));
[h,p,ci,stats]=ttest(mean(pCorCloser),mean(pCorFarther));
if h == 1
  disp(['Judgments sig. more correct for decreasing (vs increasing) spectral contrasts ',...
    '(t = ',num2str(stats.tstat,'%1.2f'),', p = ',num2str(p,'%1.2f'),')'])
end

%% Save
if saveflag
  
  save([mfilename,'_',postfix],'pFarther','dprime','dprime_M','bias','stats','meta')
  
  FontSize = 10;
  Resolution = '-r600';
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  if not(exist(mfilename,'dir'))
    mkdir(mfilename)
  end
  set(fig(1),'PaperUnits','centimeters','PaperPosition',[100,100,10,10])
  print(fig(1),Resolution,'-dpng',fullfile(mfilename,['psyFct_',postfix]))
  if length(meta.position) == 1
    set(fig(2),'PaperUnits','centimeters','PaperPosition',[100,100,8,10])
  else
    set(fig(2),'PaperUnits','centimeters','PaperPosition',[100,100,12,10])
  end
  print(fig(2),Resolution,'-dpng',fullfile(mfilename,['dprime_bias_',postfix]))
end