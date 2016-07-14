% evaluate behavioral results of EEGpilot3

saveflag = true;

postfix = 'behav'; % famili, behav, or eeg
subject = {'S01','S02','S03','S04'};

% postfix = 'eeg';
% subject = {'RB','S01','S02','S03','S04'};

%% Evaluate listener-specific results
[~,~,meta] = SpExCue_analyzeEEGpilot3_behavior([subject{1},postfix],0,0);
pFarther = nan(length(meta.Dset),length(meta.position),length(subject));
dprime = nan(length(meta.position),length(subject));
bias = dprime;
for ii = 1:length(subject)
  [pFarther(:,:,ii),stats] = SpExCue_analyzeEEGpilot3_behavior([subject{ii},postfix],0,0);
  dprime(:,ii) = stats.Pos.dprime;
  bias(:,ii) = stats.Pos.bias;
end

%% Plot
colors = 0.75*[zeros(1,3);fliplr(eye(3))];
symb = 'o<^>';
dx = [0,-.05,-.01,.05];
boxplotSettings = {'boxstyle','outline','plotstyle','traditional',...
  'colors',colors,'labels',meta.positionLabel};

% psychometric Function
fig(1) = figure;
for ii = 1:length(meta.position)
  h(ii) = errorbar(meta.Dset+dx(ii),mean(pFarther(:,ii,:),3),std(pFarther(:,ii,:),0,3),['-',symb(ii)]);
  set(h(ii),'Color',colors(ii,:),'MarkerFaceColor',colors(ii,:))
  hold on
end
XTickLabel = num2str(meta.Dset(:),'%0.2f');
set(gca,'XLim',[-1.2,1.2],'XTick',meta.Dset,'XTickLabel',XTickLabel)
xlabel('D = M_{change} - M_{onset}')
ylabel('% ´farther´ judgments')
set(gca,'YLim',[0,100])
if length(meta.position) > 1
  legend(meta.positionLabel,'Location','north')
end
% emphasize all condition
set(h(1),'LineWidth',2)
uistack(h(1), 'top')

% dprime and bias
fig(2) = figure;
subplot(121)
boxplot(dprime',boxplotSettings{:});
set(gca,'YLim',[-2.7,2.7])
ylabel('dprime')
subplot(122)
boxplot(bias',boxplotSettings{:});
set(gca,'YLim',[-0.9,0.9])
ylabel('bias towards "closer"')

%% Save
if saveflag
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