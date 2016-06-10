function SpExCue_analyzePilot4(ID)
% Analysis of SpExCue_Pilot2

if not(exist('ID','var'))
  ID = input('ID: ','s');
end

% ID = {'VB'};
% ID = ID{1};
flags.do_save = true;

% conditions = {'all','continuousNoise','noiseBurst','speech',90,0,-90}; %,'AMnoiseBurst'
% XTickLabel = {'all','cont.','discont.','speech','left','front','right'}; %,'AM (discont.)'


%% Load data
tmp = load(fullfile('Results',['SpExCue_EEGpilot4Result_' ID]));
subj = tmp.subj;

N = length(subj.resp(:)); % # trials
disp(['Number of trials: ',num2str(N,'%i')])
M = subj.M(:);
idHit = subj.hit(:);
RT = subj.RT(:);

%%
Ms = unique(M);
for m = 1:length(Ms)
  idx = M==Ms(m) & idHit;
  meanRT(m) = 1000*mean(RT(idx));
  stdRT(m) = 1000*std(RT(idx));
end

%% Plot
fig = figure;
errorbar(Ms,meanRT,stdRT)
set(gca,'XLim',[-.2,1.2],'XTick',Ms)
xlabel('Spectral Contrast, M')
ylabel('Reaction Time (ms)')

%% Print
if flags.do_save
  FontSize = 10;
  Resolution = '-r600';
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  
  fn = fullfile([mfilename,'_',ID]);
  set(fig(1),'PaperUnits','centimeters','PaperPosition',[100,100,10,8])
  print(fig(1),Resolution,'-dpng',fn)
end