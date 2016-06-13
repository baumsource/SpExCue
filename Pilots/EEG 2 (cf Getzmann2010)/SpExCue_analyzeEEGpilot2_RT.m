function SpExCue_analyzeEEGpilot2_RT(ID)
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
tmp = load(fullfile('Results',['SpExCue_EEGpilot2Result_RT_' ID]));
subj = tmp.subj;

idRespTrial = not(isnan(subj.resp(:)));
N = sum(idRespTrial); % # trials
disp(['Number of trials: ',num2str(N,'%i')])
M = subj.M(idRespTrial);
idHit = subj.hit(idRespTrial);
RT = subj.RT(idRespTrial);

%%
% Ms = unique(M);
% for m = 1:length(Ms)
%   meanRT(m) = 1000*mean(RT(idx));
%   stdRT(m) = 1000*std(RT(idx));
% end
% figure
% errorbar(Ms,meanRT,stdRT)
% set(gca,'XLim',[-.2,1.2],'XTick',Ms)

%% Plot
fig = figure;
boxplot(RT,M)
set(gca,'YLim',[0,.7],'XTickLabel',{'0','1/3','1'})
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