function SpExCue_analyzePilot2(ID)
% Analysis of SpExCue_Pilot2

% ID = {'VB'};
% ID = ID{1};
if not(exist('ID','var'))
  ID = input('ID: ','s');
end

flags.do_save = true;

conditions = {'all','continuousNoise','noiseBurst','speech',90,0,-90}; %,'AMnoiseBurst'
XTickLabel = {'all','cont.','discont.','speech','left','front','right'}; %,'AM (discont.)'
conditions = {'all','continuousNoise'}; %,'AMnoiseBurst'
XTickLabel = {'all','cont.'}; %,'AM (discont.)'


%% Load data
tmp = load(fullfile('Results',['SpExCue_pilot2Result_' ID]));
subj = tmp.subj;

N = length(subj.resp); % # trials
% M2 = subj.Mcomb(:,2);
D = diff(subj.Mcomb,1,2);
I = -diff(subj.SPL,1,2); % intensity difference
% U = abs(D);
E = zeros(length(subj.resp),1); %sign(subj.resp/max(subj.resp)-1+eps); % relative externalization, -1...closer, 1...farther
E(subj.resp == 99) = -1; % C key -> closer
E(subj.resp == 102) = 1; % F key -> further
% RT = subj.RT;
% latEc = abs(subj.pos(:,1)); % lateral eccentricity
% rPhase = subj.rphase;

for jj = 1:length(conditions)

  % select stimulus condtion
  if isscalar(conditions{jj}); %position
    idselect = subj.pos(:,1) == conditions{jj};
  elseif strcmp(conditions{jj},'all')
    idselect = true(length(subj.source),1);
  else
    idselect = ismember(subj.source,conditions{jj});
  end

  %% Psychometrics
  
  % Psychometric function
  Dset = sort(unique(D(idselect)));
  for dd = 1:length(Dset)
    idx = idselect & D == Dset(dd);
    pFarther(dd,jj) = sum(E(idx) == 1)/sum(idx);
  end
  
  % Sensitivity (hit defined as farther judgment if D > 0)
  pHit = sum(D(idselect)>0 & E(idselect)>0) / sum(D(idselect)>0);
  zHit = norminv(pHit-eps,0,1);
  pFalseAlarm = sum(D(idselect)<0 & E(idselect)>0) / sum(D(idselect)<0);
  zFalseAlarm = norminv(pFalseAlarm+eps,0,1);
  dprime(jj,1) = zHit - zFalseAlarm;
  
  % Bias
  bias(jj,1) = -0.5*(zHit+zFalseAlarm);

  %% Correlation Analysis
  y = zscore(E(:));
  X = [ones(length(E),1),zscore([D(:),I(:)])]; % ,latEc

  % exclude same-M trials
  idM = D ~= 0;
  idx = idselect & idM;
  y = y(idx);
  X = X(idx,:);

  [b(jj,:),bint(jj,:,:),r,rint,stats] = regress(y,X);
  R2(jj) = stats(1);
  p(jj) = stats(3);

  if p(jj) < .001
    plabel{jj} = '***';
  elseif p(jj) < .01
    plabel{jj} = '**';
  elseif p(jj) < .05
    plabel{jj} = '*';
  else
    plabel{jj} = '';
  end

end

%% Sensitivity (dprime) analysis for all combinations of M
Mcomb = unique(subj.Mcomb,'rows');
Mcomb = Mcomb(diff(Mcomb,1,2)>0,:);
dprime_Mcomb = nan(size(Mcomb,1),1);
for ii=1:size(Mcomb,1)
  idfarther = subj.Mcomb(:,1) == Mcomb(ii,1) & subj.Mcomb(:,2) == Mcomb(ii,2);
  idcloser = subj.Mcomb(:,2) == Mcomb(ii,1) & subj.Mcomb(:,1) == Mcomb(ii,2);
  % Sensitivity (hit defined as farther judgment if D > 0)
  pHit = .5*(sum(E(idfarther)>0) / sum(idfarther)) + ...
         .5*(sum(E(idcloser)<0) / sum(idcloser));
  zHit = norminv(pHit-eps,0,1);
  pFalseAlarm = .5*(sum(E(idfarther)<0) / sum(idfarther)) + ...
                .5*(sum(E(idcloser)<0) / sum(idcloser));
  zFalseAlarm = norminv(pFalseAlarm+eps,0,1);
  dprime_Mcomb(ii) = zHit - zFalseAlarm;
end

%% Plots & Tables

% Psychometric functions
fig(1) = figure;
% subplot(211)
hAll = plot(Dset,pFarther(:,1),'k');
hold on
hSingle = plot(Dset,pFarther(:,2:end));
set(hAll,'LineWidth',2)
set(gca,'XTick',Dset)
xlabel('M_{change} - M_{onset}')
ylabel('% ´farther´ judgments')
legend(XTickLabel,'Location','northwest')

% Tables
Ncond = length(conditions);
tab_dprime_bias = table(dprime,bias,'RowNames',XTickLabel(1:Ncond));
disp(tab_dprime_bias)

[dprime_Mcomb_sort,id_sort] = sort(dprime_Mcomb);
tab_dprime_Mcomb = table(Mcomb(id_sort,1),Mcomb(id_sort,2),dprime_Mcomb_sort,...
  'VariableNames',{'smallerM','largerM','dprime'});
disp(tab_dprime_Mcomb)

% Correlation analysis
fig(2) = figure;
bar([R2(:),b(:,2:end)])
hold on
plot([1.5,1.5],[0,1],'Color',[.5,.5,.5])
plot([4.5,4.5],[0,1],'Color',[.5,.5,.5])
text(1:Ncond,repmat(0.9,[Ncond,1]),plabel,'HorizontalAlignment','center','FontWeight','bold')
set(gca,'XTickLabel',XTickLabel,'XLim',[0.5,Ncond+.5])
% title(ID{1})
ylabel('R^2 or regression coefficient')
leg = legend('R^2','M change','SPL change');
set(leg,'Location','northoutside','Orientation','Horizontal')

%% Print
if flags.do_save
  FontSize = 10;
  Resolution = '-r600';
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  
  fn = fullfile(mfilename,[mfilename,'_',ID]);
  set(fig(1),'PaperUnits','centimeters','PaperPosition',[100,100,10,10])
  print(fig(1),Resolution,'-dpng',[fn,'_psyFct'])
  
  save(fn,'tab_dprime_bias','tab_dprime_Mcomb','pFarther','R2','b')
  
  set(fig(2),'PaperUnits','centimeters','PaperPosition',[100,100,18,10])
  print(fig(2),Resolution,'-dpng',[fn,'_regress'])
end