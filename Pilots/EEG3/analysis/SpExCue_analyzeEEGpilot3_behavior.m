function SpExCue_analyzeEEGpilot3_behavior(ID)
% Analysis of behavioral data from SpExCue_EEGpilot3

% ID = {'VB'};
% ID = ID{1};
if not(exist('ID','var'))
  ID = input('ID: ','s');
end

flags.do_save = true;
conditions = {'all',90,0,-90};
Ncond = length(conditions);
XTickLabel = {'all','left','front','right'};

%% Load data
tmp = load(fullfile('..','data',['EEG3pilot_' ID]));
subj = tmp.subj;

D = diff(subj.Mcomb,1,2);
I = -diff(subj.SPL,1,2); % intensity difference
E = subj.E; 

%% Regression models for different directions
for jj = 1:length(conditions)

  % select stimulus condtion
  if isscalar(conditions{jj}); %position
    idCond = subj.pos(:,1) == conditions{jj};
  elseif strcmp(conditions{jj},'all')
    idCond = true(length(subj.E),1);
  end

%   % exclude same-M trials
  idMdiff = D ~= 0;
  idx = idCond;% & idMdiff;
  
  % Correlation Analysis
  y = zscore(E(:));
  if var(I(:)) == 0;
    X = [ones(length(E),1),zscore(D(:))]; % ,latEc
  else
    X = [ones(length(E),1),zscore([D(:),I(:)])];
  end
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

  % Psychometrics
  
  % Psychometric function
  idselect = idCond;
  Dset = sort(unique(D(idselect)));
  for dd = 1:length(Dset)
    idx = idselect & D == Dset(dd);
    pFarther(dd,jj) = sum(E(idx) == 1)/sum(idx);
  end
  
  % Sensitivity (hit defined as farther judgment if D > 0)
  pHit = sum(D(idselect)>0 & E(idselect)>0) / sum(D(idselect)>0);
  pFalseAlarm = sum(D(idselect)<0 & E(idselect)>0) / sum(D(idselect)<0);
  pCorrect(jj,1) = nansum(subj.hit(idselect)) / sum(not(isnan(subj.hit(idselect))));
  zHit = norminv(pHit-eps,0,1);
  zFalseAlarm = norminv(pFalseAlarm+eps,0,1);
  dprime(jj,1) = zHit - zFalseAlarm;
  
  % Bias
  bias(jj,1) = -0.5*(zHit+zFalseAlarm);
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
if Ncond > 1
  hold on
  hSingle = plot(Dset,pFarther(:,2:end));
end
set(hAll,'LineWidth',2)
set(gca,'XTick',Dset,'XTickLabel',{'-1','-2/3','-1/3','0','1/3','2/3','1'})
xlabel('D = M_{change} - M_{onset}')
ylabel('% ´farther´ judgments')
legend(XTickLabel,'Location','northwest')

% Tables
Ncond = length(conditions);
tab_dprime_bias = table(dprime,bias,pCorrect,'RowNames',XTickLabel(1:Ncond));
disp(tab_dprime_bias)

[dprime_Mcomb_sort,id_sort] = sort(dprime_Mcomb);
tab_dprime_Mcomb = table(Mcomb(id_sort,1),Mcomb(id_sort,2),dprime_Mcomb_sort,...
  'VariableNames',{'smallerM','largerM','dprime'});
disp(tab_dprime_Mcomb)

%Correlation analysis
fig(2) = figure;
if var(I(:)) == 0;
  bar(R2(:))
  ylbl = 'R^2';
else
  bar([R2(:),b(:,2:end)])
  ylbl = 'R^2, regression coefficient';
  if var(I(:)) == 0;
    leg = legend('R^2','M change');
  else
    leg = legend('R^2','M change','SPL change');
  end
  set(leg,'Location','northoutside','Orientation','Horizontal')
end
hold on
plot([1.5,1.5],[0,1],'Color',[.5,.5,.5])
plot([4.5,4.5],[0,1],'Color',[.5,.5,.5])
text(1:Ncond,repmat(0.9,[Ncond,1]),plabel,'HorizontalAlignment','center','FontWeight','bold')
set(gca,'XTickLabel',XTickLabel,'XLim',[0.5,Ncond+.5])
% title(ID{1})
ylabel(ylbl)

% % Correlation analysis
% fig(2) = figure;
% bar([R2(:),b(:,2:end)])
% % hold on
% % plot([1.5,1.5],[0,1],'Color',[.5,.5,.5])
% % plot([4.5,4.5],[0,1],'Color',[.5,.5,.5])
% text(1:Ncond,repmat(0.9,[Ncond,1]),plabel,'HorizontalAlignment','center','FontWeight','bold')
% set(gca,'XTickLabel',XTickLabel,'XLim',[0.5,Ncond+.5])
% % title(ID{1})
% ylabel('R^2 or regression coefficient')
% leg = legend('R^2','M change','SPL change');
% set(leg,'Location','northoutside','Orientation','Horizontal')

%% Print
if flags.do_save
  FontSize = 10;
  Resolution = '-r600';
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  
  if not(exist(mfilename,'dir'))
    mkdir(mfilename)
  end
  fn = fullfile(mfilename,[mfilename,'_',ID]);
  set(fig(1),'PaperUnits','centimeters','PaperPosition',[100,100,10,10])
  print(fig(1),Resolution,'-dpng',[fn,'_psyFct'])
  
  save(fn,'tab_dprime_bias','tab_dprime_Mcomb','pFarther','R2','b')
  
  set(fig(2),'PaperUnits','centimeters','PaperPosition',[100,100,8,10])
  print(fig(2),Resolution,'-dpng',[fn,'_regress'])
end