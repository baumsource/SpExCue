function [pFarther,stats,meta] = SpExCue_analyzeExp1behav(ID,fnext,plotflag,saveflag,exD0flag)
% SpExCue_analyzeExp1behav - Analysis of behavioral data from SpExCue_EEGpilot3
%  Usage: [pFarther,stats,meta] = SpExCue_analyzeExp1behav(ID,fnext,plotflag,saveflag,exD0flag)
%
%  Input parameters:
%  ID:        subject ID 
%  fnext:     experimental postfix, e.g., eeg or behav
%  plotflag:  set true to plot results (default: true)
%  saveflag:  set true to save (print) results (default: false)
%  exD0flag:  set true to exclude D=0 conditions (default: true)
%
%  Output parameters:
%  pFarther:  percentage of "farther" judgments
%  stats:     .Pos: dprime, bias, pCorrect, R2, and consistency for
%                   different positions
%             .M: dprime for different combinations of spectral contrasts 
%  meta:      .Dset: set of spectral contrast differences
%             .position: azimuth 
%             .positionLabel: corresponding position labels 
%             .Mcomb: set undirected combinations of spectral contrasts 


if not(exist('ID','var'))
  ID = input('ID: ','s');
end
if not(exist('fnext','var'))
  fnext = 'behav';
end
if not(exist('plotflag','var'))
  plotflag = true;
end
if not(exist('saveflag','var'))
  saveflag = false;
end
if not(exist('exD0flag','var'))
  exD0flag = true;
end

%% Load data
tmp = load(fullfile('..','data',['Exp1behav_' ID,'_',fnext]));
subj = tmp.subj;
pos = unique(subj.pos(:,1));
Npos = length(pos);
conditions = {'left','front','right'};
PositionLabel = conditions;
PositionColors = 'bgr';
% if Npos > 1 % only one direction
%   for ii = 1:Npos
%     conditions{ii+1} = pos(ii,1);
%     if pos(ii,1) < 0
%       PositionLabel{ii+1} = 'right';
%     elseif pos(ii,1) == 0
%       PositionLabel{ii+1} = 'front';
%     elseif pos(ii,1) > 0
%       PositionLabel{ii+1} = 'left';
%     else
%       PositionLabel{ii+1} = 'none';
%     end
%   end
% end
Ncond = length(conditions);

if any(subj.Mcomb == .5)
  subj.Mcomb(subj.Mcomb == .5) = 2/3;
  Dlabels = {'1\rightarrow0','.5\rightarrow0','1\rightarrow.5','none',...
    '.5\rightarrow1','0\rightarrow.5','0\rightarrow1'};
elseif any(subj.Mcomb == 1/3)
  Dlabels = {'1\rightarrow0','1\rightarrow1/3','1/3\rightarrow0','none',...
    '0\rightarrow1/3','1/3\rightarrow1','0\rightarrow1'};
end

D = diff(subj.Mcomb,1,2);
Dset = sort(unique(D));
if exD0flag || not(any(Dset == 0))
  Dset = Dset(Dset ~= 0);
  Dlabels = Dlabels([1:3,5:7]);
end
I = -diff(subj.SPL,1,2); % intensity difference
E = subj.E; 

%% Regression models and psychometrics for different directions
R2 = nan(length(conditions),1);
p = nan(length(conditions),1);
dprime = nan(length(conditions),1);
bias = nan(length(conditions),1);
pCorrect = nan(length(conditions),1);
consistency = nan(length(conditions),1);
pFarther = nan(length(Dset),length(conditions));
idplot = false(length(conditions),1);
for jj = 1:Ncond

  % select stimulus condtion
  switch conditions{jj}
    case 'left'
      idCond = subj.pos(:,1) > 0;
    case 'front'
      idCond = subj.pos(:,1) == 0;
    case 'right'
      idCond = subj.pos(:,1) < 0;
  end
  
  if any(idCond) == 0
    continue
  end

%   % exclude same-M trials
  idMdiff = D ~= 0;
  idCondMdiff = idCond & idMdiff;
  
  % Correlation Analysis
  y = zscore(E(:));
  if var(I(:)) == 0;
    X = [ones(length(E),1),zscore(D(:))]; % ,latEc
  else
    X = [ones(length(E),1),zscore([D(:),I(:)])];
  end
  y = y(idCondMdiff);
  X = X(idCondMdiff,:);

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
  
  for dd = 1:length(Dset)
    idx = idCond & D == Dset(dd);
    pFarther(dd,jj) = 100* sum(E(idx) == 1)/sum(E(idx) ~= 0);
  end
  
  % Consistency
  Uset = unique(abs(Dset));
  NiU = zeros(length(Uset),1);
  consistencyCounter = zeros(length(Uset),1);
  for uu = 1:length(Uset)
    iU = (abs(D) == Uset(uu)) & idCondMdiff;
    NiU(uu) = sum(iU);
    consistencyCounter(uu) = abs(sum(E(iU).*sign(D(iU))));
  end
  Nchange(jj) = sum(NiU);
  consistency(jj) = 100* sum(consistencyCounter) / Nchange(jj);
  
  % Sensitivity (hit defined as farther judgment if D > 0)
  idselect = idCond;
  pHit = sum(D(idselect)>0 & E(idselect)>0) / sum(D(idselect)>0 & E(idselect)~=0);
  pFalseAlarm = sum(D(idselect)<0 & E(idselect)>0) / sum(D(idselect)<0 & E(idselect)~=0);
  pCorrect(jj,1) = nansum(subj.hit(idselect)) / sum(not(isnan(subj.hit(idselect))));
  zHit = norminv(pHit-eps,0,1);
  zFalseAlarm = norminv(pFalseAlarm+eps,0,1);
  dprime(jj,1) = zHit - zFalseAlarm;
  
  % Bias
  bias(jj,1) = -0.5*(zHit+zFalseAlarm);
  
  legendLabel{jj} = [PositionLabel{jj},': d^{\prime}=',num2str(dprime(jj,1),'%1.2f')];
  idplot(jj) = true;
end

%% Consistency evaluated 
% if Ncond > 1
%   consistency(1) = Nchange(2:end) * consistency(2:end) / Nchange(1); % OR: consistency(1) = subj.consistency;
% end

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

%% Create tables
clear stats
Ncond = length(conditions);
stats.Pos = table(dprime,bias,pCorrect,R2,consistency,'RowNames',PositionLabel(1:Ncond));

[dprime_Mcomb_sort,id_sort] = sort(dprime_Mcomb);
stats.M = table(Mcomb(id_sort,1),Mcomb(id_sort,2),dprime_Mcomb_sort,...
  'VariableNames',{'smallerM','largerM','dprime'});

%% Show plots & tables

if plotflag
  % Psychometric functions
  fig(1) = figure;
  for ii = 1:Ncond
    h(ii) = plot(Dset,pFarther(:,ii),PositionColors(ii));
    hold on
  end
  legend(h(idplot),legendLabel(idplot),'Location','northwest')
  set(gca,'XTick',Dset,'XTickLabel',Dlabels,'YLim',[-3,103])
  xlabel('Spectral contrast change (M_{onset}\rightarrowM_{change})')
  ylabel('% ´farther´ judgments')
  title(ID)
  
  % display tables
  disp(stats.Pos)
  disp(stats.M)
end

%Correlation analysis
% fig(2) = figure;
% if var(I(:)) == 0;
%   bar(R2(:))
%   ylbl = 'R^2';
% else
%   bar([R2(:),b(:,2:end)])
%   ylbl = 'R^2, regression coefficient';
%   if var(I(:)) == 0;
%     leg = legend('R^2','M change');
%   else
%     leg = legend('R^2','M change','SPL change');
%   end
%   set(leg,'Location','northoutside','Orientation','Horizontal')
% end
% hold on
% plot([1.5,1.5],[0,1],'Color',[.5,.5,.5])
% plot([4.5,4.5],[0,1],'Color',[.5,.5,.5])
% text(1:Ncond,repmat(0.9,[Ncond,1]),plabel,'HorizontalAlignment','center','FontWeight','bold')
% set(gca,'XTickLabel',XTickLabel,'XLim',[0.5,Ncond+.5])
% % title(ID{1})
% ylabel(ylbl)

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

%% Print/Save
if saveflag
  FontSize = 10;
  Resolution = '-r600';
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  
  if not(exist(mfilename,'dir'))
    mkdir(mfilename)
  end
  fn = fullfile(mfilename,[mfilename,'_',ID]);
  set(fig(1),'PaperUnits','centimeters','PaperPosition',[100,100,10,10])
  print(fig(1),Resolution,'-dpng',[fn,'_psyFct','_',fnext])
  
%   save(fn,'tabs','pFarther','R2','b')
  
%   set(fig(2),'PaperUnits','centimeters','PaperPosition',[100,100,8,10])
%   print(fig(2),Resolution,'-dpng',[fn,'_regress'])
end

meta.Mcomb = Mcomb;
meta.Dset = Dset;
meta.position = conditions;
meta.positionLabel = PositionLabel(1:length(conditions));
meta.Dlabels = Dlabels;
end