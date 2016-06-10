function SpExCue_analyzePilot3(ID)
% Analysis of SpExCue_Pilot3

% ID = {'VB'};
% ID = ID{1};
if not(exist('ID','var'))
  ID = input('ID: ','s');
end

flags.do_save = true;

conditions = {'all',90,0,-90}; %,'AMnoiseBurst'
Ncond = length(conditions);
XTickLabel = {'all','left','front','right'}; %,'AM (discont.)'

M_KEMAR = 0.75; % dummy M for KEMAR

%% Load data
tmp = load(fullfile('Results',['SpExCue_pilot3Result_' ID]));
subj = tmp.subj;

D = diff(subj.Mcomb,1,2);
I = -diff(subj.SPL,1,2); % intensity difference
E = zeros(length(subj.resp),1); %sign(subj.resp/max(subj.resp)-1+eps); % relative externalization, -1...closer, 1...farther
E(subj.resp == 99) = -1; % C key -> closer
E(subj.resp == 102) = 1; % F key -> further

%% Regression models for different directions
for jj = 1:length(conditions)

  % select stimulus condtion
  if isscalar(conditions{jj}); %position
    idCond = subj.pos(:,1) == conditions{jj};
  elseif strcmp(conditions{jj},'all')
    idCond = true(length(subj.source),1);
  end

  % exclude same-M trials
  idMdiff = D ~= 0;
  
  % KEMAR trials
  idKEMARfirst = subj.Mcomb(:,1) == M_KEMAR;
  idKEMARsecond = subj.Mcomb(:,2) == M_KEMAR;
  idKEMAR = idKEMARfirst|idKEMARsecond;
  
  % Correlation Analysis
  y = zscore(E(:));
  if var(I(:)) == 0;
    X = [ones(length(E),1),zscore(D(:))]; % ,latEc
  else
    X = [ones(length(E),1),zscore([D(:),I(:)])];
  end

  idx = idCond & idMdiff & not(idKEMAR);
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

  % KEMAR Analysis
  idx1 = idKEMARfirst & idCond;
  idx2 = idKEMARsecond & idCond;
  trueMs = unique(subj.Mcomb(idx2,1));
  for m = 1:length(trueMs)
    idx1M = idx1 & subj.Mcomb(:,2) == trueMs(m);
    idx2M = idx2 & subj.Mcomb(:,1) == trueMs(m);
    E_KEMAR(m,jj) = mean(-E(idx1M) + E(idx2M))/2;
  end

end

%% Plot Correlation analysis
fig(1) = figure;
if var(I(:)) == 0;
  bar(R2(:))
  ylbl = 'R^2';
else
  bar([R2(:),b(:,2:end)])
  ylbl = 'R^2, regression coefficient';
  leg = legend('R^2','M change','SPL change');
  set(leg,'Location','northoutside','Orientation','Horizontal')
end
hold on
plot([1.5,1.5],[0,1],'Color',[.5,.5,.5])
plot([4.5,4.5],[0,1],'Color',[.5,.5,.5])
text(1:Ncond,repmat(0.9,[Ncond,1]),plabel,'HorizontalAlignment','center','FontWeight','bold')
set(gca,'XTickLabel',XTickLabel,'XLim',[0.5,Ncond+.5])
% title(ID{1})
ylabel(ylbl)

%% Plot KEAMR analysis

fig(2) = figure;
h = plot(trueMs,E_KEMAR,'o-');
hold on
plot([-1,2],[0,0],'k:')
axis([-.1,1.1,-1.1,1.1])
set(gca,'XTick',trueMs)
set(h,'MarkerFaceColor','w')
legend(XTickLabel)
xlabel('Spectral contrast, M')
ylabel('Relative externalization')
title('Degree of externalization for KEMAR')

%% Sensitivity (dprime) analysis for all combinations of M
% Mcomb = unique(subj.Mcomb,'rows');
% Mcomb = Mcomb(diff(Mcomb,1,2)>0,:);
% dprime_Mcomb = nan(size(Mcomb,1),1);
% for ii=1:size(Mcomb,1)
%   idfarther = subj.Mcomb(:,1) == Mcomb(ii,1) & subj.Mcomb(:,2) == Mcomb(ii,2);
%   idcloser = subj.Mcomb(:,2) == Mcomb(ii,1) & subj.Mcomb(:,1) == Mcomb(ii,2);
%   % Sensitivity (hit defined as farther judgment if D > 0)
%   pHit = .5*(sum(E(idfarther)>0) / sum(idfarther)) + ...
%          .5*(sum(E(idcloser)<0) / sum(idcloser));
%   zHit = norminv(pHit-eps,0,1);
%   pFalseAlarm = .5*(sum(E(idfarther)<0) / sum(idfarther)) + ...
%                 .5*(sum(E(idcloser)<0) / sum(idcloser));
%   zFalseAlarm = norminv(pFalseAlarm+eps,0,1);
%   dprime_Mcomb(ii) = zHit - zFalseAlarm;
% end

%% Plots & Tables

% Psychometric functions
% fig(1) = figure;
% % subplot(211)
% hAll = plot(Dset,pFarther(:,1),'k');
% hold on
% hSingle = plot(Dset,pFarther(:,2:end));
% set(hAll,'LineWidth',2)
% set(gca,'XTick',Dset)
% xlabel('M_{change} - M_{onset}')
% ylabel('% ´farther´ judgments')
% legend(XTickLabel,'Location','northwest')

% Tables
% Ncond = length(conditions);
% tab_dprime_bias = table(dprime,bias,'RowNames',XTickLabel(1:Ncond));
% disp(tab_dprime_bias)

% [dprime_Mcomb_sort,id_sort] = sort(dprime_Mcomb);
% tab_dprime_Mcomb = table(Mcomb(id_sort,1),Mcomb(id_sort,2),dprime_Mcomb_sort,...
%   'VariableNames',{'smallerM','largerM','dprime'});
% disp(tab_dprime_Mcomb)
% 
% % Correlation analysis
% fig(2) = figure;
% bar([R2(:),b(:,2:end)])
% hold on
% plot([1.5,1.5],[0,1],'Color',[.5,.5,.5])
% plot([4.5,4.5],[0,1],'Color',[.5,.5,.5])
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
  print(fig(1),Resolution,'-dpng',[fn,'_regress'])
  
%   save(fn,'tab_dprime_bias','tab_dprime_Mcomb','pFarther','R2','b')
  
  set(fig(2),'PaperUnits','centimeters','PaperPosition',[100,100,10,10])
  print(fig(2),Resolution,'-dpng',[fn,'_KEMAR'])
end