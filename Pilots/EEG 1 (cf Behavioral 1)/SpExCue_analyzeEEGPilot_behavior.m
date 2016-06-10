% analysis of SpExCuePilotExp

ID = {'RBspecRove10','RBHRTF3ms','RBHRTFeq'};
flags.do_print  = false;

for jj = 1:length(ID)

tmp = load(fullfile('Results',['SpExCue_EEGpilotResult_' ID{jj}]));
subj = tmp.subj;

N = length(subj.E); % # trials
M2 = subj.Mcomb(:,2);
D = diff(subj.Mcomb,1,2);
I = -diff(subj.SPL,1,2); % intensity difference
U = abs(D);
E = subj.E; % relative externalization, -1...closer, 1...farther
RT = subj.RT;
% latEc = abs(subj.pos(:,1)); % lateral eccentricity
rPhase = [subj.rphase];
drPhase = diff(rPhase,1,2);

y = zscore(E(:));
X = [ones(length(E),1),zscore([D(:),I(:)])]; % ,rPhase(:,2),latEc

% exclude same-M trials
idx = D ~= 0;
y = y(idx);
X = X(idx,:);

[b(jj,:),bint(jj,:,:),r,rint,stats] = regress(y,X);
R2(jj) = stats(1);
p(jj) = stats(3);

% errorbar(jj+0.25*[-1,0,1],b(2:end),bint(2:end,1),bint(2:end,2))

% disp(ID{jj})
% disp('All azimuths combined.')
% percentcorrect = sum(E==sign(D))/N;
% disp(['  Percent correct externalization responses: ' num2str(percentcorrect)])
% [r.E_D,p.E_D]=corrcoef(D,E);
% disp(['  Correlation between Ext. and D: ' num2str(r.E_D) ' (p = ' num2str(p.E_D) ')'])
% [r.E_I,p.E_I]=corrcoef(I,E);
% disp(['  Correlation between Ext. and I: ' num2str(r.E_I) ' (p = ' num2str(p.E_I) ')'])
% [r.M2_RT(jj),pvalue.M2_RT(jj)]=corrcoef(M2,RT);
% disp(['  Correlation between M2 and RT: ' num2str(r.M2_RT) ' (p = ' num2str(pvalue.M2_RT) ')'])
% 
% mtx(jj,1) = percentcorrect;
% mtx(jj,2) = r.E_D;
% mtx(jj,3) = r.E_I;
% mtx(jj,4) = r.M2_RT;

% azi = unique(subj.pos(:,1));
% for ii = 1:length(azi)
%   idpos = subj.pos(:,1) == azi(ii);
%   disp(' ')
%   disp(['Only ' num2str(azi(ii)) ' azimuth.'])
%   percentcorrect = sum(E(idpos)==sign(D(idpos)))/sum(idpos);
%   disp(['  Percent correct externalization responses: ' num2str(percentcorrect)])
%   [r.E_D,p.E_D]=corrcoef(D(idpos),E(idpos));
%   disp(['  Correlation between Ext. and D: ' num2str(r.E_D) ' (p = ' num2str(p.E_D) ')'])
%   [r.E_I,p.E_I]=corrcoef(I(idpos),E(idpos));
%   disp(['  Correlation between Ext. and I: ' num2str(r.E_I) ' (p = ' num2str(p.E_I) ')'])
% end
end

p

%% Plot
fig = figure;
bar([R2(:),b(:,2:end)]);
% bar(1,R2)
% hold on
% bar(2,b(2),'FaceColor',[0 .5 .5])
% bar(3,b(3),'y')
% set(gca,'XTick',1:3,'XTickLabel',{'R^2','\beta_M','\beta_{SPL}'})
% xlabel('Listener')
set(gca,'XTick',1:length(ID),'XTickLabel',ID,'YLim',[-.4,1])
ylabel('R^2 or regression coefficient')
legend('R^2','M change','SPL change','Ripple phase','Location','northoutside')

%% Print
if flags.do_print
  FontSize = 10;
  Resolution = '-r600';
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  set(fig,'PaperUnits','centimeters','PaperPosition',[100,100,8,10])
  print(fig,Resolution,'-dpng',mfilename)
end