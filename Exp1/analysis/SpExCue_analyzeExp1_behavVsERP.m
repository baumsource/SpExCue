% export xls file for statistical analysis

change = load('SpExCue_analyzeExp1eeg_study/compAmp_Mchange_noD0');
onset = load('SpExCue_analyzeExp1eeg_study/compAmp_Onset_all.mat');
behav = load('SpExCue_analyzeExp1behav_avg_eeg.mat');
% behav = load('SpExCue_analyzeExp1behav_avg_behav.mat');
s = load('SpExCue_Exp1eeg_subjects_gender');


onset.idcond = [1,3,2]; % sucht that M = {0,0.5,1}
N1 = onset.compAmp(onset.idcond,:,1);
P2 = onset.compAmp(onset.idcond,:,2);
N1P2 = P2-N1;


change.idcond = [2,3,1,4,6,5]; % condition organization acc. to behav.meta.Dlabels
condition = change.condLbl(change.idcond);
Nsub = length(s.subject);
Ncond = length(condition);

cN1 = change.compAmp(change.idcond,:,1);
ncN1 = cN1./repmat(mean(cN1),[Ncond,1]);
zcN1 = zscore(cN1);
cP2 = change.compAmp(change.idcond,:,2);
ncP2 = cP2./repmat(mean(cP2),[Ncond,1]);
zcP2 = zscore(cP2);
cN1cP2 = cP2-cN1;

subjectRep = repmat(s.subject,[Ncond,1]);
conditionRep = repmat(condition,[Nsub,1])';
pFarther = squeeze(nanmean(behav.pFarther,2));
pCloser = 100-pFarther;
pCorrect = cat(1,pCloser(1:3,:),pFarther(4:6,:));
mpCorrect = mean(pCorrect);
zpCorrect = zscore(pCorrect);
% bias = mean(pCloser)./mean(pFarther);
bias = (mean(pCloser)-mean(pFarther))/100;
dprime = nanmean(behav.dprime);


%%
% nRangeN1 = range(N1)./abs(mean(N1));
% nRangeP2 = range(P2)./(mean(P2));
% phys = range(N1P2)./mean(N1P2);
% beh = dprime;
% [r,p] = corr(phys',beh','type','Spearman')
% figure; scatter(phys',beh','k.');hold on; xlabel('Physiology'); ylabel('Behavior')
% for ss = 1:Nsub
%   text(phys(ss),beh(ss),[' ',s.subject(ss).name,s.subject(ss).gender])
% end
%%
% SDncP2 = std(ncP2);
% [r,p] = corr(SDncP2',mpCorrect','type','Spearman') % ARO abstract
% figure; scatter(SDncP2',mpCorrect','k.');hold on; xlabel('cP2 variability'); ylabel('Percent correct')
% for ss = 1:Nsub
%   text(SDncP2(ss),mpCorrect(ss),[' ',s.subject(ss).name,s.subject(ss).gender])
% end
% 
% XLim = [0,4];
% YLim = [40,110];
% % plot(XLim,[1,1],'k:')
% % plot([0,0],YLim,'k:')
% axis([XLim,YLim])
% box on
% % nP2diff = mean(nP2(1:3,:))-mean(nP2(4:6,:));
% % [r,p] = corr(nP2diff',nanmean(behav.bias)','type','Spearman')
%%
zcP2diff = mean(zcP2(1:3,:))-mean(zcP2(4:6,:)); 
fig = figure; 
Color = 200/255*[1,0,1;0,1,1];
Symb = [char(9792);char(9794)]; % venus symbol; mars symbol
Rotation = [0,-45];
g = [s.subject.gender] == 'f';
XLim = [-1.4,1.7];
% YLim = [0.61,1.99];
YLim = [-0.16,0.39];
dx = 0.01*diff(XLim);
dy = 0.015*diff(YLim);
for ss = 1:Nsub
  idg = 2-g(ss); % 1...female, 2...male
  text(zcP2diff(ss),bias(ss),Symb(idg),...
    'HorizontalAlignment','center','Rotation',Rotation(idg),'Color',Color(idg,:))
  hold on
  
%   if g(ss) == 1 % female
%     text(zcP2diff(ss),bias(ss)-dy,'+','Color',Color(idg,:),'HorizontalAlignment','center')
%   else
%     text(zcP2diff(ss)+dx,bias(ss)+dy,'\rightarrow','Color',Color(idg,:),...
%       'Rotation',45,'HorizontalAlignment','center')
%   end
%   hold on
%   h=plot(zcP2diff(ss),bias(ss),Symb(idg));
%   if g(ss) == 1 % female
%     set(h,'MarkerFaceColor',Color(idg,:))
%   else
%     set(h,'MarkerFaceColor','w')
%   end
%   set(h,'Color',Color(idg,:))
  
%   text(zcP2diff(ss),bias(ss),[' ',s.subject(ss).name],'Color',Color(idg,:))
end
XLabel = 'cP2 looming bias (z-score)';
xlabel(XLabel); 
YLabel = 'Behavioral looming bias';
ylabel(YLabel)
dy = 5*dy;
dx = 5*dx;

TextColor = [zeros(1,3);Color];
idx = {g|not(g),g,not(g)};
for ii = 1:length(idx)
% [r,p] = corr(zcP2diff',bias','type','Spearman');
% text(XLim(1)+dx,YLim(2)-dy,...
%   ['\rho = ',num2str(r,'%1.2f'),', p = ',num2str(p,'%1.3f')],...
%   'Color','k')

  [r,p] = corr(zcP2diff(idx{ii})',bias(idx{ii})'); %,'type','Spearman'
  tmp = num2str(r,'%1.2f');
  if strcmp(tmp(1),'-') % eliminate initial zero but keep negative sign
    rstr = ['r = ',tmp([1,3:end])];
  else
%   rstr = ['\rho = ',tmp(2:end)]; 
    rstr = ['r = ',tmp(2:end)];
  end
  if p < .01
    pstr = 'p < .01';
  else
    tmp = num2str(p,'%1.2f');
    pstr = ['p = ',tmp(2:end)]; % eliminate initial zero
  end
  x = XLim(1)+dx;
  y = YLim(2)-dy-(ii-1)*0.1*range(YLim);
%   if ii > 1
%     text(x,y,Symb(ii-1),'Rotation',Rotation(ii-1))
%     text(x+dx/2,y,':')
%     x = x+dx;
%   end
  text(x,y,[rstr,', ',pstr],'Color',TextColor(ii,:))
end
% [r,p] = corr(zcP2diff(not(g))',bias(not(g))','type','Spearman');
% text(Symb(2),'Rotation',Rotation(2))
% text(XLim(1)+dx,YLim(2)-dy-0.2,...
%   [': \rho = ',num2str(r,'%1.2f'),', p = ',num2str(p,'%1.3f')],...
%   'Color',Color(2,:))

plot(XLim,[0,0],'k:')
plot([0,0],YLim,'k:')
axis([XLim,YLim])
box on
saveas(fig,'corr_cP2_bias')
saveas(fig,'corr_cP2_bias.svg')
%%
Labels = flipud(Symb);

figure
boxplot(zcP2diff,g,'colors',flipud(Color),'labels',Labels)
ylabel(XLabel)
set(gca,'YLim',XLim)
saveas(gcf,'boxplot_cP2')

figure
boxplot(bias,g,'colors',flipud(Color),'labels',Labels)
ylabel(YLabel)
set(gca,'YLim',YLim)
saveas(gcf,'boxplot_bias')

% P2ratio = mean(P2(1:3,:))./mean(P2(4:6,:));
% [r,p] = corr(mean(nP2(1:3,:)-nP2(4:6,:))',nanmean(behav.dprime)','type','Spearman')
%% Linear regression
tbl = table(zcP2diff(:),bias(:),g(:),'VariableNames',{'cP2','bias','gender'});
% tbl = table(zscore(zcP2diff(:)),zscore(bias(:)),g(:),'VariableNames',{'cP2','bias','gender'});
lme = fitlme(tbl,'bias~cP2*gender')

%% Linear regression including wihtin data
g_M = repmat(g,[3,1]);
subjectRep_M = repmat(s.subject,[Ncond/2,1]);
subject_M = {subjectRep_M.name};
bias_M = (pCloser(1:3,:)+pCloser(6:-1:4,:)-pFarther(1:3,:)-pFarther(6:-1:4,:))/100;
zcP2diff_M = zcP2(1:3,:)-zcP2(6:-1:4,:); 
tbl = table(zcP2diff_M(:),bias_M(:),g_M(:),subject_M(:),...
  'VariableNames',{'cP2','bias','gender','subject'});
lme = fitlme(tbl,'bias~cP2*gender+(1|subject)')

%% female distribution fit and male sample comparison
ii = 2; X = [zcP2diff(idx{ii})',bias(idx{ii})']; 
GM = fitgmdist(X,1);
pdf(GM,X)

ii = 3; X2 = [zcP2diff(idx{ii})',bias(idx{ii})']; 
pdf(GM,X2)

%%
% VariableNames = {'subject','condition','pCorrect','zpCorrect',...
%   'N1','P2','nN1','nP2','zN1','zP2','N1P2'};
% tbl = table(subjectRep(:),conditionRep(:),pCorrect(:),zpCorrect(:),...
%   N1(:),P2(:),nN1(:),nP2(:),zN1(:),zP2(:),N1P2(:),...
%   'VariableNames',VariableNames);
% %
% lme = fitlme(tbl,'zpCorrect~zP2+(1|subject)')
% beta = fixedEffects(lme)
% figure
% scatter(zP2(:),zpCorrect(:),'.')
% axis equal
% xlabel('cP2 (z-score)')
% ylabel('Response accuracy (z-score)')
% % lme = fitlme(tbl,'pCorrect~N1*P2+(1|subject)');
% % disp(lme)


%% export for Darrin
% behavior = mat2cell(pCorrect,ones(Ncond,1),ones(Nsub,1));
% P2cell = mat2cell(erp.P2,ones(Ncond,1),ones(Nsub,1));
% A = [VariableNames;[subject(:),condition(:),behavior(:),P2cell(:)]];
% save(mfilename,'A')
% xlswrite([mfilename,'.xlsx'],A);


% t = array2table(pCorrect');
% within = table(erp.P2,'VariableNames',{'P2'});
% rm = fitrm(t,['Var1-Var' num2str(Ncond) ' ~ 1'],'WithinDesign',within);
% ranova(rm)
%%
% subject2 = repmat(s.subjects,[Ncond/2,1]);
% dprime = behav.dprime_M([3,1,2],:);
% P2avg = (erp.P2(1:3,:)+erp.P2(6:-1:4,:))/2;
% tbl2 = table(subject2(:),dprime(:),P2avg(:),'VariableNames',{'subject','dprime','P2avg'});
% lme2 = fitlme(tbl2,'dprime~P2avg+(1|subject)');
% disp(lme2)