% export xls file for statistical analysis

erp = load('SpExCue_analyzeExp1eeg_study/P2AmpDiff_Mchange_all');
behav = load('SpExCue_analyzeExp1behav_avg_eeg.mat');
s = load('SpExCue_Exp1eeg_subjects');

Nsub = length(s.subjects);
Ncond = length(erp.condition);

subject = repmat(s.subjects,[Ncond,1]);
condition = repmat(erp.condition,[Nsub,1])';
pCorrect = cat(1,100-squeeze(nanmean(behav.pFarther(1:3,:,:),2)),squeeze(nanmean(behav.pFarther(4:6,:,:),2)));
nP2 = erp.P2./repmat(mean(erp.P2),[Ncond,1]);
% idc = [4,7,6,8,3,2];set(gca,'XTick',1:6,'XTickLabel',STUDY.condition(idc))

[r,p] = corr(std(nP2)',mean(pCorrect)','type','Spearman') % ARO abstract
%%
VariableNames = {'subject','condition','pCorrect','nP2'};
tbl = table(subject(:),condition(:),pCorrect(:),nP2(:),'VariableNames',VariableNames);

lme = fitlme(tbl,'pCorrect~nP2+(1|subject)');
disp(lme)


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
subject2 = repmat(s.subjects,[Ncond/2,1]);
dprime = behav.dprime_M([3,1,2],:);
P2avg = (erp.P2(1:3,:)+erp.P2(6:-1:4,:))/2;
tbl2 = table(subject2(:),dprime(:),P2avg(:),'VariableNames',{'subject','dprime','P2avg'});
lme2 = fitlme(tbl2,'dprime~P2avg+(1|subject)');
disp(lme2)