%% Evaluate behavioral results from SpExCue Exp. 3

fnPath = '../data';

% subjects
tmp = load('SpExCue_Exp3eeg_subjects.mat');
subjects = tmp.subject;
% subjects = subjects(2:4,:); disp('only S21, S15, and S28')

fnext = '_30deg_test.mat';
% fn = dir(fullfile(fnPath,'*_30deg_test.mat'));

seqsel = {1,2:3,1:3};
seqLbl = {{'1st syllable';'(simultaneous)'},{'2nd & 3rd syllables';'(staggered)'},'Whole sequence'};
condLbl = {'ITD','ILD','HRTF'};
pc = nan(2,height(subjects),length(seqsel));
% for ifn = 1:length(fn)
%   ID{ifn} = strrep(fn(ifn).name(1:3),'_','');
for ifn = 1:height(subjects)
  ID{ifn} = subjects.name{ifn};
  r = load(fullfile(fnPath,[subjects.name{ifn},fnext]));
  c = nan(length(r.result),length(seqsel));
  for ir = 1:length(r.result)
    for iseq = 1:length(seqsel)
      c(ir,iseq) = isequal(r.result(ir).seqtar(seqsel{iseq}),r.result(ir).response(seqsel{iseq}));
    end
  end
  for icond = 1:length(condLbl)
    I = ismember({r.result.spatialization},condLbl{icond});
    pc(icond,ifn,:) = mean(c(I,:));
  end
end
pc(:,ifn+1,:) = mean(pc,2);
ID = [ID,{'Avg'}];

%% Statistics
DV = array2table(pc(:,1:ifn,3)');
IVs = table(condLbl','VariableNames',{'cue'});
rm = fitrm(DV,['Var1-Var',num2str(length(condLbl)),' ~ 1'],'WithinDesign',IVs);
[ranovaResult,~,C,~] = ranova(rm,'WithinModel','cue');

ranovaResult.Properties.RowNames = strrep(ranovaResult.Properties.RowNames,'(Intercept):','');

% Sphericity corrections
spherCorr = epsilon(rm,C);
% Add corrected DFs to ranova table
idrep = round(0.5:0.5:length(spherCorr.GreenhouseGeisser)); % repeat iteratively
ranovaResult.DFGG = ranovaResult.DF .* ...
  reshape(spherCorr.GreenhouseGeisser(idrep),size(ranovaResult.DF));
% Add effect sizes to ranova table
SSeffect = ranovaResult.SumSq(1:2:end);
SSerror = ranovaResult.SumSq(2:2:end);
eta_pSq = nan(2*length(SSerror),1);
eta_pSq(1:2:end) = SSeffect./(SSeffect+SSerror); % effect size per (eta_partial)^2
ranovaResult.eta_pSq = eta_pSq;

disp(ranovaResult(:,[4,6,9,10]))

mc = multcompare(rm,'cue');
disp(mc)

%%
figure
hax = tight_subplot(1,3,0,.15,.1);
for iseq = 1:length(seqsel)
  axes(hax(iseq))
  plot(100*pc(:,1:end-1,iseq),':')
  hold on
  plot(100*pc(:,end,iseq),'k.-')
  set(gca,'XTick',1:3,'XTickLabel',condLbl)
  axis([0.5,3.5,41,105])
  title(seqLbl{iseq})
  if iseq == 1
    ylabel('Percent correct')
  else
    set(gca,'YTickLabel',{})
  end
  if iseq == 2
%     leg = legend(ID,'Location','south','Orientation','vertical');
    xlabel('Spatialization')
  end
end
leg = legend(ID,'Location','south','Orientation','vertical');
RB_print(gcf,[12,6],mfilename)