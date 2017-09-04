fnPath = '../data';

% subjects
tmp = load('SpExCue_Exp3eeg_subjects.mat');
subjects = tmp.subject;
% subjects = subjects(2:4,:); disp('only S21, S15, and S28')

fnext = '_30deg_test.mat';
% fn = dir(fullfile(fnPath,'*_30deg_test.mat'));

seqsel = {1,2:3,1:3};
seqLbl = {{'1st syllable';'(simultaneous)'},{'2nd & 3rd syllables';'(staggered)'},'Whole sequence'};
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
  iITD = ismember({r.result.spatialization},'ITD');
  iILD = ismember({r.result.spatialization},'ILD');
  iHRTF = ismember({r.result.spatialization},'HRTF');
  pc(1,ifn,:) = mean(c(iITD,:));
  pc(2,ifn,:) = mean(c(iILD,:));
  pc(3,ifn,:) = mean(c(iHRTF,:));
end
pc(:,ifn+1,:) = mean(pc,2);
ID = [ID,{'Avg'}];

%%
figure
hax = tight_subplot(1,3,0,.15,.1);
for iseq = 1:length(seqsel)
  axes(hax(iseq))
  plot(100*pc(:,1:end-1,iseq),':')
  hold on
  plot(100*pc(:,end,iseq),'k.-')
  set(gca,'XTick',1:3,'XTickLabel',{'ITD','ILD','HRTF'})
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
leg = legend(ID,'Location','southoutside','Orientation','vertical');
RB_print(gcf,[12,6],mfilename)