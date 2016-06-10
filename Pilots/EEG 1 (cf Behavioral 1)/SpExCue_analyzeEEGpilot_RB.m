
%% 
flags.do_print = true;
rmpath(which('eeglab'))
% addpath('/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Matlab_general')
% bdfFilename = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Pilots/EEG 1 (cf Behavioral 1)/Results/SpExCue_EEGpilot_RB.bdf';
bdfFilename = 'SpExCue_EEGpilot_RB_07032016.bdf';

YLim = [-6.5,8.5]; % Plotted amplitude range

N1latency = [.120,.220]; % seconds
P2latency = [.220,.320]; % seconds

% Trigger definitions
t = struct(... 
  'correctResponse',8,...
  'wrongResponse',9,...
  'left',000,...
  'top',100,...
  'right',200,...
  'M1_1',30,... % M1 (at stimulus onset) of 1
  'M1_p5',20,...
  'M1_0',10,...
  'M2_1',3,... % M2 (at stimulus change) of 1
  'M2_p5',2,...
  'M2_0',1,...
  'stimulusOffset',9);
trgOnset = [t.M1_1,t.M1_p5,t.M1_0]; % event triggers (as a row vector)
% trgOnset = cat(2,[trgOnset+trigVals.top,trgOnset+trigVals.right]); % only left direction tested
trgOnsetLabel = {'M: 0','M: .5','M: 1'};
trgChange = [t.M1_1  + t.M2_0 , t.M1_1  + t.M2_p5,...
             t.M1_p5 + t.M2_0 , t.M1_p5 + t.M2_1...
             t.M1_0  + t.M2_p5, t.M1_0  + t.M2_1];
trgChangeLabel = {'M: 1  -> 0 ','M: 1  -> .5','M: .5 -> 0 ','M: .5 -> 1 ','M: 0  -> .5','M: 0  -> 1 '};
changeLineStyle = {'-','-','-.','-.',':',':'};

% Load, preprocess and epoch data
tEpo = [-0.200 0.9]; % epoch time frame in sec
cf = [0.5 25]; % bandpass filter cutoffs
ERPonset = processBDF_RB(bdfFilename,trgOnset,tEpo,cf,'EOGch',32+3);
ERPchange = processBDF_RB(bdfFilename,trgChange,tEpo,cf,'EOGch',32+3);

% Analysis based on vertex electrode (#32)
CzOnset = ERPonset;
CzOnset.erp = squeeze(CzOnset.erp(:,32,:));
LepochOnset = size(CzOnset.erp,1);
CzChange = ERPchange;
CzChange.erp = squeeze(CzChange.erp(:,32,:));
LepochChange = size(CzChange.erp,1);

%% Stimulus onset
onset = zeros(LepochOnset,length(trgOnset));
for ii=1:length(trgOnset)
  idx = ceil(find(CzOnset.triggers == trgOnset(ii))/LepochOnset); % trial indices for trigger
  onset(:,ii) = mean(CzOnset.erp(:,idx),2);
end

fig(1) = figure;
fn{1} = [mfilename,'_Onset_left'];
plot(Cz.t*1e3,onset)
hold on
plot([0,0],YLim,'k:')
set(gca,'YLim',YLim,'XLim',tEpo*1e3)
legend(trgOnsetLabel,'Location','northwest')
title('Stimulus onset')
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')

%% Stimulus change
change = zeros(LepochChange,length(trgChange),3);
for ii=1:length(trgChange)
  idx = ceil(find(CzChange.triggers == trgChange(ii))/LepochChange); % trial indices for trigger
  Nidx = fix(length(idx));
  change(:,ii,1) = mean(CzChange.erp(:,idx),2);
  change(:,ii,2) = mean(CzChange.erp(:,idx(1:Nidx/2)),2); % first half of trials
  change(:,ii,3) = mean(CzChange.erp(:,idx(Nidx/2+1:end)),2); % second half of trials
end
trialSelectionLabel = {'',' (first half of trials)',' (second half of trials)'};

for ii=1:3
  fig(ii+1) = figure;
  fn{ii+1} = [mfilename,'_Change_left',trialSelectionLabel{ii}];
  for jj=1:size(change,2)
  	h = plot(Cz.t*1000,change(:,jj,ii));
    set(h,'LineStyle',changeLineStyle{jj})
    hold on
  end
  plot([0,0],YLim,'k:')
  set(gca,'YLim',YLim,'XLim',tEpo*1e3)
  legend(trgChangeLabel,'Location','northwest')
  title(['Stimulus change',trialSelectionLabel{ii}])
  xlabel('Time (ms)')
  ylabel('Amplitude (\muV)')
end

%% Print
FontSize = 8;
Resolution = '-r600';
set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
set(fig,'PaperUnits','centimeters','PaperPosition',[100,100,12,8])
if not(exist(mfilename,'dir'))
  mkdir(mfilename)
end
for ii = 1:length(fig)
  print(fig(ii),Resolution,'-dpng',fullfile(mfilename,fn{ii}))
end