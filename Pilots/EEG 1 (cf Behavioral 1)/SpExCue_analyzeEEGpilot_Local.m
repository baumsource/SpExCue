clear
ID = 'RS';
%% 
flags.do_plot = true;
flags.do_print = true;
rmpath(which('eeglab'))
SOFAstart
addpath('/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Matlab_general')
bdfPath = fullfile(filesep,'Users','rbaumgartner','Documents','ARI','ARIcloud',...
  'SpExCue','Experiments','Pilots','EEG 1 (cf Behavioral 1)','Results');
bdfFilename = ['SpExCue_EEGpilot_',ID,'.bdf'];

YLimOnset = [-4.9,4.9];
YLimChange = [-4.9,4.9];% [-6.5,8.5]; % Plotted amplitude range

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
trgOnsetLabel = {'M: 1','M: .5','M: 0'};
onsetLineStyle = {'-','-.',':'};
trgChange = [t.M1_1  + t.M2_0 , t.M1_1  + t.M2_p5,...
             t.M1_p5 + t.M2_0 , t.M1_p5 + t.M2_1...
             t.M1_0  + t.M2_p5, t.M1_0  + t.M2_1];
trgChangeLabel = {'M: 1  -> 0 ','M: 1  -> .5','M: .5 -> 0 ','M: .5 -> 1 ','M: 0  -> .5','M: 0  -> 1 '};
changeLineStyle = {'-','-','-.','-.',':',':'};

% Load, preprocess and epoch data
tEpo = [-0.200 0.9]; % epoch time frame in sec
cf = [0.5 25]; % bandpass filter cutoffs
ERP = processBDF_RB([bdfPath,filesep,bdfFilename],[trgOnset,trgChange],tEpo,cf,'EOGch',32+3);

% Analysis based on vertex electrode (#32)
Cz = ERP;
Cz.erp = squeeze(Cz.erp(:,32,:));
Lepoch = size(Cz.erp,1);

%% Stimulus onset
% onset = zeros(Lepoch,length(trgOnset));
idt = Cz.t >= -0.1 & Cz.t <= 0.1;
for ii=1:length(trgOnset)
  idx = unique(ceil(find(Cz.triggers(idt,:) == trgOnset(ii))/sum(idt))); % trial indices for trigger
  Nidx(ii) = length(idx);
  onset(:,ii) = mean(Cz.erp(:,idx),2);
end

if flags.do_plot
  fig(1) = figure;
  fn{1} = [bdfFilename(1:end-4),'_Onset_left'];
  for ii = 1:size(onset,2)
    h = plot(Cz.t*1e3,onset(:,ii));
    set(h,'LineStyle',onsetLineStyle{ii})
    hold on
  end
  plot([0,0],YLimOnset,'k:')
  set(gca,'YLim',YLimOnset,'XLim',tEpo*1e3)
  legend(trgOnsetLabel,'Location','northwest')
  title('Stimulus onset')
  xlabel('Time (ms)')
  ylabel('Amplitude (\muV)')
end

%% Stimulus change
% change = zeros(Lepoch,length(trgChange),3);
% idt = find(Cz.t >= 0,1);
for ii=1:length(trgChange)
  idx = unique(ceil(find(Cz.triggers(idt,:) == trgChange(ii))/sum(idt))); % trial indices for trigger
  Nidx(ii) = 2*fix(length(idx)/2);
  change(:,ii,1) = mean(Cz.erp(:,idx),2);
  change(:,ii,2) = mean(Cz.erp(:,idx(1:Nidx(ii)/2)),2); % first half of trials
  change(:,ii,3) = mean(Cz.erp(:,idx(Nidx(ii)/2+1:end)),2); % second half of trials
end
trialSelectionLabel = {'',' (first half of trials)',' (second half of trials)'};

if flags.do_plot
  for ii=1:3
    fig(ii+1) = figure;
    fn{ii+1} = [bdfFilename(1:end-4),'_Change_left',trialSelectionLabel{ii}];
    for jj=1:size(change,2)
      h = plot(Cz.t*1000,change(:,jj,ii));
      set(h,'LineStyle',changeLineStyle{jj})
      hold on
    end
    plot([0,0],YLimChange,'k:')
    set(gca,'YLim',YLimChange,'XLim',tEpo*1e3)
    if ii>1
      legend(trgChangeLabel,'Location','northwest')
    end
    title(['Stimulus change',trialSelectionLabel{ii}])
    xlabel('Time (ms)')
    ylabel('Amplitude (\muV)')
  end
end

%% Print
if flags.do_plot && flags.do_print
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
end