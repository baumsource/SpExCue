
%% 
flags.do_includeRespTrials = true;
flags.do_print = true;
rmpath(which('eeglab'))
addpath('/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/MATLAB_general')
bdfPath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Pilots/EEG 2 (cf Getzmann2010)/Results/';
bdfFilename = 'SpExCue_EEGpilot2Result_RT_RB.bdf';

YLim = [-4.9,6.9];

% Trigger definitions
trgOnset = fliplr([1,2,3,5]); % event triggers (as a row vector)
trgOnsetLabel = fliplr({'M=0','M=.25','M=.5','M=1'});
Color = {'k','b','r','g'};
trgChange = 110+trgOnset; % leftward motion
trgChange = cat(2,trgChange,120+trgOnset); % rightward motion
trgChange = cat(2,trgChange,130+trgOnset); % static
trgChangeR = trgChange+100; % response trials
trgChangeLabel = {'Leftward motion','Rightward motion','Static','Motion'};

% Load, preprocess and epoch data
trg = [trgOnset,trgChange,trgChangeR];
tEpo = [-0.200 0.6]; % epoch time frame in sec
cf = [0.5 25]; % bandpass filter cutoffs
ERP = processBDF_RB([bdfPath,bdfFilename],trg,tEpo,cf,'EOGch',64+(3:6));

% Analysis based on vertex electrode (#48)
Cz = ERP;
Cz.erp = squeeze(Cz.erp(:,48,:));
Lepoch = size(Cz.erp,1);

%% Stimulus onset
idt = Cz.t >= -0.1 & Cz.t <= 0.1;
onset = zeros(Lepoch,length(trgOnset));
for ii=1:length(trgOnset)
  idx = unique(ceil(find(Cz.triggers(idt,:) == trgOnset(ii))/sum(idt))); % trial indices for trigger
  onset(:,ii) = mean(Cz.erp(:,idx),2);
end

fig(1) = figure;
fn{1} = [bdfFilename(1:end-4),'_onset'];
for ii = 1:size(onset,2)
  plot(Cz.t*1000,onset(:,ii),Color{ii})
  hold on
end
plot(Cz.t*1000,zeros(length(Cz.t),1),'k-')
plot([0,0],YLim,'k:')
set(gca,'YLim',YLim)
legend(trgOnsetLabel,'Location','northwest')
title('Stimulus onset')
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')

%% Stimulus change
change = zeros(Lepoch,length(trgChange));
for ii=1:length(trgChange)
  idx = unique(ceil(find(Cz.triggers(idt,:) == trgChange(ii))/sum(idt))); % trial indices for trigger
  if isempty(idx) && flags.do_includeRespTrials % bug: all response trials
    idx = ceil(find(Cz.triggers == trgChangeR(ii))/Lepoch);
  end
  change(:,ii) = mean(Cz.erp(:,idx),2);
end

change = reshape(change,[Lepoch,length(trgOnset),3]);
change = cat(3,change,mean(change(:,:,1:2),3)); % average across leftward and rightward motion

for ii=1:size(change,3)
  fig(1+ii) = figure;
  fn{1+ii} = [bdfFilename(1:end-4),'_change',trgChangeLabel{ii}];
  for jj=1:size(change,2)
    plot(Cz.t*1000+700,change(:,jj,ii),Color{jj})
    hold on
  end
  plot(Cz.t*1000+700,zeros(length(Cz.t),1),'k-')
  plot([700,700],YLim,'k:')
  set(gca,'YLim',YLim)
  legend(trgOnsetLabel,'Location','northwest')
  title(trgChangeLabel{ii})
  xlabel('Time (ms)')
  ylabel('Amplitude (\muV)')
end

%% Print
if flags.do_print
  FontSize = 8;
  Resolution = '-r600';
  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  set(fig,'PaperUnits','centimeters','PaperPosition',[100,100,8,8])
  if not(exist(mfilename,'dir'))
    mkdir(mfilename)
  end
  for ii = 1:length(fig)
    print(fig(ii),Resolution,'-dpng',fullfile(mfilename,fn{ii}))
  end
end