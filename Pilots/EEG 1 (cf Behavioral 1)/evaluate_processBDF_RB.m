% evaluate processBDF_RB by comparison between SpExCue_analyzeEEGpilot_Local and
% SpExCue_analyzeEEGpilot_eeglab
% Problem: huge file and false triggers in SpExCue_EEGpilot_RB_07032016.bdf
% Solution: processing via Cluster and separate analysis

clear

load([mfilename,filesep,'SpExCue_analyzeEEGpilot_Cluster.mat'])
Cz.triggers = mod(Cz.triggers,100);

%%
epoch = [-.2,1.5];
trgOnset = 30:-10:10; % 10 -> M1=0, 30 -> M1=1
trgOnsetLabel = {'M: 1','M: .5','M: 0'};
trgChange = 3:-1:1; % 1 -> M2=0, 3 -> M2=1
trgChangeLabel = {'M: 1  -> 0 ','M: 1  -> .5','M: .5 -> 0 ','M: .5 -> 1 ','M: 0  -> .5','M: 0  -> 1 '};
XLim = [-200,800];
YLim = [-4.9,4.9];
% YLim = [-2.9,2.9]; % zoom

%% Onset response
idt = Cz.t >= -0.1 & Cz.t <= 0.1;
for ii=1:length(trgOnset)
  idx = unique(ceil(find(Cz.triggers(idt,:) == trgOnset(ii))/sum(idt))); % trial indices for trigger
  Nidx(ii) = length(idx);
  onset(:,ii) = mean(Cz.erp(:,idx),2);
end

onsetLineStyle = {'-','-.',':'};
fig(1) = figure;
for ii = 1:size(onset,2)
  h = plot(Cz.t*1e3,onset(:,ii));
  set(h,'LineStyle',onsetLineStyle{ii})
  hold on
end
plot([0,0],YLim,'k:')
set(gca,'YLim',YLim,'XLim',XLim,'XMinorTick','on')
legend(trgOnsetLabel,'Location','northwest')
title('Stimulus onset')
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')

%% Change response
idtc = Cz.t >= 0.5 & Cz.t <= 0.7; % time range for triggers
tBL = find(Cz.t>-0.1 & Cz.t<=0); % time range for 
Npre = fix(XLim(1)/1000*Cz.fs);
Npost = fix(XLim(2)/1000*Cz.fs)-1;
kk = 1;
for ii=1:length(trgOnset)
  id1 = ceil(find(Cz.triggers(idt,:) == trgOnset(ii))/sum(idt)); % trial indices for trigger
  id1 = unique(id1);
  tmp.triggers = Cz.triggers(:,id1);
  tmp.erp = Cz.erp(:,id1);
  for jj=mod(ii+[1,2]-1,3)+1
    id2 = ceil(find(tmp.triggers(idtc,:) == trgChange(jj))/sum(idtc)); % trial indices for trigger
    id2 = unique(id2);
    % align on triggers
    tmpERP = zeros(Npost-Npre+1,length(id2));
    for tt = 1:length(id2)
      idtcon = find(tmp.triggers(:,id2(tt)) == trgChange(jj),1);
      tmpERP(:,tt) = tmp.erp(idtcon+Npre:idtcon+Npost,id2(tt));
    end
    % baseline adjustment
    tmpERP = tmpERP-repmat(mean(tmpERP(tBL,:,:)),[size(tmpERP,1),1,1]);
    % average
    change(:,kk) = mean(tmpERP,2);
    kk = kk+1;
  end
end
change = change(:,[2,1,3,4,6,5]);

%%
changeLineStyle = {'-','-','-.','-.',':',':'};
fig(2) = figure;
t = XLim(1)/1e3:1/Cz.fs:XLim(2)/1e3;
t = t(1:length(change));
for jj=1:size(change,2)
  h = plot(t*1000,change(:,jj));
  set(h,'LineStyle',changeLineStyle{jj})
  hold on
end
plot([0,0],YLim,'k:')
set(gca,'YLim',YLim,'XLim',XLim,'XMinorTick','on')
legend(trgChangeLabel,'Location','northwest')
title('Stimulus change')
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')

%% Print
FontSize = 8;
Resolution = '-r600';
set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
set(fig,'PaperUnits','centimeters','PaperPosition',[100,100,12,8])
if not(exist(mfilename,'dir'))
  mkdir(mfilename)
end

print(fig(1),Resolution,'-dpng',fullfile(mfilename,'SpExCue_analyzeEEGpilot_Onset'))
print(fig(2),Resolution,'-dpng',fullfile(mfilename,'SpExCue_analyzeEEGpilot_Combinations'))
