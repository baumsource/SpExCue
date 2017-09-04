ERPanalysis = true;

currentDesign = 1;

precomp.redo = 'on';

designs = {...
  'attended',{'ITD','ILD','HRTF'},{'ITD','ILD','HRTF'};...
  'unattended',{'ITDu','ILDu','HRTFu'},{'ITDu','ILDu','HRTFu'};...
  'attendedVsUnattended',{{'ITD','ILD','HRTF'},{'ITDu','ILDu','HRTFu'}},{'attended','unattended'};...
%   'motion',{{'C0_left','C0_right'},{'C1_left','C1_right'},{'KEMAR_left','KEMAR_right'}},{'C = 0','C = 1','KEMAR'};...
%   'LR',{{'C0_left','C1_left','KEMAR_left'},{'C0_right','C1_right','KEMAR_right'}},{'left','right'}
};

%% Load/create STUDY
% Add paths
if not(exist('eeglab','file'))
  addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab')
  eeglab
end

experimentString = 'Exp3eeg';
datapath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp3/data/';%string with the filename

% subjects
tmp = load('SpExCue_Exp3eeg_subjects.mat');
subjects = tmp.subject;
% subjects = subjects(2:4,:); disp('only S21, S15, and S28')
Ns = height(subjects);

% configurations
tmp = load([experimentString,'_cnfg']);
cnfg = tmp.cnfg;
cond = cnfg.cond;

fnstudy = ['study_',experimentString,'.study'];
if ERPanalysis; fnstudy = strrep(fnstudy,'.','_ERP.'); end
if 0%exist(fullfile(datapath,fnstudy),'file')
  [STUDY, ALLEEG] = pop_loadstudy('filename', fullfile(datapath,fnstudy));
else
  fnX = [experimentString,'_Sxx_clean_Cxx.set'];
  commands = cell(1,height(subjects)*size(cond,1));
  fn = commands;
  ii = 0;
  for ss = 1:height(subjects)
    for cc = 1:size(cond,1)
      ii = ii+1;
      fn{ii} = strrep(fnX,'Sxx',subjects.name{ss});
      fn{ii} = strrep(fn{ii},'Cxx',cond{cc,2});
      commands{ii} = { 'index' ii 'load' fullfile(datapath,fn{ii}),...
        'subject' subjects.name{ss} 'condition' cond{cc,2}};
    end
  end

  % Create STUDY
  [STUDY,ALLEEG] = std_editset([],[],'commands',commands,'name','SpExCue Exp. 3',...
    'task', 'left- vs. right-hand /ba/-/da/-/ga/ stream attention'); 
  % Define designs
  STUDY.design(1).name = 'all';
  for dd = 1:size(designs,1)
    STUDY = std_makedesign(STUDY, ALLEEG, dd+1, 'name',designs{dd,1},...
      'variable1', 'condition', 'values1',designs{dd,2});%,'subjselect',{'S21'});
  end
  % Save
  [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG,'filename', fullfile(datapath,fnstudy));
end
%%
[STUDY] = std_selectdesign(STUDY, ALLEEG, currentDesign);

ChanLbls = {ALLEEG(1).chanlocs.labels};

%% Precomputation

% Equalize trial numbers across conditions
Ntrials = reshape([ALLEEG.trials],[length(cond),Ns]);
% Ntsub = min(Ntrials,[],1);
% kk = 1;
% for isub = 1:Ns
%   for icond = 1:length(cond)
%     it = round(1:ALLEEG(kk).trials/Ntsub(isub):ALLEEG(kk).trials);
%     ALLEEG(kk) = pop_select(ALLEEG(kk),'trial',it);
%     Ntrials(icond,isub) = ALLEEG(kk).trials;
%     kk = kk+1;
%   end
% end

% Filter the data
if strcmp(precomp.redo,'on')
  cnfg.maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
  cnfg.locutoff = 0.5;  % In Hz
  cnfg.hicutoff = 10;  % In Hz
  cnfg.transitionBandwidth = 1;  % In Hz
  KaiserWindowBeta = pop_kaiserbeta(cnfg.maxPassbandRipple);
  filtOrder = pop_firwsord('kaiser', cnfg.FS, cnfg.transitionBandwidth, cnfg.maxPassbandRipple);
  for ii = 1:length(ALLEEG)
    ALLEEG(ii) = pop_firws(ALLEEG(ii), 'fcutoff', [cnfg.locutoff,cnfg.hicutoff], 'wtype', 'kaiser',...
      'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);
  end
end

% Precomp settings
precomp.data = 'channels'; % 'components' or 'channels'
precomp.erp = 'on';
precomp.scalp = 'on';
precomp.spec = 'off';
precomp.ersp = 'off'; % 'on'
erpparams = { 'recompute',precomp.redo };
specparams = { 'freqrange', [2,25], 'timerange', [-1000 2500],'recompute',precomp.redo };
erspparams = { 'cycles', [2,10], 'freqs', [2 25],'nfreqs',20,...
  'alpha', 0.01, 'padratio', 1,'timelimits', [-1000 2500],'recompute',precomp.redo};
precompsettings = {'interp', 'on', 'recompute', precomp.redo,...
  'scalp', precomp.scalp, ...
  'erp', precomp.erp, 'erpparams', erpparams,...
  'spec', precomp.spec, 'specparams', specparams,...
  'ersp', precomp.ersp, 'itc', precomp.ersp, 'erspparams', erspparams};
[STUDY,ALLEEG] = std_precomp(STUDY, ALLEEG, precomp.data,precompsettings{:});

compLabl = {  'P1',   'N1',  'P2' };
compWin = {[30,100],[80,150],[150,250]};

%% ERPs
STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
      'timerange',[-600,1000],'averagechan','on','topotime',[]);
subsel = [subjects.name;{''}];
for ss=1:length(subsel)
  [STUDY,erpdata,erptimes,~,pcond] = std_erpplot(STUDY,ALLEEG,'channels',{'Cz'},'subject',subsel{ss});
  title(subsel{ss})
  % legend(designs{currentDesign-1,3},'Location','northwest')
  % Adjust line colors to following bar plot
  % colors = 0.9*flipud([53,42,134;55,184,157;248,250,13]/256);
  % chi = get(gca,'Children');
  % for ii = 1:length(chi)
  %   chi(ii).Color = colors(ii,:);
  % end
  % Add vertical lines to separate time windows
  hold on
  ylim = get(gca,'YLim');
  for ii = 2;%1:length(compWin)
    plot(compWin{ii}(1)*[1,1],ylim,'k:')
    plot(compWin{ii}(2)*[1,1],ylim,'k:')
  end

  ColorOrder = eye(3); ColorOrder=ColorOrder(ceil((1:6)/2),:);
  LineStyleOrder = repmat({'-','--'},1,3);
  chil = get(gca,'Children');
  chil = chil(end:-1:end-5);
  for ichil = 1:length(chil)
    set(chil(ichil),'Color',ColorOrder(ichil,:),'LineStyle',LineStyleOrder{ichil})
  end
  if isempty(subsel{ss})
    saveas(gcf,'SpExCue_Exp3_ERPwaveform_2ndSyl_avg')
  else
    saveas(gcf,['SpExCue_Exp3_ERPwaveform_2ndSyl_',subsel{ss}])
  end
end

%% Diff wave
figure
differpdata = cat(3,erpdata{1:2:end}) - cat(3,erpdata{2:2:end});
plot(erptimes,squeeze(mean(differpdata,2)))
xlabel('Time (ms)')
ylabel('Cz difference potential (µV)')
set(gca,'XLim',[-199,799])

%% ERP metrics

% peak-to-peak and N1
trange = [50,200];
idt = erptimes >= trange(1) & erptimes <= trange(2);
idtN1 = erptimes >= compWin{2}(1) & erptimes <= compWin{2}(2); 
CzP2P = nan(length(erpdata),Ns);
N1 = CzP2P;
kk = 1;
for isub = 1:Ns
  for icond = 1:length(cond)
    tmperp = mean(ALLEEG(kk).data(32,:,:),3);
    kk = kk+1;
    CzP2P(icond,isub) = range(tmperp(idt));
    N1(icond,isub) = min(tmperp(idtN1));%abs(mean(tmperp(idtN1)));
  end
end
disp(Ntrials)
CzP2PattenInc = CzP2P(1:3,:) - CzP2P(4:6,:);
N1attenInc = N1(1:3,:) - N1(4:6,:);
% for icond = 1:length(erpdata)
%   CzP2P(icond) = range(erpdata{icond}(idt));
% end
% CzP2PattenInc = CzP2P(1:2:end) - CzP2P(2:2:end);
figure;
% bar(CzP2PattenInc)
% plot(CzP2PattenInc)
% title('Peak-to-peak range 50-200 ms at Cz')
plot(N1attenInc)
title('N1 magnitude difference at Cz')
ylabel('Attended - unattended (µV)')
set(gca,'XTick',1:3,'XTickLabel',designs{currentDesign,3})
legend(subjects.name)

%% average potentials 
CzPeak = nan(length(erpdata),length(compWin));
CzLat = CzPeak;
for icomp = 1:length(compWin)
  idt = erptimes >= compWin{icomp}(1) & erptimes <= compWin{icomp}(2);
  tmp_erptimes = erptimes(idt);
  for icond = 1:length(erpdata)
    [CzPeak(icond,icomp),Ipeak] = max(abs(erpdata{icond}(idt)));
    CzLat(icond,icomp) = tmp_erptimes(Ipeak);
  end
end

% Bar plot for ERP metrics
figure
subplot(1,2,1)
bar(CzPeak')
ylabel('Peak amplitude (µV)')
set(gca,'XTickLabel',compLabl)
% legend(designs{currentDesign-1,3},'Location','northwest')
subplot(1,2,2)
bar(CzLat')
ylabel('Peak latency (ms)')
set(gca,'XTickLabel',compLabl)
if currentDesign == 1
  legend(STUDY.condition,'Location','northwest')
else
  legend(designs{currentDesign-1,3},'Location','northwest')
end

%% Inter-hemispheric differences
STUDY = pop_erpparams(STUDY, 'plotconditions', 'apart','averagechan','off');
[STUDY,erpdata,erptimes,~,pcond] = std_erpplot(STUDY,ALLEEG,'channels',ChanLbls(cnfg.eegChanNum));
chL = [1:12,14,15];
chR = 17:30;
rmsL = nan(length(erpdata),length(compWin));
rmsR = CzPeak;
for icomp = 1:length(compWin)
  idt = erptimes >= compWin{icomp}(1) & erptimes <= compWin{icomp}(2);
  for icond = 1:length(erpdata)
    rmsL(icond,icomp) = rms(max(abs(erpdata{icond}(idt,chL))));
    rmsR(icond,icomp) = rms(max(abs(erpdata{icond}(idt,chR))));
  end
end
RHD = rmsR-rmsL; % right-hemispheric dominance

% Bar plot
figure
if currentDesign == 4 % left/right
  IHAD = (rmsR-rmsL) ./ (rmsR+rmsL); % (contralateral - ipsilateral) / total
  IHAD(2,:) = -IHAD(2,:);
  bar(IHAD')
  ylabel('Inter-hemispheric contralaterality (µV)')
else
  bar(RHD')
  ylabel('Right-hemispheric dominance (µV)')
end
set(gca,'XTickLabel',compLabl)
legend(designs{currentDesign,3},'Location','northwest')

%% Topo plots
STUDY = pop_erpparams(STUDY, 'plotconditions', 'apart','averagechan','off');
compWin = {[40,80],[90,130],[140,220],[230,370]}; 
ylimTopo = 5*[-1,+1];
for ii = 1:length(compWin)
  [STUDY,erpdata] = std_erpplot(STUDY,ALLEEG,'channels',ChanLbls(cnfg.eegChanNum), ...
        'topotime',compWin{ii},'ylim',ylimTopo);
end
