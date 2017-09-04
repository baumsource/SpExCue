ERPanalysis = true;

redoStudy = true;

currentDesign = 1;

precomp.redo = 'off';

designs = {...
  'attended',{'ITD','ILD','HRTF'},{'ITD','ILD','HRTF'};...
  'unattended',{'ITDu','ILDu','HRTFu'},{'ITDu','ILDu','HRTFu'};...
  'attendedVsUnattended',{{'ITD','ILD','HRTF'},{'ITDu','ILDu','HRTFu'}},{'attended','unattended'};...
};

%% Load/create STUDY
eeglab

experimentString = 'Exp3eeg';
datapath = fullfile('..','data');%string with the filename

% subjects
tmp = load('SpExCue_Exp3eeg_subjects.mat');
subjects = tmp.subject;
% subjects = subjects(2:4,:); disp('only ...')
Ns = height(subjects);

% configurations
tmp = load([experimentString,'_cnfg']);
cnfg = tmp.cnfg;
cond = cnfg.cond;

fnstudy = ['study_',experimentString,'.study'];
if ERPanalysis; fnstudy = strrep(fnstudy,'.','_ERP.'); end
if exist(fullfile(datapath,fnstudy),'file') && not(redoStudy)
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


%% Retrieve component windows from grand average GFP
STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
      'timerange',[-600,1000],'averagechan','off','topotime',[]);
[STUDY,erpdata,erptimes] = std_erpplot(STUDY,ALLEEG,'noplot','on','channels',ChanLbls(1:32));
grandaverage = mean(cat(3,erpdata{:}),3);
gavgGFP = std(grandaverage,0,2);
figure
XLim = [-199,799];
plot(erptimes,gavgGFP,'k')
xlabel('Time (ms)')
ylabel('Grand average global field potential (µV)')
set(gca,'XLim',XLim)

% grand-average GFP shows a distinct peak in the N1 range 
% ADJUST:
compLabl = {'N1'};
compWin = {[100,180]};

% check by plotting on top of GFP
hold on
YLim = get(gca,'YLim');
for icomp = 1:length(compWin)
  plot(compWin{icomp}([1,1]),YLim,'k:')
  plot(compWin{icomp}(1+[1,1]),YLim,'k:')
  text(mean(compWin{icomp}),YLim(2)-0.05,compLabl{icomp},'HorizontalAlignment','center')
end


%% Topo plots
STUDY = pop_erpparams(STUDY, 'plotconditions', 'apart','averagechan','off');
ylimTopo = 3*[-1,+1];
for ii = 1:length(compWin)
  [STUDY,erpdata] = std_erpplot(STUDY,ALLEEG,'channels',ChanLbls(cnfg.eegChanNum), ...
        'topotime',compWin{ii},'ylim',ylimTopo);
end

% ADJUST channel of interest cluster:
chanSel = {'Pz'};
disp(['Channel(s) of interest: ', chanSel])

%% ERPs
STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
      'timerange',[-600,1000],'averagechan','on','topotime',[]);
subsel = [subjects.name;{''}];
for ss=1:length(subsel)
  [STUDY,erpdata,erptimes,~,pcond] = std_erpplot(STUDY,ALLEEG,'channels',chanSel,'subject',subsel{ss});
  ylabel([chanSel{:},' potential (µV)'])
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
  for ii = 1:length(compWin)
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
ylabel([chanSel{:},' difference potential (µV)'])
legend(STUDY.condition{1:2:end})
set(gca,'XLim',[-199,799])

%% ERP metrics

% peak-to-peak difference
% trange = [50,200];
% idt = erptimes >= trange(1) & erptimes <= trange(2);
% CzP2P = nan(length(erpdata),Ns);
% kk = 1;
% for isub = 1:Ns
%   for icond = 1:length(cond)
%     tmperp = mean(ALLEEG(kk).data(32,:,:),3);
%     kk = kk+1;
%     CzP2P(icond,isub) = range(tmperp(idt));
%   end
% end
% CzP2PattenInc = CzP2P(1:3,:) - CzP2P(4:6,:);
% plot(CzP2PattenInc)
% title('Peak-to-peak range 50-200 ms at Cz')

% N1 difference
% icomp = ismember(compLabl,'N1');
% idtN1 = erptimes >= compWin{icomp}(1) & erptimes <= compWin{icomp}(2); 
% N1 = nan(length(erpdata),Ns);
% kk = 1;
% chanNum = ismember(ChanLbls,chanSel);
% for isub = 1:Ns
%   for icond = 1:length(cond)
%     tmperp = mean(ALLEEG(kk).data(chanNum,:,:),3);
%     kk = kk+1;
%     N1(icond,isub) = min(tmperp(idtN1));%abs(mean(tmperp(idtN1)));
%   end
% end
% disp(Ntrials)
% N1attenInc = N1(1:3,:) - N1(4:6,:);
% figure;
% plot(N1attenInc)
% title(['N1 magnitude difference at ',chanSel{:}])
% ylabel('Attended - unattended (µV)')
% set(gca,'XTick',1:3,'XTickLabel',{ALLEEG(1:3).condition})
% legend(subjects.name)

%% average potentials 
STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
      'timerange',[-600,1000],'averagechan','on','topotime',[]);
[STUDY,erpdata,erptimes] = std_erpplot(STUDY,ALLEEG,'channels',chanSel,'noplot','on');
icomp = ismember(compLabl,'N1');
idt = erptimes >= compWin{icomp}(1) & erptimes <= compWin{icomp}(2);
tmp_erptimes = erptimes(idt);
CzPeak = nan(length(erpdata),Ns);
CzLat = CzPeak;
for icond = 1:length(erpdata)
  [CzPeak(icond,:),Ipeak] = max(abs(erpdata{icond}(idt,:)));
  CzLat(icond,:) = tmp_erptimes(Ipeak);
end

% Bar plot for ERP metrics
figure
subplot(1,2,1)
bar(CzPeak')
ylabel('Peak magnitude (µV)')
set(gca,'XTickLabel',subjects.name)
% legend(designs{currentDesign-1,3},'Location','northwest')
subplot(1,2,2)
bar(CzLat')
ylabel('Peak latency (ms)')
set(gca,'XTickLabel',subjects.name)
legend(STUDY.condition,'Location','northwest')


%% Inter-hemispheric differences
STUDY = pop_erpparams(STUDY, 'plotconditions', 'apart','averagechan','off');
[STUDY,erpdata,erptimes,~,pcond] = std_erpplot(STUDY,ALLEEG,'channels',ChanLbls(cnfg.eegChanNum),'noplot','on');
chL = [1:12,14,15];
chR = 17:30;
rmsL = nan(length(erpdata),length(compWin),Ns);
rmsR = rmsL;
for icomp = 1:length(compWin)
  idt = erptimes >= compWin{icomp}(1) & erptimes <= compWin{icomp}(2);
  for icond = 1:length(erpdata)
    rmsL(icond,icomp,:) = rms(max(abs(erpdata{icond}(idt,chL,:))));
    rmsR(icond,icomp,:) = rms(max(abs(erpdata{icond}(idt,chR,:))));
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
  bar(squeeze(RHD)')
  ylabel('Right-hemispheric dominance (µV)')
end
set(gca,'XTickLabel',subjects.name)
legend(STUDY.condition,'Location','northwest')
