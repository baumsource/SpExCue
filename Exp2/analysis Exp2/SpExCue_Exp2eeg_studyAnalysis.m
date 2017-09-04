currentDesign = 3;
% 1: all
% 2: 'onset'
% 3: motion
% 4: L/R


precomp.redo = 'off';

designs = {...
  'onset',{'C0_center','C1_center','KEMAR_center'},{'C = 0','C = 1','KEMAR'};...
  'motion',{{'C0_left','C0_right'},{'C1_left','C1_right'},{'KEMAR_left','KEMAR_right'}},{'C = 0','C = 1','KEMAR'};...
  'LR',{{'C0_left','C1_left','KEMAR_left'},{'C0_right','C1_right','KEMAR_right'}},{'left','right'}};

%% Load/create STUDY
datapath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Exp2/data/';%string with the filename

if exist('study_Exp2eeg','file')
  [STUDY, ALLEEG] = pop_loadstudy('filename', fullfile(datapath,'study_Exp2eeg'));
else
  fnX = 'Exp2eeg_Sxx_clean_Cxx.set';
  commands = cell(1,height(subjects)*size(cond,1));
  fn = commands;
  ii = 0;
  for ss = 1:height(subjects)
    for cc = 1:size(cond,1)
      ii = ii+1;
      fn{ii} = strrep(fnX,'Sxx',subjects(ss,:).name);
      fn{ii} = strrep(fn{ii},'Cxx',cond{cc,2});
      commands{ii} = { 'index' ii 'load' fullfile(datapath,fn{ii}),...
        'subject' subjects(ss,:).name 'condition' cond{cc,2}};
    end
  end

  % Create STUDY
  [STUDY,ALLEEG] = std_editset([],[],'commands',commands,'name','SpExCue Exp. 2',...
    'task', 'left- vs. rightward motion discrimination'); 
  % Define designs
  STUDY.design(1).name = 'all';
  for dd = 1:length(designs)
    STUDY = std_makedesign(STUDY, ALLEEG, dd+1, 'name',designs{dd,1},...
      'variable1', 'condition', 'values1',designs{dd,2}); %,'subjselect',subject
  end
  % Save
  [STUDY, ALLEEG] = std_editset(STUDY, ALLEEG,'filename', fullfile(datapath,'study_Exp2eeg'));
end
%%
[STUDY] = std_selectdesign(STUDY, ALLEEG, currentDesign);

ChanLbls = {ALLEEG(1).chanlocs.labels};

%% Precomputation
% Precomp settings
precomp.data = 'channels'; % 'components' or 'channels'
precomp.erp = 'on';
precomp.scalp = 'on';
precomp.spec = 'off';
precomp.ersp = 'off'; % 'on'
erpparams = { 'recompute',precomp.redo };
specparams = { 'freqrange', [2,20], 'timerange', [-200 800],'recompute',precomp.redo };
erspparams = { 'cycles', [2,10], 'freqs', [2 20],'nfreqs',10,...
  'alpha', 0.01, 'padratio' 1,'timelimits', [-200 800],'recompute',precomp.redo};
precompsettings = {'interp', 'on', 'recompute', precomp.redo,...
  'scalp', precomp.scalp, ...
  'erp', precomp.erp, 'erpparams', erpparams,...
  'spec', precomp.spec, 'specparams', specparams,...
  'ersp', precomp.ersp, 'itc', precomp.ersp, 'erspparams', erspparams};
[STUDY,ALLEEG] = std_precomp(STUDY, ALLEEG, precomp.data,precompsettings{:});

compLabl = {  'P1',   'N1',  'P2',    'P3'};
compWin = {[40,80],[90,130],[140,220],[230,370]};

%% ERPs
STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
      'timerange',[-200,600],'averagechan','on','topotime',[]);
[STUDY,erpdata,erptimes,~,pcond] = std_erpplot(STUDY,ALLEEG,'channels',{'Cz'});
legend(designs{currentDesign-1,3},'Location','northwest')

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

% Bar plot
figure
subplot(1,2,1)
bar(CzPeak')
ylabel('Peak amplitude (µV)')
set(gca,'XTickLabel',compLabl)
legend(designs{currentDesign-1,3},'Location','northwest')
subplot(1,2,2)
bar(CzLat')
ylabel('Peak latency (ms)')
set(gca,'XTickLabel',compLabl)

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
legend(designs{currentDesign-1,3},'Location','northwest')

%% Topo plots
STUDY = pop_erpparams(STUDY, 'plotconditions', 'apart','averagechan','off');
compWin = {[40,80],[90,130],[140,220],[230,370]}; 
ylimTopo = 5*[-1,+1];
for ii = 1:length(compWin)
  [STUDY,erpdata] = std_erpplot(STUDY,ALLEEG,'channels',ChanLbls(cnfg.eegChanNum), ...
        'topotime',compWin{ii},'ylim',ylimTopo);
end
