clear all

if not(exist('pop_loadset','file'))
  addpath('/Users/rbaumgartner/Documents/MATLAB/eeglab/')
  eeglab
end

%% Settings
filepath = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Pilots/EEG3/data';
fnX = 'SpExCue_EEGpilot3_Sxx_blinkICrej_fmax30_Cxx.set';
fnStudy = 'SpExCue_EEGpilot3_ALL_blinkICrej_fmax30_groupE.study';
subject = {'RS','RB','S01','S02','S04'}; % RS: M = 0.5; S03: reversed externalization judgements
condition = {'Closer','Farther'};
design = 3; % 1: all, 2: Mc=1/3, 3: D=1/3, 4: D=-1/3

%% Create STUDY

commands = cell(1,length(subject)*length(condition));
fn = commands;
ii = 0;
for ss = 1:length(subject)
  for cc = 1:length(condition)
    ii = ii+1;
    fn{ii} = strrep(fnX,'Sxx',subject{ss});
    fn{ii} = strrep(fn{ii},'Cxx',condition{cc});
    commands{ii} = { 'index' ii 'load' fullfile(filepath,fn{ii}),...
      'subject' subject{ss} 'condition' condition{cc} };
  end
end

% Save
[STUDY,ALLEEG] = std_editset( [], [], 'commands', commands,...
  'task', 'closer vs. farther for second of two consecutive noise bursts',...
  'name', 'SpExCue', 'filename', fnStudy,'filepath',filepath); 

% Create additional designs
STUDY = std_makedesign(STUDY, ALLEEG, 2, 'name', 'Mc0p3',...
  'variable1', 'condition', 'datselect',{'type',{12,32}},'subjselect',subject);
STUDY = std_makedesign(STUDY, ALLEEG, 3, 'name', 'D0p3',...
  'variable1', 'condition', 'datselect',{'type',{12}},'subjselect',subject(2:end));
STUDY = std_makedesign(STUDY, ALLEEG, 4, 'name', 'Dm0p3',...
  'variable1', 'condition', 'datselect',{'type',{21}},'subjselect',subject(2:end));
STUDY.currentdesign = design;

% update workspace variables and redraw EEGLAB
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
eeglab redraw

%% Compute measures

% Number of trials per subject and condition
Ntrials = zeros(length(condition),length(subject)+1);
for ii = 1:length(STUDY.design(design).cell)
  ss = ismember(subject,STUDY.design(design).cell(ii).case);
  cc = ismember(condition,STUDY.design(design).cell(ii).value{1});
  Ntrials(cc,ss) = length(STUDY.design(design).cell(ii).trials{1});
end
Ntrials(:,end) = sum(Ntrials(:,1:length(subject)),2);
NtrialsTab = array2table(Ntrials,'RowNames',condition,'VariableNames',[subject,{'All'}]);

% Channels
[STUDY ALLEEG customRes] = std_precomp(STUDY, ALLEEG, 'channels', 'interp', 'on',... % { 'Cz' 'Pz' }
	'erp', 'on', 'spec', 'off', 'ersp', 'off', 'itc', 'off', 'erspparams', ...
    { 'cycles', [2,10], 'freqs', [2 20],'nfreqs',10, 'alpha', 0.01, 'padratio' 1 });
% ICs
% [STUDY ALLEEG customRes] = std_precomp(STUDY, ALLEEG, 'components', 'interp', 'on',...
% 	'erp', 'on', 'spec', 'off', 'ersp', 'off', 'itc', 'off', 'erspparams', ...
%     { 'cycles', [2,10], 'freqs', [2 20],'nfreqs',10, 'alpha', 0.01, 'padratio' 1 });

% [STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, [], ... 
% 	{'spec' 'npca' 10 'norm' 1 'weight' 1 'freqrange' [2 20]} , ... 
% 	{'erp' 'npca' 10 'norm' 1 'weight' 1 'timewindow' [-200 800]} , ... 
% 	{'dipoles' 'norm' 1 'weight' 10} , ...
%   {'ersp' 'npca' 10 'freqrange' [2 20] 'cycles' [2 10] 'freqs' [2 20]...
% 		 'alpha' 0.01 'padratio' 4 'timewindow' [-999 1999] 'norm' 1 'weight' 1});

% [STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, [],...
%   {'erp' 'npca' 10 'norm' 1 'weight' 1 'timewindow' [-200 800]},...
%   {'scalp' 'npca' 10 'norm' 1 'weight' 1 'abso' 0 });
% STUDY = pop_clust( STUDY, ALLEEG, 'algorithm', 'kmeanscluster');   

% [STUDY, ALLEEG] = pop_precomp(STUDY, ALLEEG); % precomputed everything
              
%% Topographic analysis 
STUDY = pop_statparams(STUDY,'condstats','on','mode','eeglab','mcorrect','none','alpha',nan);
% STUDY = pop_statparams(STUDY,'condstats','on','mode','fieldtrip',...
%   'fieldtripmethod','montecarlo','fieldtripmcorrect','cluster');% 'fieldtripclusterparam',{'clusterstatistic','maxsum'}


% Scalp topographies
ChanLbls = {ALLEEG(1).chanlocs.labels};
topotime = {[100 150];[200 250];[300 600]};
% topotime = {[200 250]};
for ii = 1:length(topotime)
  [STUDY, C(ii).erpdata, ~, ~, C(ii).pcond] = std_erpplot(STUDY,ALLEEG,...
    'channels',ChanLbls(1:32), 'topotime',topotime{ii},'ylim',[-4,4]);
end

%% ERP analysis
STUDY = pop_erpparams(STUDY, 'plotconditions', 'together',...
  'timerange',[-200,800],'averagechan','on','topotime',[]);
alpha = .05;
for ii = 1%:length(C)
  if any(C(ii).pcond{1} < alpha)
    selChan = find(C(ii).pcond{1} < alpha)';
  else
    [~,selChan] = min(C(ii).pcond{1});
  end
  for jj = [23,32];%selChan
    [STUDY, C(ii).erp(jj).erpdata, C(ii).erp(jj).erptimes, ~, C(ii).erp(jj).pcond] = std_erpplot(...
      STUDY,ALLEEG,'channels',ChanLbls(jj));

    % BUG: EEGlab does not plot wave forms for fixed alpha in statparams ->
    % workaround with alpha = nan in statparams:
    tsigdiff = C(ii).erp(jj).erptimes(C(ii).erp(jj).pcond{1} < alpha); % sig. diff. time instances 
    delete(subplot(122));
    pos = get(subplot(121),'Position');
    set(subplot(121),'Position',pos.*[1,1,2,1])
    YLim = get(gca,'YLim');
    hold on
    plot(tsigdiff,(YLim(1)+.1)*ones(length(tsigdiff),1),'r.')
  end
end

% STUDY = pop_chanplot(STUDY,ALLEEG);

%% Print

disp(NtrialsTab)

FontSize = 8;
Resolution = '-r600';
set(findall(gcf,'-property','FontSize'),'FontSize',FontSize)

fn = ['analyzeEEGpilot3_eeglab_study/topo_','N1','_groupE_' STUDY.design(design).name];
set(gcf,'PaperUnits','centimeters','PaperPosition',[100,100,20,5])

% fn = ['analyzeEEGpilot3_eeglab_study/erp_','Cz','_groupE_' STUDY.design(design).name];
% set(gcf,'PaperUnits','centimeters','PaperPosition',[100,100,10,6])

print(gcf,Resolution,'-dpng',fn)