function eeglab_dataSplit(cnfg,fn,ID)

EEG = pop_loadset('filename', fn);
Ntrials = zeros(1,size(cnfg.cond,1));
EEGmem = cell(1,size(cnfg.cond,1));
for cc = 1:size(cnfg.cond,1)
  EEGsel = pop_epoch(EEG, cnfg.cond{cc,1}, cnfg.finalEpochRange);
  EEGsel = pop_rmbase(EEGsel, cnfg.baselineRange);
  EEGsel.subject = ID;
  EEGsel.condition = cnfg.cond{cc,2};
  Ntrials(cc) = EEGsel.trials;
  EEGmem{cc} = EEGsel;
end

% Equalize trial numbers
NtrialsEq = min(Ntrials);
for cc = 1:size(cnfg.cond,1)
  it = round(1:Ntrials(cc)/NtrialsEq:Ntrials(cc)); % equally spaced removal of trials
  EEGsel = pop_select(EEGmem{cc},'trial',it);
  pop_saveset(EEGsel,'filename', [fn(1:end-4),'_',cnfg.cond{cc,2},'.set']);
end