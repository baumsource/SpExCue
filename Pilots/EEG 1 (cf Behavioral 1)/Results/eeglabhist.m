% EEGLAB history file generated on the 07-Mar-2016
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename','SpExCue_EEGpilot_RB_resamp.set','filepath','/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Experiments/Pilot/Results/');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = pop_eegfiltnew(EEG, 0.5, 20, 1690, 0, [], 1);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew','/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Experiments/Pilot/Results/SpExCue_EEGpilot_RB_resamp_LP20Hz.set','overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_reref( EEG, [33 34] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew','/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Results/EEGpilot/RB_reref.set','overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_select( EEG,'nochannel',{'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8'});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew','/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Results/EEGpilot/RB_reref.set','gui','off'); 
EEG=pop_chanedit(EEG, 'load',{'/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Tools/EEG/biosemi_eloc.locs' 'filetype' 'loc'});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'extended',1,'interupt','on');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename','SpExCue_EEGpilot_RB_resamp_LP20Hz_ICAraw.set','filepath','/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Experiments/Pilot/Results/');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;
