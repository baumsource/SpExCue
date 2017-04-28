% Script to run experimental procedure for syllable segregation with
% different spatial cues: ITD-only, ILD-only, individualized HRTFs

%% Individual setting
subID='S15'; 

% un/comment as appropriate
session = 'train'; % training session, 1 block (~ 7 minutes)
% session = 'test'; % testing session (EEG-monitored) with 9 blocks (takes about 1h)

%% General settings
flags.do_debugMode = 0; % set to 1 in order to show subject screen only in a small window 
flags.do_bonus = 1;     % set to 1 in order to include bonus payment of 1 cent per correct trial
flags.do_feedback = 1;  % set to 1 in order to display feedback screen
kv.screenNumber = 1; % choose screen number for psychtoolbox

% Screen('Preference', 'SkipSyncTests', 1); % uncomment if problems with psychtoolbox sync test occur (visual timing is not critical for this experiment)

%% Start Experiment
switch session
  case 'train'
    % training session, 1 block (~ 7 minutes)
    SpExCue_BaDaGa(subID,2,'train',flags,kv)
  case 'test'
    % testing session (EEG-monitored) with 9 blocks (takes about 1h)
    SpExCue_BaDaGa(subID,18,'test',flags,kv)
end