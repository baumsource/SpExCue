function definput=arg_SpExcue(definput)

definput.keyvals.fnExtension = ''; % filename extension
definput.keyvals.azi = [-90,0,90]; % set of azimuths
definput.keyvals.ele = 0; % set of elevations
definput.keyvals.M = [1,1/3,0]; % set of spectral saliences
definput.keyvals.Nrep = 120; % repetitions
definput.keyvals.flow = 1e3; % lower cut-off frequency
definput.keyvals.fhigh = 16e3; % upper cut-off frequency
definput.keyvals.fs=48e3; % Hz
definput.keyvals.dur = 1.2; % duration of stimulus pair in sec
definput.keyvals.fadeDur = 0.05; % duration of fade-in/out in sec
definput.keyvals.xFadeDur = 0.01; % duration of cross-fade in sec
definput.keyvals.individualSignalDur=1.6; % seconds
definput.keyvals.jitter = 0.1; % temporal jitter in sec
definput.keyvals.SPL = 70; % SPL in dB
definput.keyvals.SPLrove = 0; % roving of SPL in dB
definput.keyvals.rdepth = 10; % depth (peak-to-peak) in dB of spectral magnitude ripple
definput.keyvals.rdensity = 0.25; % density of spectral magnitude ripple
definput.keyvals.screenNumber = 1; % For 3rd floor lab (when using dual screens) choose #1, otherwise use #0
definput.keyvals.NperBlock = 21; % trials per block
definput.flags.HRTFs = {'HRCeq','HRC3ms','ARI'};
definput.flags.DTF = {'','DTF'};
definput.flags.movement = {'static','moveLeftRight'};
definput.flags.source = {'continuousNoise','noiseBurst','speech','AMnoiseBurst','AMnoise','IR'};
definput.flags.stimulation = {'binaural','monaural'};
definput.flags.Mrepetition = {'repeateM','changeM'};
definput.flags.positionOrderPermutation = {'completePosOrderPermutation','halfPosOrderPermutation'}; % randomized permutation of sequential order of positions either tested completely or only the first half
definput.flags.roving = {'componentRove','pairRove','noRoving'};
definput.flags.bandpass = {'butter','noBP'};
definput.flags.familiarization = {'skipFamiliarization','familiarize'};
definput.flags.feedback = {'blockedFeedback','TbTfeedback','noFeedback','consistencyFeedback','D0detectionFeedback'}; % Give feedback on trial-by-trial basis by changing fixation dot color
definput.flags.tdt = {'TDTon','TDToff'};
definput.flags.psychtoolbox = {'','debugMode'}; % to generate a very small (useless) Psychtoolbox window but allow access to Matlab command window while executing script

end