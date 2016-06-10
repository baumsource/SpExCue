% analyze SpExCue_EEGpilot with fieldtrip


cfg.dataset      = '/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Pilots/EEG 1 (cf Behavioral 1)/Results/SpExCue_EEGpilot_DS.bdf';%string with the filename
% cfg.trl          = %Nx3 matrix with the trial definition, see FT_DEFINETRIAL
% cfg.padding      = %length (in seconds) to which the trials are padded for filtering (default = 0)
% cfg.padtype      = %string, type of padding (default: 'data' padding or 'mirror', depending on feasibility)
% cfg.continuous   = %'yes' or 'no' whether the file contains continuous data (default is determined automatic)

cfg.channel      = {'all','-EXG7','-EXG8'};%Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details

cfg.bpfilter      = 'yes';%'no' or 'yes'  bandpass filter (default = 'no')
cfg.bpfreq        = [.5,25];%bandpass frequency range, specified as [low high] in Hz

cfg.reref         = 'yes'; %'no' or 'yes' (default = 'no')
cfg.refchannel    = {'EXG1','EXG2'};%cell-array with new EEG reference channel(s), this can be 'all' for a common average reference

t = struct(... 
  'correctResponse',8,...
  'wrongResponse',9,...
  'left',000,...
  'top',100,...
  'right',200,...
  'M1_1',30,... % M1 (at stimulus onset) of 1
  'M1_p5',20,...
  'M1_0',10,...
  'M2_1',3,... % M2 (at stimulus change) of 1
  'M2_p5',2,...
  'M2_0',1,...
  'stimulusOffset',9);
trgOnset = [t.M1_1,t.M1_p5,t.M1_0]; % event triggers (as a row vector)
% trgOnset = cat(2,[trgOnset+trigVals.top,trgOnset+trigVals.right]); % only left direction tested
% trgOnsetLabel = {'M: 0','M: .5','M: 1'};
trgChange = [t.M1_1  + t.M2_0 , t.M1_1  + t.M2_p5,...
             t.M1_p5 + t.M2_0 , t.M1_p5 + t.M2_1...
             t.M1_0  + t.M2_p5, t.M1_0  + t.M2_1];

cfg.trialfun   =  'ft_trialfun_SpExCue_EEGpilot';
cfg.trialdef.eventtype      = 'M1_0';
cfg.trialdef.eventvalue     = 10; % the value of the stimulus trigger for fully incongruent (FIC).
cfg.trialdef.prestim        = .2;
cfg.trialdef.poststim       = .9;
cfg         = ft_definetrial(cfg);

[data] = ft_preprocessing(cfg);
