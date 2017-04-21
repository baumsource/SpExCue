function [stim,Obj] = SpExCue_stim( M,ID,varargin )
%SpExCue_stim creates stimuli for SpExCue project
%   Usage: [stim,Obj] = SpExCue_stim( M, ID, [pos, fs, flow, fhigh, SPL] )
%
%   Input parameters:
%     M     : spectral magnitude compression factor.
%     ID    : listener's ID
%
%   Output parameters:
%     stim :  stimulus structure. field *sig*: sound signals as cell array 
%             with dimensions M x pos; field *fs*: sampling rate.
%     Obj  :  listener-specific reference HRTFs in SOFA format
%
%   Optional key-value pairs:
%     pos   : set of positions. 1st column: azimuth, 2nd column: elevation.
%             Default is [0,0].
%     fs    : sampling rate. Default is 48kHz.
%     flow  : lower cut-off frequency. Default is 700 Hz.
%     fhigh : higher cut-off frequency. Default is 18 kHz.
%     SPL   : sound pressure level in dB (RMS of 1 represents 100 dB).
%             Default is 70 dB.
%     individualSignalDur  : signal duration in seconds. Default is 1.6s.
%     fadeDuration    : fade in/out duration in seconds. Default is 0.05s.
%
%   Optional flags:
%     noiseBurst   : Gaussian white noise burst.
%     speech  : Speech syllable.
%     DTF     : Directional transfer functions.
%     noBP    : no band-pass filtering.

% AUTHOR: Robert Baumgartner

%% Check input

% definput.keyvals.flow=700; % Hz
% definput.keyvals.fhigh=18e3; % Hz
% definput.keyvals.individualSignalDur=1.6; % seconds
% definput.keyvals.SPL=70; % dB
% definput.keyvals.fs=48e3; % Hz
% definput.keyvals.fadeDur = 0.05; % duration of fade-in/out in seconds
% definput.keyvals.pos=[0,0]; % spatial position [azimuth,elevation] in deg.

% definput.flags.HRTFs = {'HRC3ms','HRCeq','ARI'};
% definput.flags.DTF = {'','DTF'};
% definput.flags.movement = {'static','moveLeftRight'};
% definput.flags.source = {'continuousNoise','noiseBurst','speech','AMnoiseBurst','AMnoise','IR'};

definput = arg_SpExcue;
definput.keyvals.pos=[0,0]; % spatial position [azimuth,elevation] in deg.
[flags,kv]=ltfatarghelper(...
  {'pos','fs','flow','fhigh','SPL','individualSignalDur'},definput,varargin);

% position key-value pair kept for backwards compatibility
if all(kv.pos(:) == 0) && ismember('azi',varargin(cellfun('isclass',varargin,'char')))
  if length(kv.ele) == 1 && length(kv.azi) > 1
    kv.ele = repmat(kv.ele,length(kv.azi),1);
  end
  kv.pos = [kv.azi(:),kv.ele(:)];
end

%% Listener's HRTFs
HRTFpath = strrep(which('SpExCue_stim'),...
  fullfile('MATLAB_general','SpExCue_stim.m'),'HRTFs');
if flags.do_ARI
  refObj = SOFAload(fullfile(HRTFpath,[ID '_hrtf B.sofa']));
elseif flags.do_HRC3ms
  refObj = SOFAload(fullfile(HRTFpath,[ID '_3ms.sofa']));
else
  refObj = SOFAload(fullfile(HRTFpath,[ID '_eq.sofa']));
end
fsHRTF = refObj.Data.SamplingRate;

%% Optional DTF conversion
if flags.do_DTF
  refObj = SOFAhrtf2dtf(refObj);
end

%% Source signal
sigPath = strrep(HRTFpath,'HRTFs','Signals');
if flags.do_noiseBurst
  % ramped 500-ms Gaussian white noise burst
  sig = noise(round(0.5*fsHRTF),1,'white');
elseif flags.do_AMnoise % 100-Hz modulation frequency (Getzmann & Lewald, 2010)
  sig = noise(round(kv.individualSignalDur*fsHRTF),1,'white');
  t = 0:1/fsHRTF:kv.individualSignalDur-1/fsHRTF;
  AMf = 100; % modulation frequency
  MD = .12; % modulation depth
  AM = 1+MD*sin(2*pi*AMf*t'+pi);
  sig = sig.*AM;
elseif flags.do_AMnoiseBurst
  % 500-ms Gaussian white noise burst with AM
  sig = noise(round(0.5*fsHRTF),1,'white');
  t = 0:1/fsHRTF:0.5-1/fsHRTF;
  AMf = 4;
  sig = sig.*(1+0.5.*cos(2*pi*AMf*t'+pi));
elseif flags.do_speech
  tmp = load(fullfile(sigPath,'scare'));
  sig = tmp.scare(:);
  fssig = tmp.fs; % Hz
  % concatenate (commented because preprocessed) 
%   nstart = 15460;
%   sig = sig(nstart:min(nstart+round(kv.individualSignalDur*fssig),length(sig)));
  % resample acc. to HRTFs
  Q = fssig/gcd(fssig,fsHRTF);
  P = fsHRTF/gcd(fssig,fsHRTF);
  sig = resample(sig,P,Q);
% elseif flags.do_IR
%   sig = [1;zeros(round(kv.individualSignalDur*fsHRTF),1)];
else % flags.do_continuousNoise
  % long Gaussian white noise burst
  sig = noise(round(kv.individualSignalDur*fsHRTF),1,'white'); 
end

%% Fade in/out
sinRamp = sin(pi/2*(0:kv.fadeDur*fsHRTF-1)/(kv.fadeDur*fsHRTF)).^2;
sig = sig.*[sinRamp,ones(1,length(sig)-2*length(sinRamp)),fliplr(sinRamp)]';

%% HRTF filtering, M adjustment  
for ii = 1:length(M)
  Obj{ii} = SpExCue_setMagVar(refObj,M(ii),kv.flow,kv.fhigh);
  if flags.do_IR
    for jj = 1:size(kv.pos,1)
      idpos = Obj{ii}.SourcePosition(:,1) == kv.pos(jj,1) & Obj{ii}.SourcePosition(:,2) == kv.pos(jj,2);
      stim.sig{ii,jj} = shiftdim(Obj{ii}.Data.IR(idpos,:,:),2);
    end
  else
    if flags.do_static
      for jj = 1:size(kv.pos,1)
        stim.sig{ii,jj} = SpExCue_SOFAspat(sig,Obj{ii},kv.pos(jj,1),kv.pos(jj,2));
      end
    elseif flags.do_moveLeftRight    
      azi = kv.pos(1,1);
      ele = kv.pos(1,2);
      leftwards = [repmat(azi,[1,7]),azi+(0:90/5:90)];
      rightwards = [repmat(azi,[1,7]),azi-(0:90/5:90)];
      [stim.sig{ii,1}, tmp.azi, tmp.ele] = SpExCue_SOFAspat(...
        sig,Obj{ii},leftwards,ele*ones(1,length(leftwards)));
  %     figure; plot(0:kv.individualSignalDur/(length(tmp.azi)-1):kv.individualSignalDur,tmp.azi)
      stim.sig{ii,2} = SpExCue_SOFAspat(sig,Obj{ii},rightwards,ele*ones(1,length(rightwards)));
      stim.sig{ii,3} = SpExCue_SOFAspat(sig,Obj{ii},azi,ele);
    end
  end
end

%% Resample stim acc. to fs
if not(isequal(kv.fs,fsHRTF))
  P = kv.fs/gcd(kv.fs,fsHRTF);
  Q = fsHRTF/gcd(kv.fs,fsHRTF);
  for ii = 1:numel(stim.sig)
    left = resample(stim.sig{ii}(:,1),P,Q);
    right = resample(stim.sig{ii}(:,2),P,Q);
    stim.sig{ii} = [left(:),right(:)];
  end
end
stim.fs = kv.fs;

%% Band-pass filtering
% if not(flags.do_IR)
%   if kv.flow <= 100
%     filtOrder = 6;
%   else
%     filtOrder = 9;
%   end
if flags.do_butter
  filtOrder = 4;
  [b_bp,a_bp] = butter(filtOrder,[kv.flow,kv.fhigh]/(kv.fs/2));
  for ii = 1:numel(stim.sig)
    stim.sig{ii} = filter(b_bp,a_bp,stim.sig{ii});
  end
end
% else
  
% end

%% Set level (adjust all stimuli by the same factor)
currentSPL = nan(numel(stim.sig),2);
for ii = 1:numel(stim.sig)
  currentSPL(ii,:) = dbspl(stim.sig{ii});
end

adjustmentFactor = 10^((kv.SPL-mean(currentSPL(:)))/20);
for ii = 1:numel(stim.sig)
  stim.sig{ii} = adjustmentFactor*stim.sig{ii};
  if max(stim.sig{ii}(:)) > 1; 
    error(['Stimulus clips at present SPL of ',num2str(kv.SPL),' dB. Reduce SPL and try again!'])
  end
end

% 2016/04/26: correction factor for center stimuli caused SPL changes
% for ii = 1:numel(stim.sig)
%   ILD(ii) = -diff(dbspl(stim.sig{ii}));
%   if ILD(ii) > 0
%     stim.sig{ii}(:,1) = setdbspl(stim.sig{ii}(:,1),kv.SPL);
%     stim.sig{ii}(:,2) = setdbspl(stim.sig{ii}(:,2),kv.SPL-ILD(ii));
%   else % ILD <= 0
%     stim.sig{ii}(:,1) = setdbspl(stim.sig{ii}(:,1),kv.SPL+ILD(ii));
%     stim.sig{ii}(:,2) = setdbspl(stim.sig{ii}(:,2),kv.SPL);
%   end
%   if abs(ILD(ii)) < 1 % attenuate for approximately equal binaural loudness in midsagittal plane 
%     stim.sig{ii} = 0.6*stim.sig{ii};
%   end
% end

if nargout == 2
  Obj = refObj;
end

end