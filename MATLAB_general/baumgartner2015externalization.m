function [dext,tang,rang] = baumgartner2015externalization( target,template,varargin )
%BAUMGARTNER2015EXTERNALIZATION Model for sound externalization
%   Usage:    [dext,tang,rang] = baumgartner2015externalization( target,template )
%
%   Input parameters:
%     target  : binaural impulse response(s) referring to the directional 
%               transfer function(s) (DFTs) of the target sound(s).
%               Option 1: given in SOFA format -> sagittal plane DTFs will 
%               be extracted internally. 
%               Option 2: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane formated 
%               according to the following matrix dimensions: 
%               time x direction x channel/ear
%     template: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane referring to
%               the perceived lateral angle of the target sound.
%               Options 1 & 2 equivalent to *target*.
%
%   Output parameters:
%     dext    : predicted degree of externalization
   
% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

%% Check input

definput.import={'baumgartner2014'};

[flags,kv]=ltfatarghelper(...
  {'fs','S','lat','stim','fsstim','space','do','flow','fhigh',...
  'bwcoef','polsamp','rangsamp','mrsmsp','gamma'},definput,varargin);

% Print Settings
if flags.do_print 
  if flags.do_mrs
    fprintf('Settings: PSGE = %1.0f; Gamma = %1.0u; Epsilon = %1.0f deg \n',kv.do,kv.gamma,kv.mrsmsp)
  else
    fprintf('Settings: PSGE = %1.0f; Gamma = %1.0u; Epsilon = 0 deg \n',kv.do,kv.gamma)
  end
end

% HRTF format conversion
if isstruct(target) % Targets given in SOFA format
  kv.fs = target.Data.SamplingRate;
  [target,tang] = extractsp( kv.lat,target );
end

if size(target,3) == 1 && size(target,2) == 2
  target = reshape(target,size(target,1),1,2);
end

if isstruct(template) % Template given in SOFA format
  [template,rang] = extractsp( kv.lat,template );
end

if size(template,3) == 1 && size(template,2) == 2
  template = reshape(template,size(template,1),1,2);
end


%% Stimulus 
if isempty(kv.stim) 
    kv.stim = [1;0];
    kv.fsstim = kv.fs;
elseif isempty(kv.fsstim) 
    kv.fsstim = kv.fs;
end


%% DTF filtering, Eq.(1)
if ~isequal(kv.fs,kv.fsstim)
    amtdisp('Sorry, sampling rate of stimulus and HRIRs must be equal!')
    return
end

tmp = lconv(target,kv.stim);
target = reshape(tmp,[size(tmp,1),size(target,2),size(target,3)]);
    

%% Spectral Analysis, Eq.(2)

if kv.space == 1 % Standard spacing of 1 ERB
  [ireptar,fc] = auditoryfilterbank(target(:,:),kv.fs,...
      'flow',kv.flow,'fhigh',kv.fhigh);
  ireptem = auditoryfilterbank(template(:,:),kv.fs,...
      'flow',kv.flow,'fhigh',kv.fhigh);
else
  fc = audspacebw(kv.flow,kv.fhigh,kv.space,'erb');
  [bgt,agt] = gammatone(fc,kv.fs,'complex');
  ireptar = 2*real(ufilterbankz(bgt,agt,target(:,:)));  % channel (3rd) dimension resolved
  ireptem = 2*real(ufilterbankz(bgt,agt,template(:,:)));
end
Nfc = length(fc);   % # of bands

% Set back the channel dimension
ireptar = reshape(ireptar,[size(target,1),Nfc,size(target,2),size(target,3)]);
ireptem = reshape(ireptem,[size(template,1),Nfc,size(template,2),size(template,3)]);

% Averaging over time (RMS)
ireptar = 20*log10(squeeze(rms(ireptar)));      % in dB
ireptem = 20*log10(squeeze(rms(ireptem)));

if size(ireptar,2) ~= size(target,2) % retreive polar dimension if squeezed out
    ireptar = reshape(ireptar,[size(ireptar,1),size(target,2),size(target,3)]);
end
if size(ireptem,2) ~= size(template,2) % retreive polar dimension if squeezed out
    ireptem = reshape(ireptem,[size(ireptem,1),size(template,2),size(template,3)]);
end

%% Positive spectral gradient extraction, Eq.(3)

nrep.tem = zeros(size(ireptem,1)-kv.do,size(ireptem,2),size(ireptem,3)); %init
nrep.tar = zeros(size(ireptar,1)-kv.do,size(ireptar,2),size(ireptar,3)); %init
for ch = 1:size(ireptar,3)
  if kv.do == 1 % DCN inspired feature extraction
      nrep.tem(:,:,ch) = max(diff(ireptem(:,:,ch),kv.do),0);
      nrep.tar(:,:,ch) = max(diff(ireptar(:,:,ch),kv.do),0);
  elseif kv.do == 2 % proposed by Zakarauskas & Cynader (1993)
      nrep.tem(:,:,ch) = diff(ireptem(:,:,ch),kv.do);
      nrep.tar(:,:,ch) = diff(ireptar(:,:,ch),kv.do);
  else
      nrep.tem(:,:,ch) = ireptem(:,:,ch);
      nrep.tar(:,:,ch) = ireptar(:,:,ch);
  end
end

%% Comparison process, Eq.(4)

sigma=zeros(size(ireptem,2),size(ireptar,2),size(ireptem,3)); % init
for ch = 1:size(ireptar,3)
  for it = 1:size(ireptar,2)
    isd = repmat(nrep.tar(:,it,ch),[1,size(nrep.tem(:,:,ch),2),1]) - nrep.tem(:,:,ch); 
    if kv.do == 0
      sigma(:,it,ch) = sqrt(squeeze(var(isd))); % standard dev. across frequencies
    else
      sigma(:,it,ch) = mean(abs(isd)); % L1-norm across frequencies
    end
  end
end

%% Similarity estimation, Eq.(5)

si=zeros(size(sigma)); % init
for ch = 1:size(ireptar,3)
  for it = 1:size(ireptar,2)
%     si(:,it,ch) = 1+eps - (1+exp(-kv.gamma*(sigma(:,it,ch)-kv.S))).^-1; % baumgartner2014
    si(:,it,ch) = exp(-0.5*(sigma(:,it,ch)/kv.S).^2); % Langendijk
  end
end


%% Binaural weighting, Eq.(6)

if size(si,3) == 2
    binw = 1./(1+exp(-kv.lat/kv.bwcoef)); % weight of left ear signal with 0 <= binw <= 1
    si = binw * si(:,:,1) + (1-binw) * si(:,:,2);
end

dext = si;  

end