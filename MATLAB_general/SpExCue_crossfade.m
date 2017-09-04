function varargout = SpExCue_crossfade(sig1,sig2,fs,dur,tcross,durfade,fadeIOflag)
%SpExCue_crossfade creates cross-faded stimulus pair (cos^2 fade)
%   Usage: [sigpair,nM2] = SpExCue_crossfade(sig1,sig2,fs,dur,tcross,durfade,fadeIOflag)
%
%   Input parameters:
%     sig1    : first stimulus (time in first dimension)
%     sig2    : second stimulus (time in first dimension)
%     fs      : sampling rate of signals
%     dur     : overall duration of stimulus pair in sec
%     tcross  : center time of cross-fade in sec
%     durfade : duration of fades (in and out) in sec
%     fadeIOflag : flag to also fade paired signal in and out
%
%   Output parameters:
%     sigpair : cross-faded stimulus pair
%     nM2     : index of center time of cross-fade

if not(exist('fadeIOflag','var')) 
  fadeIOflag = false;
end

n1 = length(sig1); % length of first input signal
n2 = length(sig2); % length of second input signal
ntotal = round(dur*fs); % total length of output signal
ncross = round(tcross*fs); % index of cross-fade center

nfade = 2*round(durfade*fs/2)-1; % closest odd # of time indices for fade
nstop1 = ncross+(nfade-1)/2; % offset index of first input signal
nstart2 = ncross-(nfade-1)/2; % onset index of second input signal
nsig2 = ntotal-nstart2+1; % length of second stimulus part
fadein = sin(0:pi/2/(nfade-1):pi/2);
fadeout = fliplr(fadein);
  
if n1 <= ncross % short signal -> no cross-fade
  
  fadedsig1 = postpad(sig1,ntotal);
  nsig2 = ntotal-ncross+1;
  sig2 = postpad(sig2,nsig2);
  fadedsig2 = cat(1,zeros(ncross-1,size(sig2,2)),sig2); % prepad
  
else % long signal -> cross-fade
  
  fader1 = [ones(1,nstop1-nfade),fadeout];
  fadedsig1 = sig1(1:nstop1,:).*repmat(fader1(:),1,2);
  fadedsig1 = postpad(fadedsig1,ntotal);
    
  fader2 = [fadein,ones(1,nsig2-nfade)];
  fadedsig2 = sig2(n2-nsig2+1:n2,:).*repmat(fader2(:),1,2);
  fadedsig2 = cat(1,zeros(nstart2-1,2),fadedsig2); % prepad
  
end

sigpair = fadedsig1 + fadedsig2;

if fadeIOflag
  fader = [fadein,ones(1,ntotal-2*nfade),fadeout]';
  sigpair = sigpair.*repmat(fader,[1,2]);
end

varargout{1} = sigpair;
if nargout == 2
  varargout{2} = ncross;
end

end