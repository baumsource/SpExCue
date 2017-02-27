function out = SpExCue_spectralRipple(in,fs,rdepth,rdensity,rphase,flow,fhigh,fflag)
%SpExCue_spectralRipple - apply spectral ripple on stimulus
%   Usage: out = SpExCue_spectralRipple(in,rdepth,rdensity,rphase,flow,fhigh)
%
%   Input parameters:
%     in       : input signal (dimensions: time x channel)
%     fs       : sampling rate
%     rdepth   : ripple depth in dB
%     rdensity : ripples per ERB or octave
%     rphase   : ripple phase in rad
%     flow     : lower frequency bound
%     fhigh    : upper frequency bound
%     fflag    : frequency scaling flag: 'erb' or 'oct'
%
%   Output parameters:
%     out      : spectrally rippled stimulus

% flow = 700;
% fhigh = 18e3;

% Frequency representation
Nfft = size(in,1); % FFT length
tfin = fftreal(in,Nfft); % transfer function
freq = 0:fs/Nfft:fs/2; % frequency vector
idf = freq >= flow & freq <= fhigh; % rippled frequency range
mag = db(abs(tfin)); % magnitude in dB

% Determine ripple shape
if strcmp(fflag,'erb') || strcmp(fflag,'ERB')
  ftrans = freqtoerb(freq(idf)); % ERB scale
else
  ftrans = log2(freq(idf));
end
ripple = rdepth/2 * sin(2*pi*rdensity*ftrans + rphase); % magnitude ripple

% Add ripple
% Option 1: by simply adding the ripple to the magnitude response
ripmag = mag;
ripmag(idf,:,:) = ripmag(idf,:,:) + repmat(ripple(:),[1,size(mag,2),size(mag,3)]); % rippled magnitude 
riptf = 10.^(ripmag/20).*exp(1i*angle(tfin)); % rippled transfer function
out = ifftreal(riptf,Nfft); % rippled time-domain stimulus
% Option 2: by creating a corresponding min-phase filter
% firmag = zeros(size(mag,1),1);
% firmag(idf) = firmag(idf) + ripple(:);
% filttf = 10.^(firmag/20).*exp(1i*angle(tf(:,1)));
% irfilt = ifftreal(filttf,Nfft);
% irfilt = RB_minphase(irfilt,1);
% out = filter(irfilt,1,in);
end