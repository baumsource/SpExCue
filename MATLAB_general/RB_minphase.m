function IRminPhase = RB_minphase(IR,dim)
% RB_minphase - create minimum-phase filter via causal cepstrum
%
% Usage: IRminPhase = RB_minphase(IR,dim)

% RB, 2016/6/3

% Nfft = 2^nextpow2(size(IR,dim));
Nfft = size(IR,dim);

TF = fft(IR,Nfft,dim);
logTF = log(abs(TF)+eps);

cep = ifft(logTF,Nfft,dim);
Nshift = mod(dim-1,ndims(cep));
cep1 = shiftdim(cep,Nshift);
cep1(Nfft/2+2:Nfft,:,:,:,:) = 0;    % set non-causal part to zero and 
cep1(2:Nfft/2,:,:,:,:) = 2*cep1(2:Nfft/2,:,:,:,:);    % multiply causal part by 2 (due to symmetry)
cepMinPhase = shiftdim(cep1,ndims(cep)-Nshift);

logTFminPhase = fft(cepMinPhase,Nfft,dim);
TFminPhase = exp(logTFminPhase);
IRminPhase = ifft(TFminPhase,Nfft,dim);

end