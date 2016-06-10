% equalize HRTFs by measurements without listener at the center of the
% array (virtual center of listener's head)

fn = dir('*_3ms.sofa');
fnSF = 'EQfrontal1_3ms.sofa';
ObjSF = SOFAload(fnSF);

%%
for ii = 1:length(fn)
  ObjRaw = SOFAload(fn(ii).name);
  ObjEQ = ObjRaw;
  Nfft = 2*size(ObjRaw.Data.IR,3);
  SF_IRminPhase = RB_minphase(ObjSF.Data.IR,3);
  EQfiltTF = 1./fft(SF_IRminPhase,Nfft,3);
%   figure; plot(db(abs(shiftdim(EQfiltTF(:,1,:),2))))
  ObjEQ.Data.IR = ifft(fft(ObjRaw.Data.IR,Nfft,3).*EQfiltTF,Nfft,3);
  ObjEQ.Data.IR = ObjEQ.Data.IR./max(ObjEQ.Data.IR(:));
  ObjEQ.GLOBAL_History = [ObjEQ.GLOBAL_History,' Equalized by: ',fnSF,'.'];
  fnEQ = strrep(fn(ii).name,'3ms','eq');
  SOFAsave(fnEQ,ObjEQ)
end