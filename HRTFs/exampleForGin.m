function exampleForGin
% Sound example for Gin

ele = 40:20:80; 
azi = repmat(1:360,[1,length(ele)]); 
ele = repmat(ele,[360,1]); 
ele=ele(:);

Obj = SOFAload('VB_hrtf B.sofa');

sig = noise(10*Obj.Data.SamplingRate,1,'white');
[out,A,E] = SOFAspat(sig,Obj,azi,ele); out = setdbspl(out,80);
% sound(out,Obj.Data.SamplingRate)
audiowrite('exampleForGin.wav',out,Obj.Data.SamplingRate)