function exampleForDooyong
% Sound example

ele = 0; 
azi = 90:-30:-90;%[30,-30];%repmat(1:360,[1,length(ele)]); 
ele = repmat(ele,[360,1]); 
ele=ele(:);

% Obj = SOFAload('RIEC_hrir_subject_001.sofa');
% Obj = SOFAload('DS_hrtf B.sofa');
Obj = SOFAload('RB_3ms.sofa');

sig = noise(10*Obj.Data.SamplingRate,1,'white');
[out,A,E] = SOFAspat(sig,Obj,azi,ele); out = setdbspl(out,80);
sound(out,Obj.Data.SamplingRate)
% plot
figure; 
subplot(121)
sgram(out(:,1),Obj.Data.SamplingRate)
subplot(122)
sgram(out(:,2),Obj.Data.SamplingRate)
% audiowrite('exampleForDooyong.wav',out,Obj.Data.SamplingRate)