function spatializedReferenceSounds
% Sound example

ID = 'Yuqi';

azi = 90:-30:-90;
ele = zeros(1,length(azi));

fn = fullfile('..','HRTFs',[ID,'_eq.sofa']);
Obj = SOFAload(fn);

sig = noise(10*Obj.Data.SamplingRate,1,'white');
[out,A,E] = SOFAspat(sig,Obj,azi,ele); 
out = out/2/max(out(:));
sound(out,Obj.Data.SamplingRate)

% plot
figure; 
subplot(121)
sgram(out(:,1),Obj.Data.SamplingRate)
subplot(122)
sgram(out(:,2),Obj.Data.SamplingRate)

% audiowrite([mfilename,'.wav'],out,Obj.Data.SamplingRate)