% listen to concatenated stimuli with regularly ordered spectral contrast

ID = 'RB';
M = [1,0,1];%.25,.5,.75,1,1.5];
% M = [1,0,1,0];
pos = [90,0];

stim = SpExCue_stim(M,ID,'pos',pos,'noiseBurst','HRCeq','flow',1e3,'fhigh',16e3);

concatStim = stim.sig{1};
for ii = 2:length(M)
  concatStim = cat(1,concatStim,stim.sig{ii});
end

%% Spectrogram
figSgram = figure;
subplot(1,2,1)
sgram(concatStim(:,1),stim.fs,'dynrange',60,'db')
title('Left')
subplot(1,2,2)
sgram(concatStim(:,2),stim.fs,'dynrange',60,'db')
title('Right')

%% Playback
% pause
sound(concatStim,stim.fs)