% listen to concatenated stimuli with regularly ordered spectral contrast

ID = 'S06';
M = [1,0];%.25,.5,.75,1,1.5];
% M = [1,0,1,0];
pos = [0,0];
flow = 800;
fhigh = 16e3;

stim = SpExCue_stim(M,ID,'pos',pos,'speech','flow',flow,'fhigh',fhigh,'SPL',60,'HRCeq');

% Concatenate longer set of Ms
% concatStim = stim.sig{1};
% for ii = 2:length(M)
%   concatStim = cat(1,concatStim,stim.sig{ii});
% end

% OR: Create a pair as in experiment
concatStim = SpExCue_crossfade(stim.sig{1},stim.sig{2},...
          stim.fs,1.2,.6,.1);
        
% stimDTF = SpExCue_stim(M,ID,'pos',pos,'continuousNoise','flow',flow,'fhigh',fhigh,'SPL',60,'DTF');
% concatStim = SpExCue_crossfade(stimDTF.sig{1},stim.sig{1},...
%           stim.fs,1.2,.6,.1);
        
%% Spectrogram
figSgram = figure;
subplot(1,2,1)
sgram(concatStim(:,1),stim.fs,'dynrange',60,'db')
set(gca,'YLim',[flow-500,fhigh+500])
title('Left')
subplot(1,2,2)
sgram(concatStim(:,2),stim.fs,'dynrange',60,'db')
set(gca,'YLim',[flow-500,fhigh+500])
title('Right')

set(figSgram,'PaperUnits','centimeters','PaperPosition',[100,100,10,6])
print(figSgram,'-r600','-dpng',mfilename)

%% Playback
% pause
% sound(concatStim,stim.fs)
