% listen to concatenated stimuli with regularly ordered spectral contrast

ID = 'RB';
M = [0,1];
% M = [1,0,1,0];
azi = 0;
flow = 4e3;
fhigh = 16e3;

% Changes in azimuth
% M = 1;
% pos = [90,0;60,0];

stim = SpExCue_stim(M,ID,'azi',azi,'speech','flow',flow,'fhigh',fhigh,'SPL',60,'HRCeq','noBP');

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
% figSgram = figure;
% subplot(1,2,1)
% sgram(concatStim(:,1),stim.fs,'dynrange',60,'db')
% set(gca,'YLim',[flow-500,fhigh+500])
% title('Left')
% subplot(1,2,2)
% sgram(concatStim(:,2),stim.fs,'dynrange',60,'db')
% set(gca,'YLim',[flow-500,fhigh+500])
% title('Right')
% 
% set(figSgram,'PaperUnits','centimeters','PaperPosition',[100,100,10,6])
% print(figSgram,'-r600','-dpng',mfilename)

%% Playback
% pause
sound(concatStim,stim.fs)
