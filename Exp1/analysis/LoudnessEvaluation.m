% Loudness evaluation of SpExCue Exp. 1 stimuli

if ~exist('Loudness_ANSI_S34_2007','file')
  addpath('/Users/rbaumgartner/Documents/MATLAB/LoudnessToolbox 1.2')
  amtstart
  warning off
end

tmp = load('SpExCue_Exp1eeg_subjects_tab.mat');
s = tmp.subject;

M = [0,0.5,1];
HRTFs = 'HRCeq';
fs = 48.8e3;

loudness = nan(size(s,1),length(M));
loudnessLevel = loudness;
for ss = 1:size(s,1)
  
  switch s.direction{ss}
    case 'left'
      azi = 90;
      ch = 1;
    case 'front'
      azi = 0;
      ch = 1;
    case 'right'
      azi = -90;
      ch = 2;
    otherwise
      error('RB: direction not assigned!')
  end
  
  stim = SpExCue_stim( M,s.name{ss},'azi',azi,'fs',fs,'HRCeq');
  
  for m = 1:length(M)
    Sig = stim.sig{m};%20e-6*10.^(100*stim.sig{m}/20);
    [loudness(ss,m,1), ~, ~, loudnessLevel(ss,m,1)] = ...
      Loudness_ISO532B_from_sound(Sig(:,ch),fs,0);
    [loudness(ss,m,2), ~, ~, loudnessLevel(ss,m,2)] = ...
    	Loudness_ANSI_S34_2007(Sig(:,ch), fs, 0, 'head','off');
  end
  disp(ss)
end

dLoudness = loudness(:,1:2,:) - repmat(loudness(:,3,:),[1,2,1]);
dLoudnessLevel = loudnessLevel(:,1:2,:) - repmat(loudnessLevel(:,3,:),[1,2,1]);

%% Plot
figure; 
c.M0 = 0.7*[0,0,1];
c.Mi = 0.7*[0,1,0];
% subplot(2,1,1)
% b = boxplot(mtx.err(:,id),{ftstr(id),cohcstr(id),SPLstr(id)},...
%           'plotstyle','compact','medianstyle','line','colors',colors{ll},...
%           'factorgap',[],'labelverbosity','all','symbol','');
boxplot(dLoudnessLevel(:,:),{{'M1';'M1';'M2';'M2'},[0;0.5;0;.5]},...
  'Factorgap',10,'FactorSeparator',1,'Whisker',Inf,...
  'Colors',[c.Mi;c.M0;c.Mi;c.M0])
set(gca,'YLim',[-3.9,0])
hold on
% text(0.5,-0.5,'M1')
% text(3,-0.5,'M2')
ylabel('Loudness level difference (phon)')
% title('acc. to ANSI S3.4-2007 eq. to Moore et al. (JAES 1997)')
% subplot(2,1,2)
% boxplot(L,M)
xlabel('Spectral contrast')
%%
% title('acc. to ISO 532B / DIN 45631 eq. to Zwicker and Fastl (JASJ 1991)')
fn = mfilename;
saveas(gcf,fn)
print(gcf,'-dpng','-r300',fn)