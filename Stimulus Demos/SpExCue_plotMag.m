function SpExCue_plotMag
% plot Magnitude responses

if not(exist('SpExCue_stim','file'))
  addpath('/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/MATLAB_general')
end

flags.do_print = true;
stimulusFlag = 'IR'; % 'IR' 'speech'
flags.do_erb = true;

azi = 180;
% ele = 0;
M = [1,-1];
ID = 'RB';
% rdepth = 0;
% rdensity = 0.25;
% rphase = 1.5*pi;
flow = 1e2;
fhigh = 18e3;
fs = 44.1e3;

%% Generate flat-source stimulus
[subj.stim,refObj] = SpExCue_stim( M,ID,'azi',azi,'fs',fs,...
  'flow',flow,'fhigh',fhigh,'SPL',70,stimulusFlag,'ARI','noBP');


%% Add spectral ripple 
for rphase = 0%(0:.5:2)*pi;

for iplot = 1%:2

% if iplot == 2
%   for ii = 1:length(M)
%     subj.stim.sig{ii} = SpExCue_spectralRipple(subj.stim.sig{ii},subj.stim.fs,rdepth,rdensity,rphase,'oct');
%   end
% end

Nfft = 2^10;
freq = 0:subj.stim.fs/Nfft:subj.stim.fs/2; % frequency vector
frange = [flow,fhigh];
if flags.do_erb
  freq = freqtoerb(freq);
  frange = freqtoerb(frange);
end
MovAvgL = 300;
for ii = 1:length(M)
  mag{ii} = db(abs(fftreal(subj.stim.sig{ii},Nfft)));
  mag{ii} = mag{ii} + 4.2; % to center at 0dB
  tmp(:,1) = cconv(ones(MovAvgL,1)/MovAvgL,mag{ii}(:,1));
  tmp(:,2) = cconv(ones(MovAvgL,1)/MovAvgL,mag{ii}(:,2));
  mag_smoothed{ii} = tmp(MovAvgL/2:end-MovAvgL/2,:);
  mag_smoothed{ii} = mag_smoothed{ii} - 2.6; % to center at 0dB
end

%% Plot magnitude spectra

fig(iplot) = figure;
for ii = 1:length(M)
  subplot(1,2,1)
%   if strcmp(stimulusFlag,'IR')
    plot(freq,mag{ii}(:,1))
%   else   
% plot(freq,mag_smoothed{ii}(:,1))
%   end
  title('Left')
  hold on
  
  subplot(1,2,2)
%   if strcmp(stimulusFlag,'IR')
    plot(freq,mag{ii}(:,2))
%     yrange = [-10,10];
%   else
%     plot(freq,mag_smoothed{ii}(:,2))
    yrange = [-21,21];
%   end
  title('Right')
  hold on
  legendentry{ii} = ['C = ' num2str(M(ii),'%1.1f')];
end

for pp = 1:2 
  subplot(1,2,pp)
  if flags.do_erb
    XTickLabel = round(logspace(log10(flow*11/10),log10(fhigh*9/10),5));
    XTick = freqtoerb(XTickLabel);
    set(gca,'XTick',XTick,'XTickLabel',XTickLabel)
  end
  xlabel('Frequency (Hz)')
  ylabel('Magnitude (dB)')
  legend(legendentry,'Location','southwest')
%   set(gca,'XTickLabel',round(erbtofreq(get(gca,'XTick'))/10)/100)
  axis([frange,yrange])
end

%% Evaluate and plot azimuth estimation and ILDs
% lookup = itd2anglelookuptable(refObj,refObj.Data.SamplingRate,'dietz2011');
% % fig(iplot+2) = figure;
% fig(2) = figure;
% for ii = 1:length(M)
%   [phi(ii),phi_std,itd,ild,cfreqs] = wierstorf2013estimateazimuth(...
%     subj.stim.sig{ii},lookup,'fs',subj.stim.fs,'dietz2011','rms_weighting');
% %   ild = mean(ild);
% %   plot(cfreqs,ild(length(ild)-length(cfreqs)+1:end)) % sometimes first center frequency omitted
%   [exPattern,fc] = baumgartner2014spectralanalysis(subj.stim.sig{ii},'flow',flow);
%   semilogx(fc,diff(exPattern,1,2))
%   hold on
%   legendentry{ii} = [legendentry{ii},', Azi.: ' num2str(phi(ii),'%2.0f'),'\circ'];
% end
% set(gca,'XLim',[flow,fhigh])
% xlabel('Frequency (Hz)')
% ylabel('ILD (dB)')
% legend(legendentry,'Location','northeast')

end

%% Print
if flags.do_print
  FontSize = 8;
  Resolution = '-r600';

  set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
  set(fig,'PaperUnits','centimeters','PaperPosition',[100,100,14,7])

  fn = mfilename;%fullfile('.',mfilename,mfilename);
  if not(strcmp(ID,'KEMAR'))
    fn = [fn,'_',ID];
  end
  print(fig(1),Resolution,'-dpng',[fn '_flat'])
  % print(fig(2),Resolution,'-dpng',[fn '_spectRove_rphase' num2str(10*rphase/pi)])
%   print(fig(2),Resolution,'-dpng',[fn '_ILD_flat'])
  % print(fig(4),Resolution,'-dpng',[fn '_ILD_spectRove_rphase' num2str(10*rphase/pi)])
  close(fig)
end
end

