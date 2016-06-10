function SpExCue_plotMag_varID
% plot Magnitude responses

azi = [90;60;30;0;-30;-60;-90];
ele = [0;0;0;0;0;0;0];
M = 1;
ID = {'EQfrontal1','EQfrontal2','EQfrontal3'};
flow = 700;
fhigh = 18e3;
flags.do_erb = false;

fig(1) = figure;
for ii = 1:length(ID)
  %% Generate flat-source stimulus
  subj.stim = SpExCue_stim( M,ID{ii},[azi,ele],2.1*fhigh,flow,fhigh,80,'HRC3ms');

  %% Magnitude spectrum
  Nfft = 2^14;
  freq = 0:subj.stim.fs/Nfft:subj.stim.fs/2; % frequency vector
  erb = freqtoerb(freq);
  MovAvgL = 300;
  for ii = 1:length(azi)
    mag{ii} = db(abs(fftreal(subj.stim.sig{ii},Nfft)));
    tmp(:,1) = cconv(ones(MovAvgL,1)/MovAvgL,mag{ii}(:,1));
    tmp(:,2) = cconv(ones(MovAvgL,1)/MovAvgL,mag{ii}(:,2));
    mag_smoothed{ii} = tmp(MovAvgL/2:end-MovAvgL/2,:);
  end

  %% Plot magnitude spectra
  for aa = 1:length(azi)
    subplot(length(azi),2,2*(aa-1)+1)
    if flags.do_erb
      plot(erb,mag_smoothed{aa}(:,1))
    else
      plot(freq,mag_smoothed{aa}(:,1))
    end
    if aa==1; title('Left'); end
    hold on
    text(5,40,[num2str(azi(aa)),'\circ'],'FontWeight','bold')

    subplot(length(azi),2,2*(aa-1)+2)
    if flags.do_erb
      plot(erb,mag_smoothed{aa}(:,2))
    else
      plot(freq,mag_smoothed{aa}(:,2))
    end
    if aa==1; title('Right'); end
    hold on
  end

end

legend(ID,'Location','southwest')
for pp = 1:length(azi)*2 
  subplot(length(azi),2,pp)
  if flags.do_erb
    xlabel('Frequency (ERB)')
    axis([freqtoerb([flow,fhigh]),-10,40])
  else
    xlabel('Frequency (Hz)')
    axis([[flow,fhigh],-10,40])
    set(gca,'XScale','log')
  end
  ylabel('Magnitude (dB)')
end

%% Evaluate and plot azimuth estimation and ILDs
% lookup = data_wierstorf2013('itd2anglelookuptable');
% fig(2) = figure;
% for ii = 1:length(M)
%   [phi(ii),phi_std,itd,ild,cfreqs] = wierstorf2013estimateazimuth(subj.stim.sig{ii},lookup,'fs',subj.stim.fs,'dietz2011');
%   ild = mean(ild);
%   plot(cfreqs,ild(length(ild)-length(cfreqs)+1:end)) % sometimes first center frequency omitted
%   hold on
%   legendentry{ii} = [legendentry{ii},', Azi.: ' num2str(phi(ii),'%2.0f'),'\circ'];
% end
% xlabel('Frequency (Hz)')
% ylabel('ILD (dB)')
% legend(legendentry,'Location','northeast')

% end

%% Print
% FontSize = 8;
% Resolution = '-r600';
% 
% set(findall(fig,'-property','FontSize'),'FontSize',FontSize)
% set(fig,'PaperUnits','centimeters','PaperPosition',[100,100,14,20])
% 
% fn = fullfile('.',mfilename);
% print(fig(1),Resolution,'-dpng',[fn '_mag'])
% print(fig(2),Resolution,'-dpng',[fn '_ILD'])
% close(fig)
end

