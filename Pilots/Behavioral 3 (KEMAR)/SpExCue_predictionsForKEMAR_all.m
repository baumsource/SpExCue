% Externalization with KEAMR, function of azimuth

HRTFpath = strrep(which('SpExCue_stim'),...
  fullfile('MATLAB_general','SpExCue_stim.m'),'HRTFs');
KEMAR = SOFAload(fullfile(HRTFpath,['KEMAR_hrtf B.sofa']));
idhor = KEMAR.SourcePosition(:,2) == 0;
idfront = KEMAR.SourcePosition(:,1) <= 90 | KEMAR.SourcePosition(:,1) >= 270;
idsel = idhor & idfront;
KEMAR.Data.IR = KEMAR.Data.IR(idsel,:,:);
KEMAR.SourcePosition = KEMAR.SourcePosition(idsel,:);
KEMAR = SOFAupdateDimensions(KEMAR);

azi = KEMAR.SourcePosition(:,1);
lat = mod(azi+90,180)-90;
lat(azi==90) = 90;

data = data_baumgartner2014;

for ii = 1:length(data)

% ID = 'RB';

% templates = SOFAload(fullfile(HRTFpath,[ID '_hrtf B.sofa']));
templates = data(ii).Obj;
idhor = templates.SourcePosition(:,2) == 0;
idfront = templates.SourcePosition(:,1) <= 90 | templates.SourcePosition(:,1) >= 270;
idsel = idhor & idfront;
templates.Data.IR = templates.Data.IR(idsel,:,:);
templates.SourcePosition = templates.SourcePosition(idsel,:);
templates.MeasurementAudioLatency = templates.MeasurementAudioLatency(idsel,:);
templates.MeasurementSourceAudioChannel = templates.MeasurementSourceAudioChannel(idsel);
templates = SOFAupdateDimensions(templates);

%%
for aa = 1:length(lat)
  dext_kemar(aa,ii) = baumgartner2015externalization( KEMAR,templates,'lat',lat(aa) );
%   dext(aa) = max(dext_kemar);
  dext_ref(aa,ii) = baumgartner2015externalization( templates,templates,'lat',lat(aa) );
end
end
%%
dext = dext_kemar./dext_ref;
[lat,idsort] = sort(lat);
dext = dext(idsort,:);
dext = dext.^.05;
dext_med = median(dext,2);
dext_min = min(dext,[],2);
dext_max = max(dext,[],2);

%% Plot
fig = figure;
h(1) = patch([lat;flipud(lat)],[dext_max;flipud(dext_min)],0.5*ones(1,3));
hold on
h(2) = plot(lat,dext_med,'k-');
set(h(2),'LineWidth',2)
legend('Range','Median')
box on
set(gca,'XLim',[-90,90],'YLim',[0.5,1],'XTick',-90:30:90)
% h = errorbar(lat,dext_m,dext_sd,'ko-');
% set(h,'MarkerFaceColor','k')
xlabel('Lateral angle (deg)')
ylabel('Degree of externalization')

%% Save
FontSize = 10;
Resolution = '-r600';
set(findall(fig,'-property','FontSize'),'FontSize',FontSize)

if not(exist(mfilename,'dir'))
  mkdir(mfilename)
end
fn = fullfile(mfilename,mfilename);
set(fig(1),'PaperUnits','centimeters','PaperPosition',[100,100,10,8])
print(fig(1),Resolution,'-dpng',fn)