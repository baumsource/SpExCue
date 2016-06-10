% Externalization with KEAMR, function of azimuth

ID = 'RB';

HRTFpath = strrep(which('SpExCue_stim'),...
  fullfile('MATLAB_general','SpExCue_stim.m'),'HRTFs');

templates = SOFAload(fullfile(HRTFpath,[ID '_hrtf B.sofa']));
idhor = templates.SourcePosition(:,2) == 0;
idfront = templates.SourcePosition(:,1) <= 90 | templates.SourcePosition(:,1) >= 270;
idsel = idhor & idfront;
templates.Data.IR = templates.Data.IR(idsel,:,:);
templates.SourcePosition = templates.SourcePosition(idsel,:);
templates.MeasurementAudioLatency = templates.MeasurementAudioLatency(idsel,:);
templates.MeasurementSourceAudioChannel = templates.MeasurementSourceAudioChannel(idsel);
templates = SOFAupdateDimensions(templates);

KEMAR = SOFAload(fullfile(HRTFpath,['KEMAR_hrtf B.sofa']));
idhor = KEMAR.SourcePosition(:,2) == 0;
idfront = KEMAR.SourcePosition(:,1) <= 90 | KEMAR.SourcePosition(:,1) >= 270;
idsel = idhor & idfront;
KEMAR.Data.IR = KEMAR.Data.IR(idsel,:,:);
KEMAR.SourcePosition = KEMAR.SourcePosition(idsel,:);
KEMAR = SOFAupdateDimensions(KEMAR);

%%
azi = KEMAR.SourcePosition(:,1);
lat = mod(azi+90,181)-90;
for aa = 1:length(lat)
  dext_kemar(aa) = baumgartner2015externalization( KEMAR,templates,'lat',lat(aa) );
%   dext(aa) = max(dext_kemar);
  dext_ref(aa) = baumgartner2015externalization( templates,templates,'lat',lat(aa) );
end
dext = dext_kemar./dext_ref;

dext = dext.^.05;

%% Plot
[lat,idsort] = sort(lat);
dext = dext(idsort);
fig = figure;
plot(lat,dext);
xlabel('Lateral angle (deg)')
ylabel('Degree of externalization')

set(gca,'XLim',[-90,90],'YLim',[0.5,1],'XTick',-90:30:90)

%% Save
FontSize = 10;
Resolution = '-r600';
set(findall(fig,'-property','FontSize'),'FontSize',FontSize)

if not(exist(mfilename,'dir'))
  mkdir(mfilename)
end
fn = fullfile(mfilename,[mfilename,'_',ID]);
set(fig(1),'PaperUnits','centimeters','PaperPosition',[100,100,10,10])
print(fig(1),Resolution,'-dpng',fn)