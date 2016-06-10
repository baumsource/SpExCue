function HRTFmeasurement_HRC_SOFA

%% Settings
ID = 'KEMAR_mls2';

stimFlag = 'MLS'; % options: MLS or expSweep

History = ['Measurements performed in semi-anechoic booth with absorbing ',...
  'wedges on the floor and door. Recording chain: AuPMC002 binaural microphones, ',...
  'MOTU UltraLite (preamps), MOTU 24I/O (ADC). ',...
  'Measurement stimulus: ',stimFlag,...
  'Measurement script: ',mfilename];
recordingDirName = ['HRTF_',date];

speakerAngleList = [-90 -60 -30 0 30 60 90];
% speakerAngleList = [-82.5 -52.5 -22.5 7.5 37.5 67.5 82.5];
% speakerAngleList = [-75 -45 -15 15 45 75];
%speakerAngleList = [-67.5 -37.5 -7.5 22.5 52.5];

samplingRate = 44.1e3; %sample frequency (Hz)

IRonsetThreshold = -20; % dB re maximum

%% Create folders and check dependencies
if not(exist('SOFAarghelper','file'))
  SOFAdir = fullfile('sofa-api-mo-1.0.1','API_MO');
  addpath(SOFAdir);
  SOFAstart
end
if not(exist(recordingDirName,'dir'))
  mkdir(recordingDirName)
end

%% Create stimulus with Nrep repetitions

T = 0.7;
Nrep = 20;

switch stimFlag
  
  case 'MLS'

    mlsOrder = round(log2(T*samplingRate));
    sig = mseq(2,mlsOrder,0,1);
    Nrep = Nrep + 2;

  case 'expSweep' % Exponential sweep

    flow = 50;
    fhigh= 18.5e3; % some headroom for fade-out
    % T = (L-1)/samplingRate;
    t = 0:1/samplingRate:T;
    slewRate = log(fhigh/flow)/T;
    sig = sin(2*pi*flow * 1/slewRate .* (exp(slewRate*t)-1))';
%     fadeOut = cos(0:pi/2/9:pi/2).^2;
%     sig(end-9:end) = sig(end-9:end).*fadeOut(:); % fade out to reduce clicking
    
end
L = length(sig);

% Repeat Nrep times
stim = sig*ones(1,Nrep); %Nrep+10
stim = stim(:);

% Pre- and postpad zeros
Npad = 2^15;
stim = [zeros(Npad,1);stim;zeros(Npad,1)];
stimLength = length(stim);

% Scale down by factor sc to avoid clipping of speakers
sc = 0.1;
stim = sc*stim;

%% ASIO playrec settings
deviceList = playrec('getDevices');
playDeviceID =  deviceList(find(strcmp({deviceList.name},'MOTU PCI ASIO')==1)).deviceID;
playChannelList = [1 2 3 4 5 6 7];
recordChannelList = [1 2];
bufferDuration = 100e-3;
pageSize = round(bufferDuration*samplingRate);
pageBufCount = 3;
runMaxSpeed = true;

%%Playrec initialization check
if playrec('isInitialised')
    if playrec('getSampleRate')~=samplingRate
        fprintf('Changing playrec sample rate from %d to %d\n', playrec('getSampleRate'), samplingRate);
        playrec('reset');
    elseif playrec('getRecDevice')~=playDeviceID
        fprintf('Changing playrec record device from %d to %d\n', playrec('getRecDevice'), playDeviceID);
        playrec('reset');
    elseif playrec('getPlayMaxChannel')<max(playChannelList)
        fprintf('Resetting playrec to configure device to use more output channels\n');
        playrec('reset');
    end
else % initialize
    playrec('init',samplingRate,playDeviceID,playDeviceID,max(playChannelList),max(recordChannelList),pageSize)
end

%% Measurement
for nSpeaker = 1:length(speakerAngleList)
    
    disp(speakerAngleList(nSpeaker))
    inputSignal = zeros(stimLength,length(playChannelList));
    inputSignal(:,nSpeaker) = stim;
    
    %% Start Testing Loop
    
    playrec('delPage'); %clear device memory
    pageNumList = repmat(-1, [1 pageBufCount]);
    firstTimeThrough = true;
    startPoint = 1;
    endPoint = size(inputSignal,1);
    outputSignal = [];
    for startSample = startPoint:pageSize:endPoint
        endSample = min(startSample + pageSize - 1, endPoint);
        
        playSignalTemp = inputSignal(startSample:endSample,:);
        
        pageNumList = [pageNumList playrec('playrec', playSignalTemp, playChannelList,-1,recordChannelList)];
        
        if(firstTimeThrough)
            %This is the first time through so reset the skipped sample count
            playrec('resetSkippedSampleCount');
            firstTimeThrough = false;
        else
            if(playrec('getSkippedSampleCount'))
                fprintf('%d samples skipped!!\n', playrec('getSkippedSampleCount'));
                %return
                %Let the code recover and then reset the count
                firstTimeThrough = true;
            end
        end
        if(runMaxSpeed)
            while(playrec('isFinished', pageNumList(1)) == 0)
            end
        else
            playrec('block', pageNumList(1));
        end
        
        lastRecording = playrec('getRec', pageNumList(1));
        if(~isempty(lastRecording))
            % append data
            outputSignal = [outputSignal; lastRecording];
        end
        
        playrec('delPage', pageNumList(1));
        %pop page number from FIFO
        pageNumList = pageNumList(2:end);
        pause(.0001)
    end
    fn = fullfile(recordingDirName,[ID,'_',stimFlag,'_',num2str(speakerAngleList(nSpeaker)),'deg.wav']);
    wavwrite(outputSignal,samplingRate,24,fn)
    in = double(outputSignal);
    
    % Convolve input with flipped MLS.  This will yield a series of
    % estimates of the impulse response for each channel.
    % Upsample impulse response estimate by a factor of US.  This is necessary
    % for input/output device clock discrepancy correction.
    r_sig_in = fftfilt(flipud(sig),in)/L;
%     r_sig_in = r_sig_in(2*L+[1:Nrep*L],:);
    if strcmp(stimFlag,'MLS') % remove first and last iteration
        r_sig_in = r_sig_in(Npad+2*L-1 + (1:(Nrep-2)*L),:);
        Nrep = Nrep - 2;
    else
        r_sig_in = [r_sig_in;zeros(L,2)];
        r_sig_in = r_sig_in(Npad+L-1 + (1:Nrep*L),:);
    end
    h_est = zeros(L,2);
    for ch = 1:2
        h_est(:,ch) = mean(reshape(r_sig_in(:,ch),L,Nrep),2);
    end
    
    % Look for max abs value channel.  Use this to re-center h_ext
    [m,i_chmax] = max(max(abs(h_est)));
    [m,i_maxval] = max(abs(h_est(:,i_chmax)));
    h_est = circshift(h_est, [round(L/2)-i_maxval 0]);
    
%     plot(h_est(:,1),'Color',rand(1,3)); hold on
    
    % Analyze for main shift
    i_start = zeros(1,2);
    for ch = 1:2,
        [m,i_maxval] = max(abs(h_est(:,ch)));
        i_start(ch) = find(20*log10(abs(h_est(:,ch))) > 20*log10(m)+IRonsetThreshold, 1, 'first');
    end
    i_start = min(i_start);
    
    % Circular shift to 20 before i_start
    h_meas = circshift(h_est, [20-i_start 0]);
    IRlong(nSpeaker,:,:) = h_meas';
    
    % Create windowed version for anechoic
    win = [0.5-0.5*cos(pi*[1:20]'/21); ones(150-40,1); 0.5+0.5*cos(pi*[1:20]'/21)];
    h_meas_win = h_meas(1:length(win),:).*(win*ones(1,2));
    IRshort(nSpeaker,:,:) = h_meas_win';
    
end

%% Create SOFA object
Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

% general information
Obj.Data.SamplingRate = samplingRate;
Obj.GLOBAL_ListenerShortName = ID;
Obj.ListenerPosition = [0 0 0];
Obj.ListenerView = [1 0 0];
Obj.ListenerUp = [0 0 1];
Obj.GLOBAL_History=History;

% positions
azi = -speakerAngleList;
ele = 0; % elevation
r = 1.5; % radius
Obj.SourcePosition = [azi(:),repmat(ele,length(azi),1),repmat(r,length(azi),1)];

% whole IR
Obj.Data.IR = IRlong;
Obj=SOFAupdateDimensions(Obj);
SOFAsave([recordingDirName filesep ID '_' num2str(round(1000*size(IRlong,3)/samplingRate)) 'ms.sofa'],Obj);

% windowed IR
Obj.Data.IR = IRshort;
Obj=SOFAupdateDimensions(Obj);
SOFAsave([recordingDirName filesep ID '_' num2str(round(1000*size(IRshort,3)/samplingRate)) 'ms.sofa'],Obj);

%% Reset
playrec('reset')

%% Plots
% Energy-time plots (TOAs)
figure; 
subplot(121)
SOFAplotHRTF(Obj,'EtcHorizontal',1);
subplot(122)
SOFAplotHRTF(Obj,'EtcHorizontal',2);

% Magnitude responses
figure
noisefloor=-50;
hM=shiftdim(double(squeeze(Obj.Data.IR)),2);
M=(20*log10(abs(fft(hM))));
M=M(1:floor(size(M,1)/2),:,:);  % only positive frequencies
M=M-max(M(:));
M(M<noisefloor)=noisefloor;
f = 0:samplingRate/size(hM,1):(size(M,1)-1)*samplingRate/size(hM,1);
for aa=1:length(azi)
    subplot(1,length(azi),aa)
    plot(f,squeeze(M(:,aa,:)));
    hold on
    aar = length(azi)-aa+1;
    plot(f,squeeze(M(:,aar,2)),'g');
    set(gca,'XScale','log','XLim',[700,18e3])
    xlabel('Frequency (Hz)');
    if aa==1
        ylabel('Magnitude (dB)');
    end
    title(['azimuth: ',num2str(azi(aa))],'Interpreter','none');
end
legend('left','right','R mirr.','Location','southwest')
set(gcf,'Position',[166 770 1707 222])

end

function ms=mseq(baseVal,powerVal,shift,whichSeq)
%		  Maximum length sequence assuming 2,3,5 distinct values
%
%       [ms]=MSEQ(baseVal,powerVal[,shift,whichSeq])
%
%       OUTPUT:
%       ms = generated maximum length sequence, of length basisVal^powerVal-1
%
%       INPUT:
%		  baseVal  -nuber of sequence levels (2,3, or 5 allowed)
%		  powerVal -power, so that sequence length is baseVal^powerVal-1
%		  shift    -cyclical shift of the sequence
%		  whichSeq -sequence istantiation to use 
%		  (numer of sequences varies with powreVal - see the code)

% (c) Giedrius T. Buracas, SNL-B, Salk Institute
% Register values are taken from: WDT Davies, System Identification
% for self-adaptive control. Wiley-Interscience, 1970
% When using mseq code for design of FMRI experiments, please, cite:
% G.T.Buracas & G.M.Boynton (2002) Efficient Design of Event-Related fMRI 
% Experiments Using M-sequences. NeuroImage, 16, 801-813.
 
if nargin<4, whichSeq=1; end
if nargin<3, shift=1; end;

bitNum=baseVal^powerVal-1;

register=ones(powerVal,1);

if baseVal==2,
switch powerVal,
case 2, tap(1).No=[1,2];
case 3, tap(1).No=[1,3];
 		  tap(2).No=[2,3];
case 4, tap(1).No=[1,4];
 		  tap(2).No=[3,4];
case 5, tap(1).No=[2,5];
		  tap(2).No=[3,5];
		  tap(3).No=[1,2,3,5];
		  tap(4).No=[2,3,4,5];
		  tap(5).No=[1,2,4,5];
		  tap(6).No=[1,3,4,5];
case 6, tap(1).No=[1,6];
		  tap(2).No=[5,6];
		  tap(3).No=[1,2,5,6];
		  tap(4).No=[1,4,5,6];
		  tap(5).No=[1,3,4,6];
		  tap(6).No=[2,3,5,6];
case 7, tap(1).No=[1,7];
		  tap(2).No=[6,7];
		  tap(3).No=[3,7];
		  tap(4).No=[4,7];
		  tap(5).No=[1,2,3,7];
		  tap(6).No=[4,5,6,7];
		  tap(7).No=[1,2,5,7];
		  tap(8).No=[2,5,6,7];
		  tap(9).No=[2,3,4,7];
		  tap(10).No=[3,4,5,7];
		  tap(11).No=[1,3,5,7];
		  tap(12).No=[2,4,6,7];
		  tap(13).No=[1,3,6,7];
		  tap(14).No=[1,4,6,7];
		  tap(15).No=[2,3,4,5,6,7];
		  tap(16).No=[1,2,3,4,5,7];
		  tap(17).No=[1,2,4,5,6,7];
		  tap(18).No=[1,2,3,5,6,7];
case 8, tap(1).No=[1,2,7,8];
		  tap(2).No=[1,6,7,8];
		  tap(3).No=[1,3,5,8];
		  tap(4).No=[3,5,7,8];
		  tap(5).No=[2,3,4,8];
		  tap(6).No=[4,5,6,8];
		  tap(7).No=[2,3,5,8];
		  tap(8).No=[3,5,6,8];
		  tap(9).No=[2,3,6,8];
		  tap(10).No=[2,5,6,8];
		  tap(11).No=[2,3,7,8];
		  tap(12).No=[1,5,6,8];
		  tap(13).No=[1,2,3,4,6,8];
		  tap(14).No=[2,4,5,6,7,8];
		  tap(15).No=[1,2,3,6,7,8];
		  tap(16).No=[1,2,5,6,7,8];
case 9, tap(1).No=[4,9];
		  tap(2).No=[5,9];
		  tap(3).No=[3,4,6,9];
		  tap(4).No=[3,5,6,9];
		  tap(5).No=[4,5,8,9];
		  tap(6).No=[1,4,5,9];
		  tap(7).No=[1,4,8,9];
		  tap(8).No=[1,5,8,9];
		  tap(9).No=[2,3,5,9];
		  tap(10).No=[4,6,7,9];
		  tap(11).No=[5,6,8,9];
		  tap(12).No=[1,3,4,9];
		  tap(13).No=[2,7,8,9];
		  tap(14).No=[1,2,7,9];
		  tap(15).No=[2,4,7,9];
		  tap(16).No=[2,5,7,9];
		  tap(17).No=[2,4,8,9];
		  tap(18).No=[1,5,7,9];
		  tap(19).No=[1,2,4,5,6,9];
		  tap(20).No=[3,4,5,7,8,9];
		  tap(21).No=[1,3,4,6,7,9];
		  tap(22).No=[2,3,5,6,8,9];
		  tap(23).No=[3,5,6,7,8,9];
		  tap(24).No=[1,2,3,4,6,9];
		  tap(25).No=[1,5,6,7,8,9];
		  tap(26).No=[1,2,3,4,8,9];
		  tap(27).No=[1,2,3,7,8,9];
		  tap(28).No=[1,2,6,7,8,9];
		  tap(29).No=[1,3,5,6,8,9];
		  tap(30).No=[1,3,4,6,8,9];
		  tap(31).No=[1,2,3,5,6,9];
		  tap(32).No=[3,4,6,7,8,9];
		  tap(33).No=[2,3,6,7,8,9];
		  tap(34).No=[1,2,3,6,7,9];
		  tap(35).No=[1,4,5,6,8,9];
		  tap(36).No=[1,3,4,5,8,9];
		  tap(37).No=[1,3,6,7,8,9];
		  tap(38).No=[1,2,3,6,8,9];
		  tap(39).No=[2,3,4,5,6,9];
		  tap(40).No=[3,4,5,6,7,9];
		  tap(41).No=[2,4,6,7,8,9];
		  tap(42).No=[1,2,3,5,7,9];
		  tap(43).No=[2,3,4,5,7,9];
		  tap(44).No=[2,4,5,6,7,9];
		  tap(45).No=[1,2,4,5,7,9];
		  tap(46).No=[2,4,5,6,7,9];
		  tap(47).No=[1,3,4,5,6,7,8,9];
		  tap(48).No=[1,2,3,4,5,6,8,9];
case 10, tap(1).No=[3,10];
		   tap(2).No=[7,10];
		   tap(3).No=[2,3,8,10];
		   tap(4).No=[2,7,8,10];
		   tap(5).No=[1,3,4,10];
		   tap(6).No=[6,7,9,10];
		   tap(7).No=[1,5,8,10];
		   tap(8).No=[2,5,9,10];
		   tap(9).No=[4,5,8,10];
		   tap(10).No=[2,5,6,10];
		   tap(11).No=[1,4,9,10];
		   tap(12).No=[1,6,9,10];
		   tap(13).No=[3,4,8,10];
		   tap(14).No=[2,6,7,10];
		   tap(15).No=[2,3,5,10];
		   tap(16).No=[5,7,8,10];
		   tap(17).No=[1,2,5,10];
		   tap(18).No=[5,8,9,10];
		   tap(19).No=[2,4,9,10];
		   tap(20).No=[1,6,8,10];
		   tap(21).No=[3,7,9,10];
		   tap(22).No=[1,3,7,10];
		   tap(23).No=[1,2,3,5,6,10];
		   tap(24).No=[4,5,7,8,9,10];
		   tap(25).No=[2,3,6,8,9,10];
		   tap(26).No=[1,2,4,7,8,10];
		   tap(27).No=[1,5,6,8,9,10];
		   tap(28).No=[1,2,4,5,9,10];
		   tap(29).No=[2,5,6,7,8,10];
		   tap(30).No=[2,3,4,5,8,10];
		   tap(31).No=[2,4,6,8,9,10];
		   tap(32).No=[1,2,4,6,8,10];
		   tap(33).No=[1,2,3,7,8,10];
		   tap(34).No=[2,3,7,8,9,10];
		   tap(35).No=[3,4,5,8,9,10];
		   tap(36).No=[1,2,5,6,7,10];
		   tap(37).No=[1,4,6,7,9,10];
		   tap(38).No=[1,3,4,6,9,10];
		   tap(39).No=[1,2,6,8,9,10];
		   tap(40).No=[1,2,4,8,9,10];
		   tap(41).No=[1,4,7,8,9,10];
		   tap(42).No=[1,2,3,6,9,10];
		   tap(43).No=[1,2,6,7,8,10];
		   tap(44).No=[2,3,4,8,9,10];
		   tap(45).No=[1,2,4,6,7,10];
		   tap(46).No=[3,4,6,8,9,10];
		   tap(47).No=[2,4,5,7,9,10];
		   tap(48).No=[1,3,5,6,8,10];
		   tap(49).No=[3,4,5,6,9,10];
		   tap(50).No=[1,4,5,6,7,10];
		   tap(51).No=[1,3,4,5,6,7,8,10];
		   tap(52).No=[2,3,4,5,6,7,9,10];
		   tap(53).No=[3,4,5,6,7,8,9,10];
		   tap(54).No=[1,2,3,4,5,6,7,10];
		   tap(55).No=[1,2,3,4,5,6,9,10];
		   tap(56).No=[1,4,5,6,7,8,9,10];
		   tap(57).No=[2,3,4,5,6,8,9,10];
		   tap(58).No=[1,2,4,5,6,7,8,10];
		   tap(59).No=[1,2,3,4,6,7,9,10];
		   tap(60).No=[1,3,4,6,7,8,9,10];
case 11, tap(1).No=[9,11];
case 12, tap(1).No=[6,8,11,12];
case 13, tap(1).No=[9,10,12,13];
case 14, tap(1).No=[4,8,13,14];
case 15, tap(1).No=[14,15];
case 16, tap(1).No=[4,13,15,16];
case 17, tap(1).No=[14,17];
case 18, tap(1).No=[11,18];
case 19, tap(1).No=[14,17,18,19];
case 20, tap(1).No=[17,20];
case 21, tap(1).No=[19,21];
case 22, tap(1).No=[21,22];
case 23, tap(1).No=[18,23];
case 24, tap(1).No=[17,22,23,24];
case 25, tap(1).No=[22,25];
case 26, tap(1).No=[20,24,25,26];
case 27, tap(1).No=[22,25,26,27];
case 28, tap(1).No=[25,28];
case 29, tap(1).No=[27,29];
case 30, tap(1).No=[7,28,29,30];
otherwise error(sprintf('M-sequence %.0f^%.0f is not defined',baseVal,powerVal))
end;   
elseif baseVal==3,
switch powerVal,
case 2, tap(1).No=[2,1];
		   tap(2).No=[1,1];
case 3, tap(1).No=[0,1,2];
		   tap(2).No=[1,0,2];
		   tap(3).No=[1,2,2];
		   tap(4).No=[2,1,2];
case 4, tap(1).No=[0,0,2,1];
		   tap(2).No=[0,0,1,1];
		   tap(3).No=[2,0,0,1];
		   tap(4).No=[2,2,1,1];
		   tap(5).No=[2,1,1,1];
		   tap(6).No=[1,0,0,1];
		   tap(7).No=[1,2,2,1];
		   tap(8).No=[1,1,2,1];
case 5, tap(1).No=[0,0,0,1,2]; 
		   tap(2).No=[0,0,0,1,2];
		   tap(3).No=[0,0,1,2,2];
		   tap(4).No=[0,2,1,0,2];
		   tap(5).No=[0,2,1,1,2];
		   tap(6).No=[0,1,2,0,2];
		   tap(7).No=[0,1,1,2,2];
		   tap(8).No=[2,0,0,1,2];
		   tap(9).No=[2,0,2,0,2];
		   tap(10).No=[2,0,2,2,2];
		   tap(11).No=[2,2,0,2,2];
		   tap(12).No=[2,2,2,1,2];
		   tap(13).No=[2,2,1,2,2];
		   tap(14).No=[2,1,2,2,2];
		   tap(15).No=[2,1,1,0,2];
		   tap(16).No=[1,0,0,0,2];
		   tap(17).No=[1,0,0,2,2];
		   tap(18).No=[1,0,1,1,2];
		   tap(19).No=[1,2,2,2,2];
		   tap(20).No=[1,1,0,1,2];
		   tap(21).No=[1,1,2,0,2];
case 6, tap(1).No=[0,0,0,0,2,1];
		   tap(2).No=[0,0,0,0,1,1];
		   tap(3).No=[0,0,2,0,2,1];
		   tap(4).No=[0,0,1,0,1,1];
		   tap(5).No=[0,2,0,1,2,1];
		   tap(6).No=[0,2,0,1,1,1];
		   tap(7).No=[0,2,2,0,1,1];
		   tap(8).No=[0,2,2,2,1,1];
		   tap(9).No=[2,1,1,1,0,1];
		   tap(10).No=[1,0,0,0,0,1];
		   tap(11).No=[1,0,2,1,0,1];
		   tap(12).No=[1,0,1,0,0,1];
		   tap(13).No=[1,0,1,2,1,1];
		   tap(14).No=[1,0,1,1,1,1];
		   tap(15).No=[1,2,0,2,2,1];
		   tap(16).No=[1,2,0,1,0,1];
		   tap(17).No=[1,2,2,1,2,1];
		   tap(18).No=[1,2,1,0,1,1];
		   tap(19).No=[1,2,1,2,1,1];
		   tap(20).No=[1,2,1,1,2,1];
		   tap(21).No=[1,1,2,1,0,1];
		   tap(22).No=[1,1,1,0,1,1];
		   tap(23).No=[1,1,1,2,0,1];
		   tap(24).No=[1,1,1,1,1,1];
case 7, tap(1).No=[0,0,0,0,2,1,2];
		   tap(2).No=[0,0,0,0,1,0,2];
		   tap(3).No=[0,0,0,2,0,2,2];
		   tap(4).No=[0,0,0,2,2,2,2];
		   tap(5).No=[0,0,0,2,1,0,2];
		   tap(6).No=[0,0,0,1,1,2,2];
		   tap(7).No=[0,0,0,1,1,1,2];
		   tap(8).No=[0,0,2,2,2,0,2];
		   tap(9).No=[0,0,2,2,1,2,2];
		   tap(10).No=[0,0,2,1,0,0,2];
		   tap(11).No=[0,0,2,1,2,2,2];
		   tap(12).No=[0,0,1,0,2,1,2];
		   tap(13).No=[0,0,1,0,1,1,2];
		   tap(14).No=[0,0,1,1,0,1,2];
		   tap(15).No=[0,0,1,1,2,0,2];
		   tap(16).No=[0,2,0,0,0,2,2];
		   tap(17).No=[0,2,0,0,1,0,2];
		   tap(18).No=[0,2,0,0,1,1,2];
		   tap(19).No=[0,2,0,2,2,0,2];
		   tap(20).No=[0,2,0,2,1,2,2];
		   tap(21).No=[0,2,0,1,1,0,2];
		   tap(22).No=[0,2,2,0,2,0,2];
		   tap(23).No=[0,2,2,0,1,2,2];
		   tap(24).No=[0,2,2,2,2,1,2];
		   tap(25).No=[0,2,2,2,1,0,2];
		   tap(26).No=[0,2,2,1,0,1,2];
		   tap(27).No=[0,2,2,1,2,2,2];
otherwise error(sprintf('M-sequence %.0f^%.0f is not defined',baseVal,powerVal))
end;   
elseif baseVal==5,
switch powerVal,
case 2, tap(1).No=[4,3];
		   tap(2).No=[3,2];
		   tap(3).No=[2,2];
		   tap(4).No=[1,3];
case 3, tap(1).No=[0,2,3];
		   tap(2).No=[4,1,2];
		   tap(3).No=[3,0,2];
		   tap(4).No=[3,4,2];
		   tap(5).No=[3,3,3];
		   tap(6).No=[3,3,2];
		   tap(7).No=[3,1,3];
		   tap(8).No=[2,0,3];
		   tap(9).No=[2,4,3];
		   tap(10).No=[2,3,3];
		   tap(11).No=[2,3,2];
		   tap(12).No=[2,1,2];
		   tap(13).No=[1,0,2];
		   tap(14).No=[1,4,3];
		   tap(15).No=[1,1,3];
case 4, tap(1).No=[0,4,3,3];
		   tap(2).No=[0,4,3,2];
		   tap(3).No=[0,4,2,3];
		   tap(4).No=[0,4,2,2];
		   tap(5).No=[0,1,4,3];
		   tap(6).No=[0,1,4,2];
		   tap(7).No=[0,1,1,3];
		   tap(8).No=[0,1,1,2];
		   tap(9).No=[4,0,4,2];
		   tap(10).No=[4,0,3,2];
		   tap(11).No=[4,0,2,3];
		   tap(12).No=[4,0,1,3];
		   tap(13).No=[4,4,4,2];
		   tap(14).No=[4,3,0,3];
		   tap(15).No=[4,3,4,3];
		   tap(16).No=[4,2,0,2];
		   tap(17).No=[4,2,1,3];
		   tap(18).No=[4,1,1,2];
		   tap(19).No=[3,0,4,2];
		   tap(20).No=[3,0,3,3];
		   tap(21).No=[3,0,2,2];
		   tap(22).No=[3,0,1,3];
		   tap(23).No=[3,4,3,2];
		   tap(24).No=[3,3,0,2];
		   tap(25).No=[3,3,3,3];
		   tap(26).No=[3,2,0,3];
		   tap(27).No=[3,2,2,3];
		   tap(28).No=[3,1,2,2];
		   tap(29).No=[2,0,4,3];
		   tap(30).No=[2,0,3,2];
		   tap(31).No=[2,0,2,3];
		   tap(32).No=[2,0,1,2];
		   tap(33).No=[2,4,2,2];
		   tap(34).No=[2,3,0,2];
		   tap(35).No=[2,3,2,3];
		   tap(36).No=[2,2,0,3];
		   tap(37).No=[2,2,3,3];
		   tap(38).No=[2,1,3,2];
		   tap(39).No=[1,0,4,3];
		   tap(40).No=[1,0,3,3];
		   tap(41).No=[1,0,2,2];
		   tap(42).No=[1,0,1,2];
		   tap(43).No=[1,4,1,2];
		   tap(44).No=[1,3,0,3];
		   tap(45).No=[1,3,1,3];
		   tap(46).No=[1,2,0,2];
		   tap(47).No=[1,2,4,3];
		   tap(48).No=[1,1,4,2];
otherwise error(sprintf('M-sequence %.0f^%.0f is not defined',baseVal,powerVal))
end;
end;

ms=zeros(bitNum,1);
if isempty(whichSeq), whichSeq=ceil(rand(1)*length(tap)); 
else, 
   if whichSeq>length(tap) | whichSeq<1
      disp(sprintf(' wrapping arround!'));
      whichSeq=rem(whichSeq,length(tap))+1;
   end;
end;

weights=zeros(1,powerVal);
if baseVal==2,
	weights(tap(whichSeq).No)=1;
elseif baseVal>2,
   weights=tap(whichSeq).No;
end;
  
%weights
  
for i=1:bitNum
   % calculating next digit with modulo powerVal arithmetic
   %   register, (tap(1).No)

  %ms(i)=rem(sum(register(tap(whichSeq).No)),baseVal);
  ms(i)=rem(weights*register+baseVal,baseVal);
  % updating the register
  
  register=[ms(i);register(1:powerVal-1)];
end

ms=ms(:);
if ~isempty(shift),
   shift=rem(shift, length(ms));
   ms=[ms(shift+1:end); ms(1:shift)];
end;

if baseVal==2,     
  ms=ms*2-1;
elseif baseVal==3, 
  ms(ms==2)=-1;
elseif baseVal==5, 
  ms(ms==4)=-1;
  ms(ms==3)=-2;
else
  error('wrong baseVal!');
end;

end