function ERP = processBDF_RB(bdfFilename,varargin)
% ERP = processBDF_RB(bdfFilename,trg,tEpo,cf)
%   This function processes the continuous multi-channel EEG data in .bdf
%   form, and returns a bandpass filtered, trial epoched version that can
%   be used for additional epoch-based data analysis.
%
%   Users can load a .bdf file and define the epochs around one or multiple
%   event trigger(s) [trg] and a time range [tEpo].  Filtering is
%   implemented using a zero-phase FIR bandpass filter with passbands
%   defined by [cf].  The output of this function returns a MATLAB struct
%   labeled ERP with fields summarizing processing steps.
%
%   All input parameters can be left blank, however, the function will
%   query the user to enter in specific parameters or use default
%   parameters, as defined below.
%
%   IMPORTANT NOTE:
%       You need to have the following functions in your MATLAB path
%           openbdf.m
%           readbdf.m
%           BPF.m
%       Additionally, make sure EEGLab is not in your path.  It will
%       conflict with the openbdf.m and readbdf.m functions.
%
%--------------------------------------------------------------------------
%
% INPUT VARIABLES:
%   bdfFilename : full path and filename of the .bdf file.  If no file is
%                 provided, a GUI will pop-up and instruct the user to
%                 select the .bdf file
%           trg : event trigger(s) as a row vector.  If none provided, user
%                 will be queried at the command window.  Answers must be
%                 entered as a single number or a vector using the
%                 appropriate syntax (e.g. [1 2 5])
%          tEpo : epoch time range [min max] in seconds.  If none provided, 
%                 user will be queried at the command window.Answers must 
%                 be entered as a vector using the appropriate syntax
%                 (e.g. [-0.100 0.500] seconds, peristimulus)
%            cf : bandpass filter frequency cutoffs (Hz)
%                 (default = [1 30] Hz)
%
% OUTPUT VARIABLE:
%   Data structure, ERP, with the following fields:
%       subjID : subject-specific identification code
%          erp : EEG data matrix [n x channels x trials]
%             trialType : full(trialType); % list of trial type by event trigger
%              triggers : list of all within-epoch trigger events
%                diodes : list of all within-epoch diode events
%                     t : epoch time vector (seconds)
%           preStimTime : pre-stimulus time range for epoch
%          postStimTime : post-stimulus time range for epoch
%                    fs : data sampling frequency (Hz)
%           eyeblinkRej : used SSP eyeblink rejection [0=no, 1=yes]
%  eyeblinkRejThreshold : threshold level for eyeblink detection (µV)
%    baselineCorrection : [0=no, 1=yes]
%             filterObj : bandpass filter object
%         filterCutOffs : [lo hi] filter passband cutoff frequencies (Hz)
%        dataDimensions : '[n x channels x trials]';
%
% Revision history:
%   2015-04-01: adjusted epoching scheme SCB
%   2016-04-25: generalized version to apply broadly to many EEG
%               applications (for WRNMMC)
%   2016-05-02: fixed variables removed and handling of default values
%               solved via SOFAarghelper

definput.keyvals.cf = [1 30];
definput.keyvals.tEpo = [];
definput.keyvals.trg = [];
definput.keyvals.mastRef = [];
definput.keyvals.EOGch = [];
definput.keyvals.minEOGthresh = 25;
definput.keyvals.recordTrigger = 254;
definput.keyvals.ERPthresh = 70;
% definput.keyvals.fs = [];
definput.flags.removeEyeblinks = {'removeEyeblinks','keepEyeblinks'};
definput.flags.componentAnalysis = {'ssp','ica'};
definput.flags.baselineCorrect = {'baselineCorrect','noBaselineCorrect'};    % baseline correction of prestimulus average
definput.flags.caching = {'caching','redo','nocaching'};

posdepnames = {'trg','tEpo','cf'};
[flags,kv]  = SOFAarghelper(posdepnames,definput,varargin);

if not(exist('bdfFilename','var'))
  [fname,pname] = uigetfile('*.bdf');
  bdfFilename = [pname,fname];
end
if isempty(kv.trg)
  kv.trg = input('Which trigger(s) do you want to epoch? ');
end
if isempty(kv.tEpo)
  kv.tEpo = input('Enter peristimulus epoch times: ');
end

% Save/cache data to:
[subjDir,subjID,~] = fileparts(bdfFilename);
cachingDir = mfilename('fullpath');
savename = fullfile(cachingDir,[subjID,'.mat']);
if exist(savename,'file') && flags.do_caching
  cache = load(savename); % 'data','kv','flags','triggers','diodes','fs'
  data = cache.data;
  triggers = cache.triggers;
  diodes = cache.diodes;
  fs = cache.fs;
else
  
  %% FROM .bdf FILES
  % Open and read the bdf file
  fprintf(['Loading file: ',bdfFilename,'...\n'])
  h = openbdf_RB(bdfFilename);
  d = readbdf_RB(h,1:h.Head.NRec);

  if (h.Head.NS<=5)
      fprintf('4-channel setup...\n')
      channels = [1:2];
  elseif (h.Head.NS>5 & h.Head.NS<=41)
      fprintf('32-channel setup...\n');
      channels = [1:32];
  elseif (h.Head.NS>41 & h.Head.NS<=73)
      fprintf('64-channel setup...\n');
      channels = [1:64];
  end

  % Set EXTERNAL ELECTRODES indices: mastoid & EOG
  if isempty(kv.mastRef)
    kv.mastRef = length(channels) + [1 2];
  end
  fprintf('Reference channels: %1u, %1u\n',kv.mastRef(1),kv.mastRef(2));
  if isempty(kv.mastRef)
    kv.EOGch = length(channels) + 3:5;
  end
  fprintf('EOG channels: %1u:%1u\n',kv.EOGch(1),kv.EOGch(end));

  fs = h.Head.SampleRate(1);
  rawBDF = d.Record([channels,kv.mastRef,kv.EOGch],:)';

  % Get and organize event triggers
  triggers = mod(d.Event.triggers,256);
  diodes = d.Event.buttons(:,[4 3 2 1]);
  diodes = cleanupDiodes_RB(max(0,cat(1,zeros(1,4),diff(diodes))));

  %% Process

%   % Resample
%   if not(isempty(kv.fs))
%     Q = fs/gcd(fs,kv.fs);
%     P = kv.fs/gcd(fs,kv.fs);
%     rawBDF = resample(rawBDF,P,Q);
%     triggers = sparse(resample(full(triggers),P,Q)); % TRIGGERS NOT MAINTAINED
%     diodes = resample(diodes,P,Q);
%     fs = kv.fs;
%   end
  
  % Bandpass filter continuous raw dataset
  nOrder = 2*fs; % number of FIR taps
  fprintf('Filtering data (%1.0f to %1.0f Hz)...\n',kv.cf(1),kv.cf(2));
  [data,Hd] = bandpassfilter(rawBDF,fs,nOrder,kv.cf(1),kv.cf(2));

  % Re-reference to mastoid
  data = data-repmat(mean(data(:,kv.mastRef),2),1,size(data,2));

  % Determine if there are multiple recording blocks:
  %   You can get discontinuity artifacts that will affect the eyeblink
  %   removal processing
  OnOff = [0;diff(mod(d.Event.triggers,256))];
  trRec = find(OnOff==kv.recordTrigger);
  for k = 2:length(trRec)
      data(trRec(k)-fix(2*fs):trRec(k)+fix(2*fs),:) = 0;
  end

  if strcmp(bdfFilename,'/Users/rbaumgartner/Documents/ARI/ARIcloud/SpExCue/Experiments/Pilots/EEG 2 (cf Getzmann2010)/Results/SpExCue_EEGpilot2Result_RT_RB.bdf')
    data = data(1:1.322e6,:);
  end

  % Zero-out transient response of FIR filter
  data(1:fix(2*fs),:) = 0;

  % Remove eyeblinks using SSP (Uusatilo & Ilmoniemi)
  if flags.do_removeEyeblinks
    fprintf('Removing eyeblinks...\n');

    if flags.do_ssp
      % EOG signals
      if length(kv.EOGch) < 4 % only external below-eye electrode
        EOG = abs(data(:,1)-data(:,kv.EOGch(1)));
      else
        EOG = abs(data(:,kv.EOGch(1))-data(:,kv.EOGch(2)));
      end
      EOG(1:2*fs) = 0; % zero-out transient response of FIR filter

      % determine and apply threshold
      [f,x] = ecdf(EOG);
      cdfThresh = x(find(f>0.95,1));
      if(cdfThresh>kv.minEOGthresh)
          kv.EOGthresh = cdfThresh;
      else
          kv.EOGthresh = kv.minEOGthresh;
      end
      EOGsamples = EOG>kv.EOGthresh;

      % Remove strongest principal component
      blinkData = data(EOGsamples,:);
      blinkCov = cov(blinkData);
      [E,S] = eig(blinkCov);
      [maxEigVec,IdMaxEigVec] = max(diag(S));
      C = E(:,IdMaxEigVec); % eigenvector with max eigenvalue
      P = eye(size(data,2)) - C*C';
      data = data*P';

    elseif flags.do_ica
      error('ICA-based artifact removal not yet implemented.')
    end
  end
  
  % Cache data
  if not(flags.do_nocaching)
    mkdir(cachingDir)
    save(savename,'data','kv','flags','triggers','diodes','fs')
  end
  cache.kv = kv;
  cache.flags = flags;
end

%% Epoching

OnOff = [0;diff(triggers)];

% Find trigger and diode event samples
% di1 = find(diodes(:,1));   
% di2 = find(diodes(:,2));
% di3 = find(diodes(:,3));
% di4 = find(diodes(:,4));

% Find start of epoch trigger(s)
%   Note: you can enter multiple "start of trial" triggers
trCue = []; % initialize variable
for k = 1:length(kv.trg)
    trCue = cat(1,trCue,find(OnOff==kv.trg(k)));
end

trCue = sort(trCue);
trialType = OnOff(trCue);

% Epoch to desired trigger events
tPre = kv.tEpo(1); 
tPost = kv.tEpo(2);
fprintf('Epoching %1.3f sec to +%1.3f sec...\n',kv.tEpo);
preroll = fix(tPre*fs);
postroll = fix(tPost*fs);
    
for k = 1:length(trCue)
    erp(:,:,k) = data(trCue(k)+preroll:trCue(k)+postroll,:);
    epDiodes(:,:,k) = diodes(trCue(k)+preroll:trCue(k)+postroll,:);
    epTrigs(:,k) = triggers(trCue(k)+preroll:trCue(k)+postroll,:);
end

t = (0:size(erp,1)-1)/fs;
t = t+(tPre);

% Correct for baseline offset
if flags.do_baselineCorrect
  fprintf('Removing baseline re: first %1.0f ms...\n',-tPre*1000);
  tBL = find(t<0.000);
  erp = erp-repmat(mean(erp(tBL,:,:)),[size(erp,1),1,1]);
end

% Apply absolute ERP amplitude threshold
if isscalar(kv.ERPthresh)
  fprintf('Removing epochs exceeding absolute ERP threshold of %1.0fe-6 V...\n',kv.ERPthresh);
  idthresh = max(max(abs(erp))) < kv.ERPthresh;
  erp = erp(:,:,idthresh);
  epTrigs = epTrigs(:,idthresh);
  epDiodes = epDiodes(:,:,idthresh);
end

%% Compile ERP data structure
ERP.subjID = subjID;
ERP.erp = erp; % EEG data [n x channels x trials]
ERP.trialType = full(trialType); % list of trial type by event trigger
ERP.triggers = epTrigs; % list of all within-epoch trigger events
ERP.diodes = epDiodes; % list of all within-epoch diode events
ERP.t = t; % epoch time vector (seconds)
ERP.preStimTime = tPre; % pre-stimulus time range for epoch
ERP.postStimTime = tPost; % post-stimulus time range for epoch
ERP.fs = fs; % data sampling frequency (Hz)
ERP.eyeblinkRej = flags.do_removeEyeblinks; % [0=no, 1=yes]
% ERP.eyeblinkRejThreshold = kv.EOGthresh; % threshold level for eyeblink detection (µV)
ERP.baselineCorrection = flags.do_baselineCorrect; % [0=no, 1=yes]
% ERP.filterObj = Hd; % filter object
ERP.filterCutOffs = kv.cf; % filter passband cutoff frequencies (Hz)
ERP.dataDimensions = '[n x channels x trials]';
ERP.flags = flags;
ERP.keyvals = kv;
ERP.cacheSettings.flags = cache.flags;
ERP.cacheSettings.kv = cache.kv;

end

function [DAT,H1]=openbdf_RB(FILENAME)
% EDF=openedf(FILENAME)
% Opens an EDF File (European Data Format for Biosignals) in MATLAB (R)
% <A HREF="http://www.medfac.leidenuniv.nl/neurology/knf/kemp/edf.htm">About EDF</A> 

%	Copyright (C) 1997-1998 by Alois Schloegl
%	a.schloegl@ieee.org
%	Ver 2.20 	18.Aug.1998
%	Ver 2.21 	10.Oct.1998
%	Ver 2.30 	 5.Nov.1998
%
%	For use under Octave define the following function
% function s=upper(s); s=toupper(s); end;

% V2.12    Warning for missing Header information  
% V2.20    EDF.AS.* changed
% V2.30    EDF.T0 made Y2K compatible until Year 2090

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
% Name changed Sept 6,2002  T.S. Lorig

SLASH='/';		% defines Seperator for Subdirectories
BSLASH=setstr(92);

cname=computer;
if cname(1:2)=='PC' SLASH=BSLASH; end;

fid=fopen(FILENAME,'r','ieee-le');          
if fid<0 
	fprintf(2,['Error LOADEDF: File ' FILENAME ' not found\n']);  
	return;
end;

EDF.FILE.FID=fid;
EDF.FILE.OPEN = 1;
EDF.FileName = FILENAME;

PPos=min([max(find(FILENAME=='.')) length(FILENAME)+1]);
SPos=max([0 find((FILENAME=='/') | (FILENAME==BSLASH))]);
EDF.FILE.Ext = FILENAME(PPos+1:length(FILENAME));
EDF.FILE.Name = FILENAME(SPos+1:PPos-1);
if SPos==0
	EDF.FILE.Path = pwd;
else
	EDF.FILE.Path = FILENAME(1:SPos-1);
end;
EDF.FileName = [EDF.FILE.Path SLASH EDF.FILE.Name '.' EDF.FILE.Ext];

H1=setstr(fread(EDF.FILE.FID,256,'char')');     %
EDF.VERSION=H1(1:8);                     % 8 Byte  Versionsnummer 
%if 0 fprintf(2,'LOADEDF: WARNING  Version EDF Format %i',ver); end;
EDF.PID = deblank(H1(9:88));                  % 80 Byte local patient identification
EDF.RID = deblank(H1(89:168));                % 80 Byte local recording identification
%EDF.H.StartDate = H1(169:176);         % 8 Byte		
%EDF.H.StartTime = H1(177:184);         % 8 Byte		
EDF.T0=[str2num(H1(168+[7 8])) str2num(H1(168+[4 5])) str2num(H1(168+[1 2])) str2num(H1(168+[9 10])) str2num(H1(168+[12 13])) str2num(H1(168+[15 16])) ];

% Y2K compatibility until year 2090
if EDF.VERSION(1)=='0'
        if EDF.T0(1) < 91
                EDF.T0(1)=2000+EDF.T0(1);
        else
                EDF.T0(1)=1900+EDF.T0(1);
        end;
else ;
        % in a future version, this is hopefully not needed   
end;

EDF.HeadLen = str2num(H1(185:192));  % 8 Byte  Length of Header
% reserved = H1(193:236);	         % 44 Byte		
EDF.NRec = str2num(H1(237:244));     % 8 Byte  # of data records
EDF.Dur = str2num(H1(245:252));      % 8 Byte  # duration of data record in sec
EDF.NS = str2num(H1(253:256));       % 8 Byte  # of signals

EDF.Label = setstr(fread(EDF.FILE.FID,[16,EDF.NS],'char')');		
EDF.Transducer = setstr(fread(EDF.FILE.FID,[80,EDF.NS],'char')');	
EDF.PhysDim = setstr(fread(EDF.FILE.FID,[8,EDF.NS],'char')');	

EDF.PhysMin= str2num(setstr(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));	
EDF.PhysMax= str2num(setstr(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));	
EDF.DigMin = str2num(setstr(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));	%	
EDF.DigMax = str2num(setstr(fread(EDF.FILE.FID,[8,EDF.NS],'char')'));	%	

% check validity of DigMin and DigMax
if (length(EDF.DigMin) ~= EDF.NS)
        fprintf(2,'Warning OPENEDF: Failing Digital Minimum\n');
        EDF.DigMin = -(2^15)*ones(EDF.NS,1);
end
if (length(EDF.DigMax) ~= EDF.NS)
        fprintf(2,'Warning OPENEDF: Failing Digital Maximum\n');
        EDF.DigMax = (2^15-1)*ones(EDF.NS,1);
end
if (any(EDF.DigMin >= EDF.DigMax))
        fprintf(2,'Warning OPENEDF: Digital Minimum larger than Maximum\n');
end  
% check validity of PhysMin and PhysMax
if (length(EDF.PhysMin) ~= EDF.NS)
        fprintf(2,'Warning OPENEDF: Failing Physical Minimum\n');
        EDF.PhysMin = EDF.DigMin;
end
if (length(EDF.PhysMax) ~= EDF.NS)
        fprintf(2,'Warning OPENEDF: Failing Physical Maximum\n');
        EDF.PhysMax = EDF.DigMax;
end
if (any(EDF.PhysMin >= EDF.PhysMax))
        fprintf(2,'Warning OPENEDF: Physical Minimum larger than Maximum\n');
        EDF.PhysMin = EDF.DigMin;
        EDF.PhysMax = EDF.DigMax;
end  
EDF.PreFilt= setstr(fread(EDF.FILE.FID,[80,EDF.NS],'char')');	%	
tmp = fread(EDF.FILE.FID,[8,EDF.NS],'char')';	%	samples per data record
EDF.SPR = str2num(setstr(tmp));	%	samples per data record

fseek(EDF.FILE.FID,32*EDF.NS,0);

EDF.Cal = (EDF.PhysMax-EDF.PhysMin)./ ...
    (EDF.DigMax-EDF.DigMin);
EDF.Off = EDF.PhysMin - EDF.Cal .* EDF.DigMin;
tmp = find(EDF.Cal < 0);
EDF.Cal(tmp) = ones(size(tmp));
EDF.Off(tmp) = zeros(size(tmp));

EDF.Calib=[EDF.Off';(diag(EDF.Cal))];
%EDF.Calib=sparse(diag([1; EDF.Cal]));
%EDF.Calib(1,2:EDF.NS+1)=EDF.Off';

EDF.SampleRate = EDF.SPR / EDF.Dur;

EDF.FILE.POS = ftell(EDF.FILE.FID);
if EDF.NRec == -1   % unknown record size, determine correct NRec
  fseek(EDF.FILE.FID, 0, 'eof');
  endpos = ftell(EDF.FILE.FID);
  EDF.NRec = floor((endpos - EDF.FILE.POS) / (sum(EDF.SPR) * 2));
  fseek(EDF.FILE.FID, EDF.FILE.POS, 'bof');
  H1(237:244)=sprintf('%-8i',EDF.NRec);      % write number of records
end; 

EDF.Chan_Select=(EDF.SPR==max(EDF.SPR));
for k=1:EDF.NS
	if EDF.Chan_Select(k)
	    EDF.ChanTyp(k)='N';
	else
	    EDF.ChanTyp(k)=' ';
	end;         
	if findstr(upper(EDF.Label(k,:)),'ECG')
	    EDF.ChanTyp(k)='C';
	elseif findstr(upper(EDF.Label(k,:)),'EKG')
	    EDF.ChanTyp(k)='C';
	elseif findstr(upper(EDF.Label(k,:)),'EEG')
	    EDF.ChanTyp(k)='E';
	elseif findstr(upper(EDF.Label(k,:)),'EOG')
	    EDF.ChanTyp(k)='O';
	elseif findstr(upper(EDF.Label(k,:)),'EMG')
	    EDF.ChanTyp(k)='M';
	end;
end;

EDF.AS.spb = sum(EDF.SPR);	% Samples per Block
bi=[0;cumsum(EDF.SPR)]; 

idx=[];idx2=[];
for k=1:EDF.NS, 
	idx2=[idx2, (k-1)*max(EDF.SPR)+(1:EDF.SPR(k))];
end;
maxspr=max(EDF.SPR);
idx3=zeros(EDF.NS*maxspr,1);
for k=1:EDF.NS, idx3(maxspr*(k-1)+(1:maxspr))=bi(k)+ceil((1:maxspr)'/maxspr*EDF.SPR(k));end;

%EDF.AS.bi=bi;
EDF.AS.IDX2=idx2;
%EDF.AS.IDX3=idx3;


DAT.Head=EDF;  
DAT.MX.ReRef=1;

%DAT.MX=feval('loadxcm',EDF);

return;
end

function [DAT,S]=readbdf_RB(DAT,Records,Mode)
% [DAT,signal]=readedf(EDF_Struct,Records)
% Loads selected Records of an EDF File (European Data Format for Biosignals) into MATLAB
% <A HREF="http://www.medfac.leidenuniv.nl/neurology/knf/kemp/edf.htm">About EDF</A> 
%
% Records1	List of Records for Loading
% Mode		0 	default
%		1 	No AutoCalib
%		Mode+2	Concatanated (channels with lower sampling rate if more than 1 record is loaded)

%	Version 2.11
%	03.02.1998
%	Copyright (c) 1997,98 by Alois Schloegl
%	a.schloegl@ieee.org	
                      
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program has been modified from the original version for .EDF files
% The modifications are to the number of bytes read on line 53 (from 2 to
% 3) and to the type of data read - line 54 (from int16 to bit24). Finally the name
% was changed from readedf to readbdf
% T.S. Lorig Sept 6, 2002
% Arash Fazl added the lines to read the event off the last channel 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3; Mode=0; end;
 
EDF=DAT.Head; 
RecLen=max(EDF.SPR);

S=NaN*zeros(RecLen,EDF.NS);
DAT.Record=zeros(length(Records)*RecLen,EDF.NS);
DAT.Valid=uint8(zeros(1,length(Records)*RecLen));
DAT.Idx=Records(:)';
        
for nrec=1:length(Records),
% nrec
	NREC=(DAT.Idx(nrec)-1);
	if NREC<0; fprintf(2,'Warning READEDF: invalid Record Number %i \n',NREC);end;

	fseek(EDF.FILE.FID,(EDF.HeadLen+NREC*EDF.AS.spb*3),'bof');
	[s, count]=fread(EDF.FILE.FID,EDF.AS.spb,'bit24');
try
	S(EDF.AS.IDX2)=s;
catch
  keyboard; 
end
    
	%%%%% Test on  Over- (Under-) Flow
%	V=sum([(S'==EDF.DigMax(:,ones(RecLen,1))) + (S'==EDF.DigMin(:,ones(RecLen,1)))])==0;
	V=sum([(S(:,EDF.Chan_Select)'>=EDF.DigMax(EDF.Chan_Select,ones(RecLen,1))) + ...
	       (S(:,EDF.Chan_Select)'<=EDF.DigMin(EDF.Chan_Select,ones(RecLen,1)))])==0;
	EDF.ERROR.DigMinMax_Warning(find(sum([(S'>EDF.DigMax(:,ones(RecLen,1))) + (S'<EDF.DigMin(:,ones(RecLen,1)))]')>0))=1;
	%	invalid=[invalid; find(V==0)+l*k];
                             
	if floor(Mode/2)==1
		for k=1:EDF.NS,
			DAT.Record(nrec*EDF.SPR(k)+(1-EDF.SPR(k):0),k)=S(1:EDF.SPR(k),k);
		end;
	else
		DAT.Record(nrec*RecLen+(1-RecLen:0),:)=S;
	end;

	DAT.Valid(nrec*RecLen+(1-RecLen:0))=V;
end;
if rem(Mode,2)==0	% Autocalib
	DAT.Record=[ones(RecLen*length(Records),1) DAT.Record]*EDF.Calib;
end;                   

DAT.Record=DAT.Record';


% get the events bits out 
events_ful=DAT.Record(end,:);
events_ful= events_ful-min(events_ful(:)); % take out the baseline
events_ind=find(events_ful);
events=events_ful(events_ind);

MSB_bit=events<0; % if this is negative, the Most Significant Bit of the event was one
ctemp=zeros(numel(events),23);
for n = 1:23, ctemp(:,n) = bitget(events,n);end %get the bit value of all 23 bits
gg=([ ctemp , MSB_bit' ]); % the 24th bit is the MSB
% clear ctemp MSB_bit events
tt=gg(:,8:-1:1); 
tt=num2str(tt);
tt=bin2dec(tt);
triggers=zeros(numel(events_ful),1);
triggers(events_ind)=tt;
DAT.Event.triggers=sparse(triggers); % the decimal trigger code
buttons=zeros(numel(events_ful),16); % grabs all 16 bits (edit: SB 2015-Feb-26)
buttons(events_ind,:)=gg(:,16:-1:1);
DAT.Event.buttons=buttons;
return; 
end

function [y,Hd] = bandpassfilter(x, fs, norder, cf1, cf2)

% y = randn(1,10*fs)*0.1; % make 10 seconds noise
% n = 10000;
% n = 1000;
% n = 100;
n = norder;

Ny = fs/2;
Wn = [cf1 cf2]/(Ny);
Hd = fir1(n,Wn);

y = flipud(fftfilt(Hd,flipud(fftfilt(Hd,x))));
end

function [diodes,tagIt] = cleanupDiodes_RB(diodes)

nEpochs = size(diodes,3);
nTriggers = size(diodes,2);
nSamples = size(diodes,1);

thresh = 1000;

for k = 1:nEpochs
    for m = 1:nTriggers
        
        wTrigs = find(diodes(:,m,k));
        n = find(diff(find(diodes(:,m,k)))<=thresh);
        
        if(sum(n)>0)
            diodes(wTrigs(n+1),m,k) = 0;
            tagIt(k,m) = 1;
        end
            
    end
end

end