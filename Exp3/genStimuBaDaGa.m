function stimu = genStimuBaDaGa(fs,azi,Nrep,subID,hrtfPath,lvl )
% genStimuBaDaGa - two differen-sex talkers at +/-30 degrees
% onsets of leading stream: 0s, 0.4s, 1.0s
% onsets of lagging stream: 0s, 0.6s, 1.2s
%
% binary (logical: 1/0) trigger code: 
%   (1)   distractor/target, 
%   (2-3) ITD (y/n), ILD (y/n), 
%   (4)   L/R, 
%   (5)   gender(f/m), 
%   (6-7) 3rd syllable, 2nd syllable
%
% trigger of cue: 250 (decimal)
%
% Requirements:
% 1) SOFA API
% 2) RB_minphase (contained in MATLAB_general)

offset2Cue = 0.2;
offsetCue2Target = 0.8; % time offset between cue and target in seconds
onsets = [0.6,0.4]; % onsets of 2nd syllable for lagging and leading stream, respectively

spatLbl = {'ITD','ILD','HRTF'};
trgSpat = {'10','01','11'}; 
trgValCue = 250;
Nspat = length(spatLbl);

% outputRMS = 0.1; % RMS of output (important to equalize level between different spatialization methods)
lvlRove = 0; % level roving in dB

Ntalkers = 3;
samPath = 'SpeechSamples';
talkerPath = {'BaDaGaHannah','BaDaGaLouis','BaDaGaScott'};
gender = 'fm';
cvFn = {'ba','da','ga'};
Ncv = length(cvFn);
cv = cell(length(cvFn),length(talkerPath));
for itP = 1:Ntalkers
  for icv = 1:Ncv
    [sig,ofs] = audioread(fullfile(samPath,talkerPath{itP},[cvFn{icv},'.wav'])); 
    cv{icv,itP} = sig/sqrt(mean(sig.^2)); % normalize intensities 
%     plot((cv{icv,itP}).^2)
  end
end

cue = cell(Nspat,length(azi));
 
%% HRTF-based

% Load HRTF
if not(exist('SOFAload','file'))
  if IsWin
    addpath('C:\Experiments\Robert\sofa\API_MO')
    SOFAstart
  else
    amtstart
  end
end
if not(exist('RB_minphase','file'))
  if IsWin
    addpath('C:\Experiments\Robert\Experiments_GIT\MATLAB_general')
  else
    addpath('../MATLAB_general')
  end
end
if ~exist('subID','var') || isempty(subID)
  subID = 'KEAMR';
end
if ~exist('hrtfPath','var') || isempty(hrtfPath)
  hrtfPath = fullfile('..','HRTFs');
end
fnHRTF = fullfile(hrtfPath,[subID,'_eq.sofa']);
HRTF = SOFAload(fnHRTF);
if HRTF.Data.SamplingRate ~= ofs
    for cc=1:numel(cv)
        cv{cc} = resample(cv{cc},HRTF.Data.SamplingRate,ofs);
    end
    ofs = HRTF.Data.SamplingRate;
end

% Remove ITD by minimizing the phase
ILD = HRTF;
ILD.Data.IR = RB_minphase(ILD.Data.IR,3);


% Determine ITD
ITD = HRTF;
ITD.Data.IR = zeros(size(HRTF.Data.IR));
corrcoeff=zeros(HRTF.API.M,HRTF.API.R);
for ii=1:ILD.API.M
    for jj=1:ILD.API.R
        [c,lag]=xcorr(squeeze(HRTF.Data.IR(ii,jj,:)),squeeze(ILD.Data.IR(ii,jj,:)),HRTF.API.N*4-1,'none');
        [corrcoeff(ii,jj),idx]=max(abs(c));
        corrcoeff(ii,jj)=corrcoeff(ii,jj)/sum(HRTF.Data.IR(ii,jj,:).^2);
        toaEst(ii,jj)=lag(idx);
        ITD.Data.IR(ii,jj,lag(idx)) = 1;
    end
end
itd = -diff(toaEst,1,2)/fs*1e6;

% [~,toa] = ziegelwanger2014(HRTF);
% itd = -diff(toa.toa,1,2)/fs*1e6;


idR(1) = find(HRTF.SourcePosition(:,1) == azi(1));
idR(2) = find(HRTF.SourcePosition(:,1) == azi(2));
itd = itd(idR);

% Spatialize
SPAT = {ITD,ILD,HRTF};
nOnsets = round(onsets*fs);
h = cell(Ncv,length(azi),Ntalkers,2,length(SPAT));
for ispat = 1:length(SPAT)
  for iazi=1:length(azi)
      for itone=1:2
          for cc=1:numel(cv)
              tmp=SOFAspat(cv{cc},SPAT{ispat},azi(iazi),0);
              tmp=resample(tmp,fs,ofs);
              tmp=[tmp;zeros(nOnsets(itone)-length(tmp),2)];
              h{mod(cc-1,3)+1,iazi,ceil(cc/3),itone,ispat}=tmp; % cv, azi, sex, lead/lag
          end
      end
  end
  cue(ispat,:) = h(1,:,3,2,ispat);
end


%% ITD-based



% % itd = azi/90*700; % microsec
% % itd = 3*0.07/343*sin(pi*azi/180)*1e6;
% for i=1:length(itd)
%     for itone=1:2
%         for cc=1:numel(cv)
%             x{mod(cc-1,3)+1,i,ceil(cc/3),itone} = genbdg(fs,ofs,itd(i),tonelen(itone),cv{cc},cv{1}); % resampling and ITD
%         end
%     end
% end
% cue(1,:) = x(1,:,3,2);
% 
% xh = {x,h(:,:,:,:,1),h(:,:,:,:,2)};

xh = {h(:,:,:,:,1),h(:,:,:,:,2),h(:,:,:,:,3)};

%% Sequences
nOffset2Cue = round(fs*offset2Cue);
nOffsetCue2Target = round(fs*offsetCue2Target);
trialnumber=0;
for iazi=1:length(azi) % azi
    for tarf=1:2 % sex
        for tarfast=1:2 % lead/lag
            for ii=1:Nrep
                toneseqtar=randi(3,1,3);
                toneseqdis=randi(3,1,3);
                for ispat=1:Nspat
                    x = xh{ispat};
                    if tarfast==1
                        targetsound=[x{toneseqtar(1),iazi,tarf,tarfast};x{toneseqtar(2),iazi,tarf,3-tarfast};x{toneseqtar(3),iazi,tarf,3-tarfast};];
                        distractsound=[x{toneseqdis(1),3-iazi,3-tarf,3-tarfast};x{toneseqdis(2),3-iazi,3-tarf,3-tarfast};x{toneseqdis(3),3-iazi,3-tarf,3-tarfast};];
                        leading='distract';% leading means faster smaller interval 0.4 0.4
                        temp=zeros(size(targetsound));
                        temp(1:size(distractsound,1),:)=distractsound;
                        distractsound=temp;
                    else
                        targetsound=[x{toneseqtar(1),iazi,tarf,tarfast};x{toneseqtar(2),iazi,tarf,tarfast};x{toneseqtar(3),iazi,tarf,tarfast};];
                        distractsound=[x{toneseqdis(1),3-iazi,3-tarf,3-tarfast};x{toneseqdis(2),3-iazi,3-tarf,tarfast};x{toneseqdis(3),3-iazi,3-tarf,tarfast};];
                        leading='target';
                        temp=zeros(size(distractsound));
                        temp(1:size(targetsound,1),:)=targetsound;
                        targetsound=temp;
                    end
                    frame=cat(1,zeros(nOffset2Cue,2),...
                      cue{ispat,iazi},...
                      zeros(nOffsetCue2Target-length(cue{ispat,iazi}),2),...
                      distractsound+targetsound);
                    deltaLvl=lvlRove*(rand-0.5);
                    frame = db2mag(lvl+deltaLvl-100)*frame/mean(sqrt(mean(frame.^2))); % equalize overall intensity
                    trialnumber=trialnumber+1;
                    
                    trg0 = trgValCue;
                    % target sound triggers
                    trgT = ['0',trgSpat{ispat},num2str(iazi-1),num2str(2-tarf)];
                    trgT1 = bin2dec([trgT,'00']); % 1st syllable (simultaneous onset)
                    trgT2 = bin2dec([trgT,'01']); % 2nd target syllable (1st staggered)
                    trgT3 = bin2dec([trgT,'10']); % 3rd target syllable (2nd staggered)
                    % distractor sound triggers
                    trgD = ['1',trgSpat{ispat},num2str(2-iazi),num2str(tarf-1)];
                    trgD2 = bin2dec([trgD,'01']); % 2nd distractor syllable (1st staggered)
                    trgD3 = bin2dec([trgD,'10']); % 3rd distractor syllable (2nd staggered)
                    triggers = sortrows(fliplr([...
                      trg0, 1+nOffset2Cue;...
                      trgT1,1+nOffset2Cue+nOffsetCue2Target;...
                      trgT2,1+nOffset2Cue+nOffsetCue2Target+nOnsets(tarfast);...
                      trgT3,1+nOffset2Cue+nOffsetCue2Target+nOnsets(tarfast)+nOnsets(2);...
                      trgD2,1+nOffset2Cue+nOffsetCue2Target+nOnsets(3-tarfast);...
                      trgD3,1+nOffset2Cue+nOffsetCue2Target+nOnsets(3-tarfast)+nOnsets(2)]));
                    
                    stimu(trialnumber)=struct('sound',frame,'relativeLevel',deltaLvl,...
                      'targetlocation',azi(iazi),'targetgender',gender(tarf),...
                      'seqtar',toneseqtar,'leading',leading,'distractseq',toneseqdis,...
                      'spatialization',spatLbl{ispat},'fs',fs,'triggers',triggers);
                      
                end
            end
        end
    end
end

% for ii=1:numel(cue)
%     cue{ii} = outputRMS*cue{ii}/mean(sqrt(mean(cue{ii}.^2)));
% end

end


