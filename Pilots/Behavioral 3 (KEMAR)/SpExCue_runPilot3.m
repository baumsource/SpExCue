function SpExCue_runPilot3(ID)
%SpExCue_runPilot3 - Experimental routine for pilot experiment of SpExCue
%project testing the effect of spectral magnitude compression of HRTFs on 
%sound externalization 
%
% Changes from Pilot to Pilot2: M may remain the same, different stimuli
% Changes from Pilot2 to Pilot3: KEMAR as M=.75 condtion.

% AUTHOR: Robert Baumgartner

%% Experimental variables
azi = 90;%[90,0,-90];
ele = 0;%[0,0,0];
M = [1,.75,.5,.25,0]; % NOTE that 0.75 will be replaced by KEMAR
Nrep = 6; % repetitions
flow = 1e3; % lower cut-off frequency
fhigh = 16e3; % upper cut-off frequency
dur = 1.6; % duration of stimulus pair in sec
durfade = 0.05; % duration of fade-in/out in sec
jitter = 0.1; % temporal jitter in sec
SPL = 70; % SPL in dB
SPLrove = 0; % roving of SPL in dB
rdepth = 0;
rdensity = 0.25;
flags.do_familiarize = true;
flags.do_monaural = false;
flags.do_includeMrepetition = false;
flags.do_analyze = true;
flags.SPLrove = 'componentRove'; % 'pairRove' or 'componentRove'
flags.sourceSignal = {'continuousNoise'}; %'speech','noiseBurst','AMnoiseBurst','continuousNoise'

%% Check availability of dependent functions
if not(exist('SpExCue_stim','file'))
  addpath(fullfile('..','MATLAB_general'))
end
if not(exist('fftreal','file'))
  amtstart
end

%% Black screen and display settings
fig = figure('Units','normalized','Position',[0 0 1 1],'MenuBar','none');
ax = axes('Units','normalized', 'Position',[0 0 1 1]);
axis off;
set(gcf, 'Color', [0 0 0]);

set(gca,'DefaulttextHorizontalAlignment','center')
set(gca,'DefaulttextFontSize',25)
set(gca,'DefaulttextColor','w')

x0 = 0.5;
y0 = 0.5;

%% Enter listener ID
if not(exist('ID','var'))
  ID = inputdlg('Listener ID','Enter ID of listener',1,{''});
  subj.ID = ID{1};
else
  subj.ID = ID;
end

%% Listener instruction
t = text(x0,y0,{...
  'In this experiment you will hear pairs of two different sounds.';... 
  'Only focus on the SPATIAL LOCATION of the sounds';...
  'and try to ignore other differences like timbre.';...
  'Sometimes it might sound as coming from a specific location outside your head and';...
  'sometimes it might sound as coming from close to your skin or even inside your head.';...
  '';...
  'During the experiment you will be asked to perform the following task:';...
  'If the SECOND sound appears to be CLOSER than the first sound, press the *C* key.';...
  'If the SECOND sound appears to be FARTHER than the first sound, press the *F* key.';...
  '';...
  'Press any key to hear some exemplary sounds before starting the experiment.';...
  'Use these exemplary sounds to adjust to comfortable loudness.';...
  });

%% Stimulus generation
Nsource = length(flags.sourceSignal);
pos = [azi(:),ele(:)];
Npos = length(azi);
pos = pos(randperm(Npos),:);
fs = 48828; % sampling rate of TDT
for ii = 1:Nsource
  if ii == 1
    subj.stim = SpExCue_stim( M,subj.ID,pos,fs,flow,fhigh,SPL,flags.sourceSignal{ii});
    kemar = SpExCue_stim( 1,'KEMAR',pos,fs,flow,fhigh,SPL,flags.sourceSignal{ii});
    subj.stim.sig(2,:) = kemar.sig; % replace M=.75 with KEMAR
  else
    tmp = SpExCue_stim( M,subj.ID,pos,fs,flow,fhigh,SPL,flags.sourceSignal{ii});
    kemar = SpExCue_stim( 1,'KEMAR',pos,fs,flow,fhigh,SPL,flags.sourceSignal{ii});
    tmp.sig(2,:) = kemar.sig; % replace M=.75 with KEMAR
    subj.stim.sig = cat(3,subj.stim.sig,tmp.sig);
  end
end

%% Monaural?
if flags.do_monaural
  for ii = 1:numel(subj.stim.sig)
    ILD(ii) = -diff(dbspl(subj.stim.sig{ii}));
    if ILD(ii) > 0
      subj.stim.sig{ii}(:,2) = 0;
    else
      subj.stim.sig{ii}(:,1) = 0;
    end
  end
end

%% Listener familiarization
if flags.do_familiarize
  waitforbuttonpress;

  for famDur = dur %[2,1]*dur

    delete(t)
%     if famDur == 2*dur
%       t = text(x0,y0,'First a bit slower.');
%     else
%       t = text(x0,y0,'Now, at testing speed.');
%     end
    t = text(x0,y0,'Stimulus familiarization.');
    pause(3)

    for jj = 1:Npos

      delete(t)
      if pos(jj,1) > 0
        azilabel = 'left';
      elseif pos(jj,1) < 0
        azilabel = 'right';
      else
        azilabel = 'front';
      end
      t = text(0.5,0.5,['Sounds from ' azilabel '.']);
      pause(2)
      delete(t)
      t = text(0.5,0.5,'+');
      pause(1)

      for isig = 1:Nsource
        % M=1 -> M=0
        sound(SpExCue_crossfade(...
          subj.stim.sig{1,jj,isig},subj.stim.sig{length(M),jj,isig},...
          subj.stim.fs,famDur,famDur/2,0.05),subj.stim.fs);
        pause(famDur+dur)
%         % M=0.25 -> M=0
%         sound(SpExCue_crossfade(...
%           subj.stim.sig{2,jj,isig},subj.stim.sig{3,jj,isig},...
%           subj.stim.fs,famDur,famDur/2,0.05),subj.stim.fs);
%         pause(famDur+dur)
        % M=0 -> M=1
        sound(SpExCue_crossfade(...
          subj.stim.sig{length(M),jj,isig},subj.stim.sig{1,jj,isig},...
          subj.stim.fs,famDur,famDur/2,0.05),subj.stim.fs);
        pause(famDur+dur)
      end

    end
  end
end

delete(t)
t = text(x0,y0,{...
  'These were examples of the different sounds you will hear during the experiment.';...
  'Now you are ready to start the experiment.';...
  '';...
  'Recall the task:';...
  'If the SECOND sound appears to be CLOSER than the first sound, press the *C* key.';...
  'If the SECOND sound appears to be FARTHER than the first sound, press the *F* key.';...
  '';...
  'Press any key to start the experiment!';...
  });
waitforbuttonpress;

%% Test procedure

if flags.do_includeMrepetition
  % Pilot2: define all possible combinations of M
  iM = repmat(1:length(M),[length(M),1]);
  iMtransp = transpose(iM);
  iM = [iM(:),iMtransp(:)];
  Mcomb = M(iM); % combinations of M
  Ncomb = size(Mcomb,1);
else
  % Pilot1: define all possible combinations of two different M
  C = nchoosek(M,2);
  Mcomb = [C;fliplr(C)]; % combinations of M
  C = nchoosek(1:length(M),2);
  iM = [C;fliplr(C)]; % corresponding combinations of indices of M
  Ncomb = size(Mcomb,1);
end

Ntotal = Npos*Nrep*Ncomb*Nsource;

% test all combinations Nrep times with randomized order
ii = 0; % incremental counter
subj.Mcomb = nan(Ntotal,2); % M combination
subj.SPL = nan(Ntotal,1); % SPL
if strcmp(flags.SPLrove,'componentRove')
  subj.SPL = nan(Ntotal,2); % SPL
end
subj.resp = nan(Ntotal,1); % externalization response
subj.RT = nan(Ntotal,1); % reaction time
subj.rphase = nan(Ntotal,1); % spectral ripple phase
subj.source = cell(Ntotal,1); % spectral ripple phase
 
for irep = 1:Nrep
  for jj = 1:Npos

    delete(t)
    if pos(jj,1) > 0
      azilabel = 'left';
    elseif pos(jj,1) < 0
      azilabel = 'right';
    else
      azilabel = 'front';
    end
    t = text(0.5,0.5,['Sounds from ' azilabel '.']);
    pause(2)

    delete(t)
    t = text(0.5,0.5,'+');
    pause(1)
  

    % randomized presentation order
    for cc = randperm(Ncomb*Nsource);
      iC = mod(cc-1,Ncomb)+1; % fast index
      isig = ceil(cc/Ncomb); % slow index
      ii = ii+1;
      subj.pos(ii,:) = [azi(jj),ele(jj)];
      subj.Mcomb(ii,:) = Mcomb(iC,:);
      subj.source(ii) = flags.sourceSignal(isig);
      sig1 = subj.stim.sig{iM(iC,1),jj,isig};
      sig2 = subj.stim.sig{iM(iC,2),jj,isig};

      % level roving of stimulus components
      if strcmp(flags.SPLrove,'componentRove')
        dSPL = SPLrove*(rand(1,2)-0.5);
        subj.SPL(ii,:) = SPL + dSPL;
        sig1 = 10^(dSPL(1)/20)*sig1;
        sig2 = 10^(dSPL(2)/20)*sig2;
      end

      % combine stimulus pairs with temporal jitter of crossfade
      dt = jitter*(rand-0.5);
      sigpair = SpExCue_crossfade(sig1,sig2,...
        subj.stim.fs,dur,dur/2+dt,durfade);

      % level roving of stimulus pair
      if strcmp(flags.SPLrove,'pairRove')
        dSPL = SPLrove*(rand-0.5);
        subj.SPL(ii) = SPL + dSPL;
        sigpair = 10^(dSPL/20)*sigpair;
      end

      % spectral roving
      subj.rphase(ii) = 2*pi*rand;
      sigpair = SpExCue_spectralRipple(sigpair,subj.stim.fs,rdepth,rdensity,subj.rphase(ii),flow,fhigh,'oct');

      tic;
      sound(sigpair,subj.stim.fs)
      waitforbuttonpress; % wait for listener's response
      subj.RT(ii) = toc-dur/2+dt;
      subj.resp(ii) = get(gcf, 'CurrentCharacter'); 

      delete(t); t = text(0.5,0.5,'o');
      pause(0.5)
      delete(t); t = text(0.5,0.5,'+');
      
      pause(1 + jitter*(rand-0.5))
    end

  end
  
  delete(t)
  t = text(x0,y0,{[num2str(irep),' of ',num2str(Nrep),' blocks completed.'];...
    'Press any key to continue!'});
  waitforbuttonpress;
  
end

delete(t)
t = text(x0,y0,{'Experiment successfully completed. Thank you!'});

%% Save results
fnpath = 'Results';
fn = fullfile(fnpath,['SpExCue_pilot3Result_' subj.ID]);
if flags.do_monaural
  fn = [fn '_mon'];
end
% if strcmp(flags.sourceSignal,'speech')
%   fn = [fn '_speech'];
% end
  
if not(exist(fnpath,'dir'))
  mkdir(fnpath)
end
save(fn,'subj')

pause(2)
close(fig)

%% Analyze results
if flags.do_analyze
  SpExCue_analyzePilot3(subj.ID)
end

% end