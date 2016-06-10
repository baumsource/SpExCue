function SpExCue_runPilot(ID)
%SpExCue_runPilot - Experimental routine for pilot experiment of SpExCue
%project testing the effect of spectral magnitude compression of HRTFs on 
%sound externalization 

% AUTHOR: Robert Baumgartner

%% Experimental variables
azi = -90;%[30,0,-30];
ele = 0;%[0,30,0];
M = [1,0.5,0];
Nrep = 3; % repetitions
flow = 700; % lower cut-off frequency
fhigh = 18000; % upper cut-off frequency
dur = 1.5; % duration of stimulus pair in sec
durfade = 0.05; % duration of fade-in/out in sec
jitter = 0.1; % temporal jitter in sec
SPL = 70; % SPL in dB
SPLrove = 10; % roving of SPL in dB
rdepth = 10;
rdensity = 0.25;
flags.do_monaural = false;
flags.SPLrove = 'componentRove'; % 'pairRove' or 'componentRove'

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
  'and try to ignore other differences like intensity or timbre.';...
  'Sometimes it might sound as coming from a specific location outside your head and';...
  'sometimes it might sound as coming from close to your skin or even inside your head.';...
  '';...
  'During the experiment you will be asked to perform the following task:';...
  'If the SECOND sound appears to be CLOSER than the first sound, press the *C* key.';...
  'If the SECOND sound appears to be FARTHER than the first sound, press the *F* key.';...
  'Always respond as fast and accurate as you can while keeping your eyes focused';...
  'on the cross appearing at the center of the screen.';...
  '';...
  'Press any key to hear some exemplary sounds before starting the experiment.';...
  'Use these exemplary sounds to adjust to comfortable loudness.';...
  });

%% Stimulus generation and position permutation
pos = [azi(:),ele(:)];
Npos = length(azi);
pos = pos(randperm(Npos),:);
fs = 48828; % sampling rate of TDT
subj.stim = SpExCue_stim( M,subj.ID,pos,fs,flow,fhigh,SPL );

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
waitforbuttonpress;

for famDur = [2,1]*dur
  
  delete(t)
  if famDur == 2*dur
    t = text(x0,y0,'First a bit slower.');
  else
    t = text(x0,y0,'Now, at testing speed.');
  end
  pause(3)
  
  for jj = 1:Npos

    delete(t)
    if pos(jj,1) > 0
      azilabel = 'left';
    elseif pos(jj,1) < 0
      azilabel = 'right';
    else
      azilabel = 'above';
    end
    t = text(0.5,0.5,['Sounds from front ' azilabel '.']);
    pause(2)
    delete(t)
    t = text(0.5,0.5,'+');
    pause(1)

    sound(SpExCue_crossfade(...
      subj.stim.sig{1,jj},subj.stim.sig{length(M),jj},...
      subj.stim.fs,famDur,famDur/2,0.05),subj.stim.fs);
    pause(famDur+dur)
    sound(SpExCue_crossfade(...
      subj.stim.sig{2,jj},subj.stim.sig{3,jj},...
      subj.stim.fs,famDur,famDur/2,0.05),subj.stim.fs);
    pause(famDur+dur)
    sound(SpExCue_crossfade(...
      subj.stim.sig{3,jj},subj.stim.sig{1,jj},...
      subj.stim.fs,famDur,famDur/2,0.05),subj.stim.fs);
    pause(famDur+dur)

  end
end

delete(t)
t = text(x0,y0,{...
  'These were examples of the different sounds you will hear during the experiment.';...
  'Now you are ready to start the experiment.';...
  '';...
  'Recall the task description:';...
  'If the SECOND sound appears to be CLOSER than the first sound, press the *C* key.';...
  'If the SECOND sound appears to be FARTHER than the first sound, press the *F* key.';...
  'Always respond as fashelp tdtt and accurate as you can while keeping your eyes focused';...
  'on the cross appearing at the center of the screen.';...
  '';...
  'Press any key to start the experiment!';...
  });
waitforbuttonpress;

%% Test procedure

% define all possible combinations of two different M
C = nchoosek(M,2);
Mcomb = [C;fliplr(C)]; % combinations of M
C = nchoosek(1:length(M),2);
iM = [C;fliplr(C)]; % corresponding combinations of indices of M
Ncomb = size(Mcomb,1);
Ntotal = Npos*Nrep*Ncomb;

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
 
for irep = 1:Nrep
  for jj = 1:Npos

    delete(t)
    if pos(jj,1) > 0
      azilabel = 'left';
    elseif pos(jj,1) < 0
      azilabel = 'right';
    else
      azilabel = 'above';
    end
    t = text(0.5,0.5,['Sounds from front ' azilabel '.']);
    pause(2)

    delete(t)
    t = text(0.5,0.5,'+');
    pause(1)
  

    % randomized presentation order
    for iC = randperm(Ncomb);
      ii = ii+1;
      subj.pos(ii,:) = [azi(jj),ele(jj)];
      subj.Mcomb(ii,:) = Mcomb(iC,:);
      sig1 = subj.stim.sig{iM(iC,1),jj};
      sig2 = subj.stim.sig{iM(iC,2),jj};

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
      sigpair = SpExCue_spectralRipple(sigpair,subj.stim.fs,rdepth,rdensity,subj.rphase(ii),'oct');

      tic;
      sound(sigpair,subj.stim.fs)
      waitforbuttonpress; % wait for listener's response
      subj.RT(ii) = toc-dur/2+dt;
      subj.resp(ii) = get(gcf, 'CurrentCharacter'); 

      pause(1.5 + jitter*(rand-0.5))
    end

  end
  
end

delete(t)
t = text(x0,y0,{...
  'Experiment successfully completed. Thank you!';...
  });

%% Save results
fn = fullfile('Results',['SpExCue_pilotResult_' subj.ID]);
if flags.do_monaural
  fn = [fn '_mon'];
end
save(fn,'subj')

pause(2)
close(fig)
end