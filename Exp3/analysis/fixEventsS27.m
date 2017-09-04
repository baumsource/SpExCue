function EEG = fixEventsS27(EEG)
% EEG = fixEventsS27(EEG)
%
% old binary (logical: 1/0) trigger code: 
%   (5,6) ITD (y/n), ILD (y/n), 
%   (4)   L/R, 
%   (3)   gender(f/m), 
%   (2)   lead/lag, 
%   (1)   1st/2nd syllable
%
% new binary (logical: 1/0) trigger code: 
%   (7)   distractor/target, 
%   (5,6) ITD (y/n), ILD (y/n), 
%   (4)   L/R, 
%   (3)   gender(f/m), 
%   (2,1) 3rd syllable, 2nd syllable
%
% new trigger of cue: 250 (decimal)
%
% response triggers old->new: 
%   100->251 (right), 101->252 (wrong), 102->253 (no response)

type = [EEG.event.type];
latency = [EEG.event.latency];
urevent = [EEG.event.urevent];

%% update response triggers
type(type==100) = 251;
type(type==101) = 252;
type(type==102) = 253;

%% add cue trigger 0.8 s before 1st syllable
dlatCue = round(0.8*EEG.srate);
for ievent = 17:2:63 % all first-syllable events
  idt = type == ievent;
  latcue = latency(idt) - dlatCue;
  latency = [latency,latcue];
  typeCue = 250*ones(1,sum(idt));
  type = [type,typeCue];
end

%% add 2nd-syllable distractor triggers
dlatTD = round(0.2*EEG.srate);
for ievent = 16:2:62 % all 2nd-syllable events
  idt = type == ievent;
  if ievent/4 == round(ievent/4) % target lagging
    latD = latency(idt) - dlatTD; % distractor occurs earlier
  else % target leading
    latD = latency(idt) + dlatTD; % distractor occurs later
  end
  latency = [latency,latD];
  typeDbin = ['1',dec2bin(ievent,6)]; % add distractor bin
  for ib = 4:5 % invert side and gender
    typeDbin(ib) = num2str(1-str2num(typeDbin(ib)));
  end
  typeD = bin2dec(typeDbin) * ones(1,sum(idt));
  type = [type,typeD];
end

%% add 3rd-syllable trigger 0.4 s after 2nd syllable and convert trigger format
typePreConv23 = type;
dlatSyl2to3 = round(0.4*EEG.srate);
for ievent = 16:2:126 % all second-syllable events
  idt = typePreConv23 == ievent;
  lat3 = latency(idt) + dlatSyl2to3;
  latency = [latency,lat3];
  typeBinCode = dec2bin(ievent);
  typeNew2 = bin2dec([typeBinCode(1:end-2),'01']);
  typeNew3 = bin2dec([typeBinCode(1:end-2),'10']);
  type(idt) = typeNew2;
  type3 = typeNew3*ones(1,sum(idt));
  type = [type,type3];
end

%% Convert first-syllable trigger format
for ievent = 17:2:63
  idt = typePreConv23 == ievent;
  typeBinCode = dec2bin(ievent);
  typeNew1 = bin2dec([typeBinCode(1:end-2),'00']);
  type(idt) = typeNew1;
end

%% Sort by latency and r
lt = sortrows([latency;type]');
for ii = 1:size(lt,1)
  EEG.event(ii).latency = lt(ii,1);
  EEG.event(ii).type = lt(ii,2);
  EEG.event(ii).urevent = ii;
end

end