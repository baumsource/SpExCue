function info = SpExCue_Exp3_trigCode(trig)
% info = SpExCue_Exp3_trigCode(trig)
%
% binary (logical: 1/0) trigger code: 
%   (1)   distractor/target, 
%   (2-3) ITD (y/n), ILD (y/n), 
%   (4)   L/R, 
%   (5)   gender(f/m), 
%   (6-7) 3rd syllable, 2nd syllable
%
% trigger for cue: 250 (decimal)
% trigger for correct answer: 251 (decimal)
% trigger for wrong answer: 252 (decimal)
% trigger for no answer: 253 (decimal)

if trig == 250
  info = 'cue';
  disp(info)
  return
elseif trig == 251
  info = 'correct';
  disp(info)
  return
elseif trig == 252
  info = 'wrong';
  disp(info)
  return
elseif trig == 253
  info = 'no response';
  disp(info)
  return
elseif trig == 0
  info = 'N/A';
  disp(info)
  return
end

b = dec2bin(trig);
b = [repmat('0',1,7-length(b)),b]; % fill up with zeros

switch b(1)
  case '1'
    stream = 'distractor: ';
  otherwise
    stream = 'target:     ';
end
  
switch b(2:3)
  case '10'
    spat = 'ITD, ';
  case '01'
    spat = 'ILD,';
  otherwise 
    spat = 'HRTF,';
end

switch b(4)
  case '1'
    side = 'left, ';
  otherwise
    side = 'right,';
end

switch b(5)
  case '1'
    gender = 'female,';
  otherwise
    gender = 'male,  ';
end

switch b(6:7)
  case '10'
    syllable = '3rd syllable.';
  case '01'
    syllable = '2nd syllable.';
  otherwise
    syllable = '1st syllable.';
end

info = [stream,spat,side,gender,syllable];
disp(info)
end