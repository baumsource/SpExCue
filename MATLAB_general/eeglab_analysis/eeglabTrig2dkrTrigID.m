%  trigID = eeglabTrig2dkrTrigID(trigs)
%
% This function converts the trigger values resulting from the BDF export
% with 'eeglab'.  The function simply  takes only the last 8 bits and
% converts it back to the integer trigger ID specified in my scripts.
%
% INPUT:
%   trigs = [4041, 3939, 3874];
%   trigID = eeglabTrig2dkrTrigID(trigs)
% 
% OUTPUT:
%   trigID = [201; 99; 34];  
%
%  The output is a column vector regardless of the input vector being a
%  column or row of numbers.

function trigID = eeglabTrig2dkrTrigID(trigs)

tmp = dec2bin(trigs);
trigID = bin2dec(tmp(:,(end-7):end));