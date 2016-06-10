% The only difference between this preprocessing script and the one without
% 'ICA' in the filename is that ICA is that ICA is also computed in this
% script. 

% [paramsOut] = SpExCue_EEGLABpreprocessing(subjID, expDate)
% [paramsOut] = SpExCue_EEGLABpreprocessing(subjID, expDate, FS)
% [paramsOut] = SpExCue_EEGLABpreprocessing(subjID, expDate, FS, [locutoff, hicutoff])
% [paramsOut] = SpExCue_EEGLABpreprocessing(subjID, expDate, FS, [locutoff, hicutoff], transitionBW, maxRipple )
%
%  This function runs some basic preprocessing on the raw BDF file.
%  Namely, it: 
%   (1) re-references the data to the two external channels [33,34]
%   (2) removes the extra external channels [36:40], recall ch35 is ECOG
%   (3) defines the channel locations per an .locs file
%   (4) re-samples the data at defined sampling rate (DEFAULT: 512)
%   (5) filters the data given default or user input parameters
%
%  The eeglab data files are saved after Steps 2, 4, and 5.  Upone
%  completion of the preprocessing, this function also saves the parameters
%  of the preprocessing, i.e. 'paramsOut' as a .mat file with a timestamped
%  filename.
%
%  INPUTS:
%   subjID          STRING for the subject ID, e.g. '013'
%   expDate         STRING for the date of the experiment in the filename
%                   of the BDF file, e.g. '01112015'.
%
% Optional key value pairs:
%   FS = 100;
%   cutoffs = [0.5,25]; % filter cutoff frequencies in Hz.
%   transitionBandwidth = 1;  % In Hz
%   maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
%   experimentString = 'EEGpilot';  % String to label experiment number for the saved eeglab dataset files
%   dirBDF = '/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Experiments/Pilot/Results/';  % directory with the .bdf/.set files
%   dirLOCS = '/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Tools/EEG/';       % directory with the eeglab locations file
%   referenceChannels = [33,34]; % Indicies of reference channels in BDF file
%   nChansLocsFilename = 'biosemi_eloc.locs';  % filename of the eeglab locations file
%
% Optional flags:
%   ICA = {'','ICA'}; % run ICA or not
%
%  OUTPUT:
%   paramsOut       A structure with all of the important preprocessing
%                   settings.  NOTE: this structure is automatically saved
%                   upon completion of the function. The filename of the
%                   saved structure is timestamped to help identify the
%                   resulting eeglab data files with the .mat strucutre.
%
%
%  Created: Darrin K. Reed (Jan. 14, 2016)
%  Modified from 'dkrEEGLAB_exp2_preprocessICA.m': 
%       Darrin K. Reed (Jan 31, 2016)... added the ICA analysis
% Modified by Robert Baumgartner for SpExCue (Mar 8, 2016)
                
%--------------------------------------------------------------------------

function [paramsOut] = SpExCue_EEGLABpreprocessing_Cluster(subjID, expDate, varargin)

definput.keyvals.FS = 100;
definput.keyvals.cutoff = [0.5,25];
definput.keyvals.transitionBandwidth = 1;  % In Hz
definput.keyvals.maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
definput.keyvals.experimentString = 'EEGpilot';  % String to label experiment number for the saved eeglab dataset files
definput.keyvals.dirBDF = '.';  % directory with the .bdf/.set files
definput.keyvals.dirLOCS = '.';       % directory with the eeglab locations file
definput.keyvals.referenceChannels = [33,34]; % Indicies of reference channels in BDF file
definput.keyvals.nChansLocsFilename = 'biosemi_eloc.locs';  % filename of the eeglab locations file
definput.flags.ICA = {'','ICA'};

posdepnames = {'FS','cutoff','transitionBandwidth','maxPassbandRipple'};
[flags,kv]=ltfatarghelper(posdepnames,definput,varargin);

% switch nargin
%     case 3
%         FS = 512;  % In samps/sec
%         locutoff = 0.5;  % In Hz
%         hicutoff = 20;  % In Hz
%         transitionBandwidth = 1;  % In Hz
%         maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
%     
%     case 4
%         FS = varargin{1};
%         locutoff = 0.5;  % In Hz
%         hicutoff = 20;  % In Hz
%         transitionBandwidth = 1;  % In Hz
%         maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
% 
%     case 5
%         FS = varargin{1};
%         locutoff = varargin{2}(1);
%         hicutoff = varargin{2}(2);
%         transitionBandwidth = 1;  % In Hz
%         maxPassbandRipple = 0.0002;  % The default value used when calling 'pop_firwsord.m' GUI is 0.0002
%         
%     case 6
%         FS = varargin{1};
%         locutoff = varargin{2}(1);
%         hicutoff = varargin{2}(2);
%         transitionBandwidth = varargin{3}(1);
%         maxPassbandRipple = varargin{3}(2);
%         
%     otherwise
%         error('DKR: Invalid number of input parameters.')
% end

%**************************************************************************
% Open EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%--------------------------------------------------------------------------
% Load BDF file
bdfFilepath = fullfile(kv.dirBDF,subjID,'');  % NOTE: this path is VERY SPECIFIC TO THE USAGE
bdfFilename = ['SpExCue_', kv.experimentString, '_', subjID, '_', expDate, '.bdf'];  % NOTE: this filename is VERY SPECIFIC TO THE USAGE
EEG = pop_biosig(fullfile(bdfFilepath, bdfFilename), 'ref',kv.referenceChannels ,'refoptions',{'keepref' 'off'});
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
EEG = eeg_checkset(EEG);

%--------------------------------------------------------------------------
% Remove unused channels
EEG = pop_select(EEG,'nochannel',{'EXG4', 'EXG5' 'EXG6' 'EXG7' 'EXG8'});
EEG = pop_chanedit(EEG, 'load',{fullfile(kv.dirLOCS,kv.nChansLocsFilename) 'filetype' 'autodetect'});  % Set channel locations
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset(EEG);
fnprefix = [kv.experimentString '_subj', subjID];
[ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'savenew',fullfile(bdfFilepath,[fnprefix,'_raw']),'gui','off'); 

%--------------------------------------------------------------------------
% Resample the data to a specified sampling rate (DEFAULT: 512 samps/sec)
EEG = pop_resample(EEG, kv.FS);
fn = [fnprefix,'_resamp' num2str(kv.FS)];
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'savenew',fullfile(bdfFilepath,fn),'gui','off'); 

%--------------------------------------------------------------------------
% Filter the data
  
KaiserWindowBeta = pop_kaiserbeta(kv.maxPassbandRipple);% Calculate the filter parameters
filtOrder = pop_firwsord('kaiser', kv.FS, kv.transitionBandwidth, kv.maxPassbandRipple);

EEG = pop_firws(EEG, 'fcutoff', kv.cutoff, 'ftype', 'bandpass', 'wtype', 'kaiser', 'warg', KaiserWindowBeta, 'forder', filtOrder, 'minphase', 0);
fn = [fn,'_filt', num2str(kv.cutoff(2))];
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'savenew',fullfile(bdfFilepath,fn),'gui','off'); 

%--------------------------------------------------------------------------
% Run the ICA analysis
if flags.do_ICA
  [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
  EEG = eeg_checkset( EEG );
  EEG = pop_runica(EEG, 'extended',1,'interupt','on');
  [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
  fn = [fn,'_ICAraw'];
  EEG = pop_saveset(EEG,  'filename', fn, 'filepath', bdfFilepath);
end

%--------------------------------------------------------------------------
% Close EEGlab
close all

%--------------------------------------------------------------------------
% Output and save the data parameters 
paramsOut.timestamp = clock;
paramsOut.subjID = subjID;
paramsOut.expDate = expDate;
paramsOut.kv = kv;
paramsOut.flags = flags;
% paramsOut.FS = FS;
% paramsOut.locutoff = locutoff;
% paramsOut.hicutoff = hicutoff;
% paramsOut.transitionBandwidth = transitionBandwidth;
% paramsOut.maxPassbandRipple = maxPassbandRipple;
paramsOut.filterOrder = filtOrder;
paramsOut.KaiserWindowBeta = KaiserWindowBeta;

filename = fullfile(bdfFilepath, [kv.experimentString, '_subj', subjID, '_paramsOut', num2str(paramsOut.timestamp(1)), '_', num2str(paramsOut.timestamp(2)), '_', num2str(paramsOut.timestamp(3)), '_', num2str(paramsOut.timestamp(4)), '_', num2str(paramsOut.timestamp(5)), '.mat']);
save(filename, 'paramsOut');

end

function [flags,keyvals,varargout]  = ltfatarghelper(posdepnames,definput,arglist,callfun)
%LTFATARGHELPER  Parse arguments for LTFAT
%   Usage: [flags,varargout]  = ltfatarghelper(posdepnames,definput,arglist,callfun);
%
%   Input parameters:
%      posdepnames : Names of the position dependant parameters.
%      definput    : Struct to define the allowed input
%      arglist     : Commandline of the calling function (varargin)
%      callfun     : Name of calling function (optional)
%
%   Output parameters:
%      flags       : Struct with information about flags.
%      keyvals     : Struct with key / valtemplateues.
%      varargout   : The position dependant pars. properly initialized
%
%   [flags,keyvals]=LTFATARGHELPER(posdepnames,definput,arglist) assists in
%   parsing input parameters for a function in LTFAT. Parameters come in
%   four categories:
%  
%    Position dependant parameters. These must not be strings. These are
%     the first parameters passed to a function, and they are really just a
%     short way of specifying key/value pairs. See below.
%
%    Flags. These are single string appearing after the position dependant
%     parameters.
%
%    Key/value pairs. The key is always a string followed by the value,
%     which can be anything.
%
%    Expansions. These appear as flags, that expand into a pre-defined list
%     of parameters.  This is a short-hand way of specifying standard sets of
%     flags and key/value pairs.
%
%   The parameters are parsed in order, so parameters appearing later in
%   varargin will override previously set values.
%
%   The following example for calling LTFATARGHELPER is taken from DGT:
%  
%     definput.keyvals.L=[];
%     definput.flags.phase={'freqinv','timeinv'};
%     [flags,kv]=ltfatarghelper({'L'},definput,varargin);
%
%   The first line defines a key/value pair with the key 'L' having an
%   initial value of [] (the empty matrix).
%
%   The second line defines a group of flags by the name of phase.  The
%   group phase contains the flags 'freqinv' and 'timeinv', which can
%   both be specified on the command line by the user. The group-name
%   phase is just for internal use, and does not appear to the user. The
%   flag mentioned first in the list will be selected by default, and only
%   one flag in a group can be selected at any time. A group can contain as
%   many flags as desired.
%  
%   The third line is the actual call to LTFATARGHELPER which defines the
%   output flags and kv.  The input {'L'} indicates that the value of
%   the parameter 'L' can also be given as the very first value in
%   varargin.
%
%   The output struct kv contains the key/value pairs, so the value
%   associated to 'L' is stored in kv.L.
%
%   The output struct flags contains information about the flags choosen
%   by the user. The value of flags.phase will be set to the selected flag
%   in the group phase and additionally, the value of flags.do_timeinv
%   will be 1 if 'timeinv' was selected and 0 otherwise, and similarly for
%   'freqinv'. This allows for easy checking of selected flags.
%
%   See also: ltfatgetdefaults, ltfatsetdefaults
%
%   Url: http://ltfat.github.io/doc/ltfatarghelper.html

% Copyright (C) 2005-2015 Peter L. Soendergaard <peter@sonderport.dk>.
% This file is part of LTFAT version 2.1.1
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

persistent TF_CONF;

if isempty(TF_CONF)

%  basepath=which('ltfatarghelper');
%  % Kill the function name and comp from the path.
%  basepath=basepath(1:end-22);
%  % add the base path
%  addpath(basepath);
%  ltfatstart;

  TF_CONF.fundefs = struct;
end;

if ischar(posdepnames)
  % Special interface needed for ltfatsetdefaults and ltfatgetdefaults,
  % activated when first argument is a string.

  % First input  argument, posdepnames, is a string, one of the options
  % in the "switch" section below
  % Second input argument, definput,    is a function name to get or set
  % Third  input argument, arglist ,    is a cell-array with options to set.
  
  switch(lower(posdepnames))
   case 'get'
    if isfield(TF_CONF.fundefs,definput)
      flags=TF_CONF.fundefs.(definput);
    else
      flags={};
    end;
   case 'set'
    TF_CONF.fundefs.(definput)=arglist;
   case 'all'
    flags=TF_CONF.fundefs;
   case 'clearall'
    TF_CONF.fundefs=struct; 
  end;
  return
end;

if nargin<4
  f=dbstack;  
  callfun=f(2).name;
end;

nposdep=numel(posdepnames);

% Resolve import specifications BEFORE adding our own specifications.
if isfield(definput,'import')
  for imp = definput.import;
    definput=feval(['arg_',imp{1}],definput);
  end;
end;

if isfield(definput,'flags')
  defflags=definput.flags;
else
  defflags=struct;
end;

if isfield(definput,'keyvals')
  defkeyvals=definput.keyvals;
else
  defkeyvals=struct;
end;

if isfield(definput,'groups')
  groups=definput.groups;
else
  groups=struct;
end;

total_args = numel(arglist);

% Determine the position of the first optional argument.
% If no optional argument is given, return nposdep+1
first_str_pos = 1;
while first_str_pos<=total_args && ~ischar(arglist{first_str_pos}) 
  first_str_pos = first_str_pos +1;    
end;

% If more than nposdep arguments are given, the first additional one must
% be a string
if (first_str_pos>nposdep+1)
  error('%s: Too many input arguments',upper(callfun));
end;

n_first_args=min(nposdep,first_str_pos-1);

keyvals=defkeyvals;      

% Copy the given first arguments
for ii=1:n_first_args
  keyvals.(posdepnames{ii})=arglist{ii};
end;

% Initialize the position independent parameters.
% and create reverse mapping of flag -> group
flagnames=fieldnames(defflags);
flags=struct;
% In order for flags to start with a number, it is necessary to add
% 'x_' before the flag when the flags are used a field names in
% flagreverse. Externally, flags are never used a field names in
% structs, so this is an internal problem in ltfatarghelper that is
% fixed this way.
flagsreverse=struct;
for ii=1:numel(flagnames)
  name=flagnames{ii};
  flaggroup=defflags.(name);
  flags.(name)=flaggroup{1};
  for jj=1:numel(flaggroup)
    flagsreverse.(['x_', flaggroup{jj}])=name;
    flags.(['do_',flaggroup{jj}])=0;
  end;
  flags.(['do_',flaggroup{1}])=1;
end;

%Get the rest of the arguments
restlist = arglist(first_str_pos:end);

%Check for default arguments
if isfield(TF_CONF.fundefs,callfun)
  s=TF_CONF.fundefs.(callfun);
  restlist=[s,restlist];
end;

% Check for import defaults
if isfield(definput,'importdefaults')
  % Add the importdefaults before the user specified arguments.
  restlist=[definput.importdefaults,restlist];
end;

while ~isempty(restlist)
  argname=restlist{1};
  restlist=restlist(2:end);  % pop
  found=0;
  
  % Is this name a flag? If so, set it
  if isfield(flagsreverse,['x_',argname])
    % Unset all other flags in this group
    flaggroup=defflags.(flagsreverse.(['x_',argname]));
    for jj=1:numel(flaggroup)
      flags.(['do_',flaggroup{jj}])=0;
    end;
    
    flags.(flagsreverse.(['x_',argname]))=argname;
    flags.(['do_',argname])=1;
    found=1;
  end;
  
  % Is this name the key of a key/value pair? If so, set the value.
  if isfield(defkeyvals,argname)      
    keyvals.(argname)=restlist{1};
    restlist=restlist(2:end);
    found=1;
  end;
  
  % Is this name a group definition? If so, put the group in front of the parameters
  if isfield(groups,argname)
    s=groups.(argname);
    restlist=[s,restlist];
    found=1;
  end;
  
  % Is the name == 'argimport'
  if strcmp('argimport',argname)   
    fieldnames_flags= fieldnames(restlist{1});  
    fieldnames_kvs  = fieldnames(restlist{2});        
    for ii=1:numel(fieldnames_flags)
      importname=fieldnames_flags{ii};
      flags.(importname)=restlist{1}.(importname);
    end;
    for ii=1:numel(fieldnames_kvs)
      importname=fieldnames_kvs{ii};
      keyvals.(importname)=restlist{2}.(importname);
    end;      
    restlist=restlist(3:end);
    found=1;
  end;
  
  if found==0
    if ischar(argname)
      error('%s: Unknown parameter: %s',upper(callfun),argname);
    else
      error('%s: Parameter is not a string, it is of class %s',upper(callfun),class(argname));          
    end;      
  end;
  
  %ii=ii+1;
end;

% Fill varargout

varargout=cell(1,nposdep);
for ii=1:nposdep
    varargout(ii)={keyvals.(posdepnames{ii})};
end;

end