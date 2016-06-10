classdef tdt < handle
% tdtObject = tdt(paradigmType, sampleRate, scaling, trigDuration, ...
%                 noiseAmp, figNum=99999)
%
% Creates a new tdt object.
% 
% Inputs:
% ----------------------------------------------------------------------------
%
% paradigmType: 'playback_1channel' or 'playback_2channel'. 
%    playback_1channel: one-channel stimuli up to 8.3842E6 samples
%    playback_2channel: two-channel stimuli up to 4.1921E6 samples
%
% requestedSampleRate: must be 48, 24, or 12, for 48828.125 Hz, 24414.0625 Hz,
% or 12207.03125 Hz, respectively. Please make note of the non-standard sample
% rates.
%
% scaling: controls the bounds defining full scale, in volts. Specify as a 2
% element vector for different scaling per channel. On the RP2, this should not
% exceed 10 V. If one number is specified, then this scaling value is used for
% both channels. Note: Even when paradigmType=='playback_1channel', a
% two-element vector can be specified; this is because the noise in the
% opposite ear can still be played.
% 
% trigDuration: the duration, in seconds, that each event signal should last.
% Default: 5E-3 s
%
% noiseAmp: the RMS amplitude of the noise relative to 1V RMS, in dB; e.g.,
% noiseRMSVoltage = 10^(noiseAmpInDecibels / 20). The actual output RMS depends
% on the value set for "scaling". Accepts -Inf as an input for no noise (the
% default).
%
% figNum: by default, creates the ActiveX figure as figure number 99999;
% specify an integer argument if for some reason you want another value. There
% is no good reason to change this setting, unless you have a figure number
% 99999 used for something else in your code (or this value conflicts with some
% other dummy figure in another software package).
%
%
% Outputs:
% ----------------------------------------------------------------------------
%
% tdtObject: an object of type "tdt", with the following properties and
%   methods:
%
%   Properties:
%
%       sampleRate - the real sample rate at which the RP2/RP2.1 operates
%
%       bufferSize - the maximum number of samples that can be handled by the
%       circuit without further input from the user. The exact value will vary
%       depending on the paradigm type.
%
%       channel1Scale / channel2Scale - the scaling value x mapping floating
%       point values between [-1,1] to [-x,x] for channel 1/2 (in Volts)
%
%       noise1RMS / noise2RMS - the RMS value of the background noise (in V)
%
%       status: a status string describing the current state of the circuit
%        
%       stimSize: size (in samples) of the stimulus loaded on the TDT
%       
%       nChans: the number of playback channels (either 1 or 2)
%
%       trigDuration: the duration of a digital event sent via the digital out
%       port on the RP2.1
%
%   User-facing methods; (type "help tdtObj.<function_name>" for a full
%   description, where <function_name> is one of the following:
%
%       load_stimulus(audioData, [triggerInfo = [1,1] ])
%
%       play([stopAfter = obj.stimSize])
%
%       play_blocking([stopAfter = obj.stimSize])
%
%       pause()
%
%       rewind()
%
%       reset()
%
%       send_event(integerEventValue)
%
%       get_current_sample([consistencyCheck = true])
%
%       set_noise_level(noiseLevelVector)
%
%
% Version 1.2 (2015-06-09) 
% Auditory Neuroscience Lab, Boston University
% Contact: lennyv_at_bu_dot_edu

    properties(SetAccess = 'private', GetAccess='public')
        sampleRate
        channel1Scale
        channel2Scale
        noise1RMS
        noise2RMS
        status
        stimSize
        nChans
        trigDuration
    end
    
    properties(Access='private')
        RP
        f1
        bufferSize
    end

    methods
        function obj = tdt(paradigmType, requestedSampleRate, scaling, ...
                           trigDuration, noiseAmpDB, figNum)
          
            %%% sample rate check
            if nargin < 2
                error('Desired sample rate must be specified.')
            end

            if 48 == requestedSampleRate 
                rateTag = 3;
            elseif 24 == requestedSampleRate
                rateTag = 2;
            elseif 12 == requestedSampleRate 
                rateTag = 1;
            else
                error('invalid sample rate specified (must be 48, 24, 12)')
            end

            %%% voltage scaling
            if nargin < 3 
                error('Scaling must be specified.')
            end
           
            % by default, set scaling equal on both channels
            if length(scaling) < 2
                scaling(2) = scaling(1);
            end

            if nargin < 4 || isempty(trigDuration)
               trigDuration = 5E-3; % s
            end

            %%% control the background noise amplitude (dbFS)
            if nargin < 5 || isempty(noiseAmpDB)
                noiseAmpDB = [-Inf, -Inf];
            end
            
            if length(noiseAmpDB) < 2
                noiseAmpDB(2) = noiseAmpDB(1);
            end
            
            noiseAmpVolts = 10.^(noiseAmpDB ./ 20);
            

            if nargin < 6 || isempty(figNum)
                figNum = 99999;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Start ActiveX controls and hides the figure at the start of
            % each block
            obj.f1 = figure(figNum);
            set(obj.f1,'Position', [5 5 30 30], 'Visible', 'off');
            obj.RP = actxcontrol('RPco.x', [5 5 30 30], obj.f1);

            % gigabit isn't supported anymore (as of tdt version 70)
            obj.RP.ConnectRP2('USB', 1);

            %Clears all the Buffers and circuits on that RP2
            obj.RP.ClearCOF;
            %Loads the appropriate circuit, with a quick binary check to ensure
            %file versions are correct
            if strcmpi(paradigmType, 'playback_1channel')
                obj.nChans = 1;
                obj.bufferSize = 8.3842E6;
            elseif strcmpi(paradigmType, 'playback_2channel')
                obj.nChans = 2;
                obj.bufferSize = 4.1921E6;
            else
                error('paradigm type is currently unsupported.')
            end
            
            load(['bin/' paradigmType '.mat'])
            fileID = fopen(['bin/' paradigmType '.rcx']);
            temp = fread(fileID, Inf, 'int32=>int32');
            fclose(fileID);
            if any(size(temp) ~= size(binInfo)) || any(temp ~= binInfo)
                error('Version mismatch between .m and .rcx files.')
            end
            obj.RP.LoadCOFsf(['bin/' paradigmType '.rcx'], rateTag);

            % store some relevant info in the object itself

            % sample rate
            obj.sampleRate = obj.RP.GetSFreq();
            
            % trigger duration (fixed)
            obj.trigDuration = trigDuration;
            obj.RP.SetTagVal('triggerDuration', 1000 * obj.trigDuration);

            % scaling factors and RMS levels
            obj.channel1Scale = single(scaling(1));
            obj.channel2Scale = single(scaling(2));
            obj.RP.SetTagVal('chan1Scaler', obj.channel1Scale);
            obj.RP.SetTagVal('chan2Scaler', obj.channel2Scale);

            % noise
            obj.noise1RMS = single(noiseAmpVolts(1)) * ...
                            obj.channel1Scale;
            obj.noise2RMS = single(noiseAmpVolts(2)) * ...
                            obj.channel2Scale;

            % as per TDT documentation on "GaussNoise" component
            if (obj.noise1RMS > 2.1) || (obj.noise2RMS > 2.1)
                error('Noise RMS is too high. Clipping will occur.')
            end
            obj.RP.SetTagVal('chan1NoiseAmp', obj.noise1RMS);
            obj.RP.SetTagVal('chan2NoiseAmp', obj.noise2RMS);
            obj.RP.SetTagVal('chan1NoiseSeed', randi(2^15));
            obj.RP.SetTagVal('chan2NoiseSeed', randi(2^15));

            % zero tag the relevant buffers
            obj.RP.ZeroTag('audioChannel1');
            if obj.nChans == 2
                obj.RP.ZeroTag('audioChannel2');
            end
            obj.RP.ZeroTag('triggerIdx');
            obj.RP.ZeroTag('triggerVals');

            % now attempt to actually run the circuit
            obj.RP.Run;
            if obj.RP.GetStatus ~= 7
                obj.RP.close();
                error('TDT connection error. Try rebooting the TDT.');
            end
            
            % do an "initial reset" of the buffers to fix indexing on
            % source buffers
            obj.RP.SoftTrg(3);
            
            % display some status information to the user
            fprintf('Channel 1, [-1.0, 1.0] --> [-%2.4f, %2.4f] V\n', ...
                 obj.channel1Scale, obj.channel1Scale);
            fprintf('Channel 2, [-1.0, 1.0] --> [-%2.4f, %2.4f] V\n', ...
                 obj.channel2Scale, obj.channel2Scale);
            fprintf('Channel 1, masking noise RMS (V) = %2.4f\n', ...
                 obj.channel1Scale * obj.noise1RMS);
            fprintf('Channel 2, masking noise RMS (V) = %2.4f\n', ...
                 obj.channel2Scale * obj.noise2RMS);
             
            obj.stimSize = 0;
            obj.status = sprintf('No stimulus loaded.');
        end


        function load_stimulus(obj, audioData, triggerInfo)
        % tdt.load_stimulus(audioData, triggerInfo) 
        %
        % function to load stimulus and triggers to TDT circuit.
        %
        % audioData: a 1D or 2D column array specifying audio data ** See note
        % 1
        %
        % triggerInfo: an n x 2 array specifying index and value tuples to send
        % a digital "word" value at the specified sample of playback. ** see
        % note 2
        %
        % note 1: audioData must be limited to [-1, 1], and must be in sample x
        % channel format (the default for Matlab); it will be converted to TDT-
        % friendly format in this function.
        %
        % This function will downconvert the arrays to single-precision prior
        % to writing to the TDT if they are not already stored as single
        % precision.
        %
        % note 2: Trigger samples should be specified using Matlab index style,
        % i.e., the first sample of audio is sample 1. Permissible trigger
        % values will vary by device; e.g., on the RP2, values should be <=
        % 255, corresponding to the 8 bit precision of the digital output. If a
        % value is > 255, only the least significant 8 bits are used. Duration 
        % should be specified in seconds. It should be obvious that the values
        % need to be non-negative.
        %
        % last updated: 2015-06-09, LAV, lennyv_at_bu_dot_edu
            
            
            %%%%%%%%%%%%%%%%%%%%
            % input validation %
            %%%%%%%%%%%%%%%%%%%%

            if nargin < 3
                triggerInfo = [];
            else
                if (size(triggerInfo, 2) ~= 2) || ...
                    (length(size(triggerInfo)) ~= 2) || ...
                    any(triggerInfo(:) < 0)
            
                    error(['triggerInfo must be specified as ',...
                           '[idx, val], array, and the values '...
                           'should all be positive.'])
                end
            end
            
            % convert down to single precision floating point, 
            % since that's what the TDT natively uses for DAC
            if ~isa(audioData, 'single')
                audioData = single(audioData);
            end

            % stimulus size checks

            if size(audioData, 1) > obj.bufferSize
                error(['Stimulus should be <= %d samples long. ' ...
                       'Shorten the stimulus and try again.'], ...
                       obj.bufferSize)
            end
            
            if size(audioData, 2) ~= obj.nChans
                error(['Number of columns in audioData should ' ...
                       'match number of channels specified (%d)'], ...
                       obj.nChans)
            end
            
            if size(triggerInfo, 1) > 2197
                error(['Circuit can only support a maximum of 2197 ', ...
                       'trigger values. Reduce the number of ', ...
                       'triggers specified in triggerInfo.'])
            end

            if ~isempty(triggerInfo)
                triggerIdx = int32(triggerInfo(:, 1));
                triggerVals = int32(triggerInfo(:, 2));
            else
                % send a single trigger of value 1 at start of playback
                triggerIdx = int32(1);
                triggerVals = int32(1);
            end
            
            if any(triggerVals < 0)
                error('Trigger values should be non-negative.')
            end
            
            if (any(triggerIdx > size(audioData, 1)))
                error('Trigger index must be smaller than stimulus size.')
            end
            
            if any(triggerIdx < 1)
                error('Trigger index should be positive.')
            end

            if any(abs(audioData) > 1)
                error('All audio data must be scaled between -1.0 and 1.0.')
            end

            %%%%%%%%%%%%%%%%%%%%%
            % write data to TDT %
            %%%%%%%%%%%%%%%%%%%%%
           
            % hack - the WriteTagVEX methods don't like single value inputs
            % also correct for 1 sample difference in index
            triggerIdx = [triggerIdx - 1; -1];
            triggerVals = [triggerVals; 0];
            
            % reset buffer indexing and zeroTag everything
            obj.reset_buffers(true)
            % size + 1 is intentional on next line
            obj.RP.SetTagVal('stimSize', size(audioData, 1)+1); 
            obj.stimSize = size(audioData, 1);
            
            fprintf('Writing to channel 1 buffer...\n')
            curStatus = obj.RP.WriteTagVEX('audioChannel1', 0, 'F32',...
                                           audioData(:, 1));
            if ~curStatus
                error('Error writing to audioChannel1 buffer.')
            end
            
            if obj.nChans == 2
                fprintf('Writing to channel 2 buffer...\n')
                curStatus = obj.RP.WriteTagVEX('audioChannel2', 0, 'F32',...
                                            audioData(:, 2));
                if ~curStatus
                    error('Error writing to audioChannel2 buffer.')
                end
            end
            
            fprintf('Writing to triggerIdx buffer...\n')
            curStatus = obj.RP.WriteTagVEX('triggerIdx', 0, 'I32',...
                                           triggerIdx);
            if ~curStatus
                error('Error writing to triggerIdx buffer.')
            end
            
            fprintf('Writing to triggerVals buffer...\n')
            curStatus = obj.RP.WriteTagVEX('triggerVals', 0, 'I32',...
                                           triggerVals);
            if ~curStatus
                error('Error writing to triggerVals buffer.')
            end
            
            fprintf('Stimulus loaded.\n')
        end


        function play(obj, stopAfter)
        % tdt.play(stopAfter)
        %
        % Plays the contents of the audio buffers on the TDT.
        %
        % Inputs:
        % --------------------------------------------------------------------
        % stopAfter - the sample number at which playback should cease. If not
        % specified, playback will continue until the end of the stimulus is 
        % reached.
        %
        % last updated: 2015-03-11, LAV, lennyv_at_bu_dot_edu

            if obj.stimSize == 0
                error(['No stimulus loaded.'])
            end

            if nargin < 2
                stopAfter = obj.stimSize;
            end
            if stopAfter < obj.get_current_sample()
                error(['Buffer index already passed desired stop point. ' ...
                       'Did you mean to rewind the buffer first?'])
            end
            stat = obj.RP.SetTagVal('stopSample', stopAfter);
            if ~stat
                error('Error setting stop sample.')
            end
            obj.RP.SoftTrg(1);
            obj.status = sprintf('playing then stopping at buffer index %d',...
                                 stopAfter);
        end


        function pause(obj)
        % tdt.pause()
        %
        % Pauses playback on the TDT.
        %
        % last updated: 2015-04-03, LAV, lennyv_at_bu_dot_edu
            
            stat = obj.RP.SetTagVal('stopSample', 0);
            if ~stat
                error('Error setting stop sample.')
            end
            pause(0.02);
            currentSample = obj.get_current_sample();
            obj.status = sprintf('stopped at buffer index %d', currentSample);
        end
       

        function play_blocking(obj, stopAfter)
        % tdt.play_blocking(stopAfter)
        %
        % Plays the contents of the audio buffers on the TDT and holds up
        % Matlab execution while doing so.
        %
        % Inputs:
        % --------------------------------------------------------------------
        % stopAfter - the sample number at which playback should cease. If not
        % specified, playback will continue until the end of the stimulus is 
        % reached.
        %
        % version added: 1.1
        % last updated: 2015-04-06, LAV, lennyv_at_bu_dot_edu

            if obj.stimSize == 0
                error(['No stimulus loaded.'])
            end

            if nargin == 1
                stopAfter = obj.stimSize;
            end
            
            if stopAfter < obj.get_current_sample()
                error(['Buffer index already passed desired stop point. ' ...
                       'Did you mean to rewind the buffer first?'])
            end
            
            stat = obj.RP.SetTagVal('stopSample', stopAfter);
            if ~stat
                error('Error setting stop sample.')
            end
            
            obj.RP.SoftTrg(1);
            fprintf('Playing stimulus in blocking mode...')
            try
                currentSample = obj.RP.GetTagVal('chan1BufIdx');
                while currentSample < stopAfter 
                    currentSample = obj.RP.GetTagVal('chan1BufIdx');
                    pause(0.1);
                end
            catch ME
                obj.pause()
                fprintf(['\n' obj.status]);
                throw(ME);
            end
            fprintf('done.\n')
            obj.rewind();
        end


        function rewind(obj)
        % tdt.rewind()
        %
        % Rewinds the buffer without clearing it. Useful when new audio
        % does not need to be loaded into the TDT.
        %
        % last updated: 2015-03-11, LAV, lennyv_at_bu_dot_edu
        
            obj.reset_buffers(false);
            currentSample = obj.get_current_sample();
            obj.status = sprintf('stopped at buffer index %d', currentSample);
        end
       

        function reset(obj)
        % tdt.reset()
        %
        % Rewinds the buffer and sets all values in the buffer to 0. 
        %
        % last updated: 2015-03-11, LAV, lennyv_at_bu_dot_edu

            obj.reset_buffers(true);
            obj.stimSize = 0;
            currentSample = obj.get_current_sample();
            obj.status = sprintf('stopped at buffer index %d', currentSample);
        end
        
        function send_event(obj, eventVal)
        % tdt.send_event(eventVal)
        %
        % Sends an arbitrary integer event to the digital out port on the TDT.
        % Timing will not be sample-locked in any way.
        %
        % last updated: 2015-03-11, LAV, lennyv_at_bu_dot_edu
        
            statusVal = obj.RP.SetTagVal('arbitraryEvent', eventVal);
            if ~statusVal
                error('Event could not be written.')
            end
            pause(0.01);
            obj.RP.SoftTrg(4);
        end
       

        function set_noise_level(obj, noiseAmpDB)
        % tdt.set_noise_level(noiseAmpDB)
        %
        % Controls the noise levels in each ear. Specify noiseAmpDB as a scalar
        % to set both channels to the same value, or a 2-element vector. Set 
        % the value of a channel to -Inf to turn noise off.
        %
        % Version added: 1.1
        %
        % last updated: 2015-03-30, LAV, lennyv_at_bu_dot_edu
        
            if length(noiseAmpDB) < 2
                noiseAmpDB(2) = noiseAmpDB(1);
            end
            
            noiseAmpVolts = 10.^(noiseAmpDB ./ 20);

            newNoise1RMS = single(noiseAmpVolts(1)) * ...
                               obj.channel1Scale;
            newNoise2RMS = single(noiseAmpVolts(2)) * ...
                                   obj.channel2Scale;
                               
            obj.RP.SetTagVal('chan1NoiseAmp', newNoise1RMS);
            obj.noise1RMS = newNoise1RMS;
            obj.RP.SetTagVal('chan2NoiseAmp', newNoise2RMS);
            obj.noise2RMS = newNoise2RMS;
            pause(0.1);
        end
        

        function [currentSample1, trigBufSample1] = get_current_sample(obj, checks)
        % [audioIdx, triggerIdx] = tdt.get_current_sample(checks)
        %
        % Gets the current buffer position for the audio stimuli (output 1) and
        % for triggers (output 2).
        %
        % An error is raised if the audio buffers or the trigger buffers become
        % misaligned.
        %
        % last updated: 2015-04-03, LAV, lennyv_at_bu_dot_edu
            if nargin == 1
                checks = true;
            end
        
            currentSample1= obj.RP.GetTagVal('chan1BufIdx');
            if checks && obj.nChans == 2
                currentSample2 = obj.RP.GetTagVal('chan2BufIdx');
                if currentSample1 ~= currentSample2
                    obj.reset_buffers(false);
                    error(['Audio buffers are misaligned (%d/%d.).',...
                        'Buffers reset, but not cleared.'], ...
                        currentSample1, currentSample2)
                end
            end
            
            trigBufSample1 = obj.RP.GetTagVal('trigIdxBufferIdx');
            if checks
                trigBufSample2 = obj.RP.GetTagVal('trigValBufferIdx');
                if trigBufSample1 ~= trigBufSample2
                    obj.reset_buffers(false);
                    error(['Trigger buffers are misaligned (%d/%d.)',...
                        'Buffers reset, but not cleared.'],...
                        trigBufSample1, trigBufSample2)
                end
            end
        end
    end
   

    methods(Access='private')
        function reset_buffers(obj, clearBuffer)
        % tdt.reset_buffers(clearBuffer)
        %
        % Resets and optionally zero-tags the buffers in the circuit. Not meant
        % to be called by the end user.
        %
        % Inputs:
        % --------------------------------------------------------------------
        % 
        % clearBuffer: boolean. If true, will zero-tag (i.e., erase) buffers by
        % setting them to 0. Otherwise just resets the all buffer indexing to
        % 0.
        %
        % last updated: 2015-03-11, LAV, lennyv_at_bu_dot_edu           
            
            obj.RP.SoftTrg(2);
            pause(0.01);
            
            if clearBuffer
                obj.RP.ZeroTag('audioChannel1');
                obj.RP.ZeroTag('audioChannel2');
                obj.RP.ZeroTag('triggerIdx');
                obj.RP.ZeroTag('triggerVals');
                obj.RP.SetTagVal('stopSample', 0);
                obj.RP.SetTagVal('stimSize', 1);
                obj.stimSize = 0;
            end
            
            obj.RP.SoftTrg(3);
            pause(0.01);
            currentSample = obj.get_current_sample();
            if currentSample ~= 0
                error('Buffer rewind error.');
            end
            
            obj.status = sprintf('stopped at buffer index %d',...
                                 currentSample);
        end
        
        function delete(obj)
        % tdt.delete()
        %
        % cleanly back out and close the TDT when the object is deleted. Not
        % meant to be called by the user.
        %
        % last updated: 2015-03-11, LAV, lennyv_at_bu_dot_edu

            obj.reset_buffers(true);
            obj.RP.Halt;
            pause(0.01);
            obj.RP.ClearCOF;
            close(obj.f1);
            obj.status = sprintf('Not connected.');
        end
        
    end
end
