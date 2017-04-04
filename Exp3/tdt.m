classdef tdt < handle
% ANL TDT interface code, v1.6 (2016-04-22)
%
% tdtObject = tdt(paradigmType, sampleRate, scaling, ...
%                 'triggerDuration', 0.005, 'buttonHoldDuration', 0.2, ...
%                 'xorVal', 0, 'figNum', 9999)
%
% Creates a new object to interact with the TDT RP2.1 or the TDT RZ6.
%
% Inputs:
% ----------------------------------------------------------------------------
%   Required
%   '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%   paradigmType: character string, with one of the following values:
%      'playback_1channel':
%          one-channel stimuli up to 8.38E6 (RP2.1) or 1.67E7 (RZ6) samples
%      'playback_2channel':
%          two-channel stimuli up to 4.19E6 (RP2.1) or 8.38E6 (RZ6) samples
%      'playback_2channel_16bit':
%          two-channel stimuli up to 8.38E6 (RP2.1) or 1.67E7 (RZ6)
%          samples, but with slight loss of precision due to 32-->16 bit
%          conversion during transfer
%
%   sampleRate: must be 48, 24, or 12, for 48828.125 Hz, 24414.0625 Hz, or
%   12207.03125 Hz, respectively. If using an RZ6, can use 97 to get a
%   sample rate of 97656.250 Hz. Please make note of the non-standard
%   sample rates.
%
%   scaling: controls the bounds defining full scale, in volts. Specify as
%   a 2 element vector for different scaling per channel. An error is
%   raised if this exceeds 10 V. If one number is specified, then this
%   scaling value is used for both channels when paradigm type is
%   "playback_2channel". When paradigm type is "playback_1channel" and only
%   one number is specified, this value is used for channel 1; the opposite
%   channel is set to 0. Thus, for monaural playback, use paradigmType =
%   "playback_1channel" and specify a 1-channel scaler. For diotic
%   playback, use paradigmType = "playback_1channel" and specify a
%   2-channel scaler.
%   
%   Optional Keyword Arguments (new in v1.6)
%   '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%   triggerDuration: the duration, in seconds, that each event signal should
%   last. Default: 5E-3 s
%
%   buttonHoldDuration: the duration, in s, that is required before a
%   continuous button press registers as a new event. e.g., if a button is
%   held down for 0.3s, buttonHoldDuration = 0.2, two button presses will
%   be registered. If buttonHoldDuration = 0.090, four button presses
%   should be registered (t = 0, t = 0.090, t = 0.180, t = 0.270). Defaults
%   to 0.2s.
%
%   xorVal: an integer value that button box inputs will be xor-ed with;
%   helps accomodate button boxes with logic normally high and those with
%   logic normally low. i.e., if a 4-button box is logic high when no
%   buttons are pressed, it will report "15" (sum(2.^[0, 1, 2, 3]))...if
%   button 1 is pressed, it will report "14" (sum(2.^[1, 2, 3])). To
%   "invert" the logic and store the correct values, set xor to 15. You may
%   have to play with this value depending on the configuration of your
%   digital inputs to the TDT. If not specificed, autoconfiguration is
%   attempted.
%
%   figNum: by default, creates the ActiveX figure as figure number 99999;
%   specify an integer argument if for some reason you want another value.
%   There is no good reason to change this setting, unless you have a
%   figure number 99999 used for something else in your code (or this value
%   conflicts with some other dummy figure in another software package).
%   
%   Specifying positional arguments is still supported for the sake of
%   backwards compatibility with existing experiments using earlier
%   versions of this code. Arguments 4-7 should then be specified in the
%   order: triggerDuration (seconds), buttonHoldDuration (seconds), xorVal
%   (int), figNum (int).
%
%
% Outputs:
% ----------------------------------------------------------------------------
%   tdtObject: an object of type "tdt", with the following properties and
%     methods:
%
%     Properties:
%
%       sampleRate - the real sample rate at which the RP2/RP2.1/RZ6
%       operates
%
%       bufferSize - the maximum number of samples that can be handled by
%       the circuit without further input from the user. The exact value
%       will vary depending on the paradigm type.
%
%       channel1Scale / channel2Scale - the scaling value x mapping
%       floating point values between [-1,1] to [-x,x] for channel 1/2 (in
%       Volts)
%
%       status: a status string describing the current state of the circuit
%        
%       stimSize: size (in samples) of the stimulus loaded on the TDT
%       
%       nChans: the number of playback channels (either 1 or 2)
%    
%       paradigmType: the selected paradigm type used to generate this
%       object
%
%       triggerDuration: the duration of a digital event sent via the digital
%       out port on the RP2.1
%
%       buttonHoldDuration: the button hold duration (see description of
%       input argument)
%
%
%     User-facing methods; (type "help <tdtObj>.<function_name>" for a full
%     description, where <tdtObj> is the name of the tdt object generated
%     by calling tdt(...), and <function_name> is one of the following
%     function names:
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
%       get_button_presses()
%
%       get_current_sample([consistencyCheck = true])
%
%     e.g.: 
%       >> myTDT = tdt('playback_1channel', 48, [1, 1]); 
%       >> help myTDT.get_button_presses 
%       >> help myTDT.play_blocking
%
%   See the example.m file to see more examples of intended/potential usage.
% ----------------------------------------------------------------------------
% Release notes:
%   As of v1.6, less frequently used arguments are now able to be specified
%   using Matlab keyword arguments: 'triggerDuration', 'buttonHoldDuration', 
%   'xorVal', and 'figNum'. These optional arguments are numeric, and
%   default to the same values as in earlier versions. Positional arguments
%   are still supported.
%
%   TDT RZ6 support added as of v1.6 (2016-04-20). The setup for the RZ6
%   assumes bytes A + C are set up as outputs, and byte B is set up as an
%   input. Stimulus scaling is handled through the programmable attenuator
%   component (see TDT System3 documentation for details). Also note that
%   if comparing output from RP2/RP2.1 and RZ6, you should see a 1-sample
%   difference in event value output timing. This is due to the slightly
%   different ADC group delays on the two systems (30 samples for
%   RP2/RP2.1, 31 samples for the RZ6).
%
%   As of v1.3, paramters relating to TDT-generated noise have been
%   dropped; 1) addition of the button press detection pushes the circuits
%   over the RP2.1's capabilities when noise is being generated by the TDT
%   as well; 2) it makes more sense for users to generate their own noise
%   and combine it with their stimuli in Matlab itself.
% -----------------------------------------------------------------------------
% Version 1.6 (2016-04-21) Auditory Neuroscience Lab, Boston University
% Contact: lennyv_at_bu_dot_edu

    properties(SetAccess = 'private', GetAccess='public')
        sampleRate
        deviceType
        channel1Scale
        channel2Scale
        status
        stimSize
        nChans
        triggerDuration
        buttonHoldDuration
        paradigmType
    end
    
    properties(Access='private')
        RP
        f1
        hiddenFigure
        bufferSize
        xorVal
    end

    methods
        function obj = tdt(paradigmType, requestedSampleRate, scaling, ...
                           varargin)
                       
            rateTag = obj.parse_inputs(paradigmType, requestedSampleRate,...
                                       scaling, varargin{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % try connecting to an RZ6 via Optibit first...if that fails,
            % then fall back to RP2. As per activeX manual, 'GB' is correct
            % for optibit interfaces.
            useRZ6 = false;
            deviceFound = false;
            for devNum = 0:4
                if obj.RP.ConnectRZ6('GB', devNum);
                    useRZ6 = true;
                    deviceFound = true;
                    fprintf(1, 'RZ6 found (device number: %d)!\n', devNum);
                    break
                end
            end
            
            if ~useRZ6
                for connType = {'USB', 'GB'};
                    for devNum = 0:4
                        if obj.RP.ConnectRP2('USB', devNum);
                            fprintf(1, ['RP2.1 found ',...
                                    '(%s, device number: %d)!\n'], ...
                                    char(connType), devNum);
                            deviceFound = true;
                            break
                        end
                    end
                    if deviceFound
                        break
                    end
                end
            end
            if ~useRZ6 && rateTag == 4
                error(['Cannot use sample rate setting 97 with RP2.1.', ...
                       'Use 48, 24, or 12 instead.']);
            end

            %Clears all the Buffers and circuits on the processor
            obj.RP.ClearCOF;
            %Loads the appropriate circuit, with a quick binary check to ensure
            %file versions are correct
            if useRZ6
                if strcmpi(obj.paradigmType, 'playback_1channel')
                    obj.nChans = 1;
                    obj.bufferSize = 1.67E7;
                elseif strcmpi(obj.paradigmType, 'playback_2channel')
                    obj.nChans = 2;
                    obj.bufferSize = 8.38E6;
                elseif strcmpi(obj.paradigmType, 'playback_2channel_16bit')
                    obj.nChans = 2;
                    obj.bufferSize = 1.67E7;
                else
                    error('paradigm type is currently unsupported.')
                end
                circuitName = ['bin/' obj.paradigmType '_button_rz6'];
            else
                if strcmpi(obj.paradigmType, 'playback_1channel')
                    obj.nChans = 1;
                    obj.bufferSize = 8.38E6;
                elseif strcmpi(obj.paradigmType, 'playback_2channel')
                    obj.nChans = 2;
                    obj.bufferSize = 4.19E6;
                elseif strcmpi(obj.paradigmType, 'playback_2channel_16bit')
                    obj.nChans = 2;
                    obj.bufferSize = 8.38E6;
                else
                    error('paradigm type is currently unsupported.')
                end
                circuitName = ['bin/' obj.paradigmType '_button'];
            end
            
            load([circuitName '.mat'], 'binInfo')
            fileID = fopen([circuitName '.rcx']);
            temp = fread(fileID, Inf, 'int32=>int32');
            fclose(fileID);
            if any(size(temp) ~= size(binInfo)) || any(temp ~= binInfo)
                error('Version mismatch between .m and .rcx files.')
            end
            
            % load circuit into TDT memory
            obj.RP.LoadCOFsf([circuitName '.rcx'], rateTag);

            % trigger duration (fixed)
            obj.RP.SetTagVal('triggerDuration', ...
                             1000 * obj.triggerDuration);
            % button hold time (fixed)
            obj.RP.SetTagVal('buttonHoldDuration', ...
                             1000 * obj.buttonHoldDuration);

            if useRZ6
                % for the RZ6, use the full scale of 10V, and bring it back
                % into range using the programmable attenuator. The
                % attenuation is in dB, and as per documentation,
                % attenuation will be achieved using combination of digital
                % and analog scaling depending on the magnitude
                obj.RP.SetTagVal('chan1Scaler', 10);
                obj.RP.SetTagVal('chan2Scaler', 10);
                obj.RP.SetTagVal('chan1Att', 20*log10(10 / obj.channel1Scale));
                obj.RP.SetTagVal('chan2Att', 20*log10(10 / obj.channel2Scale));
            else
                % for the RP2.1, just use the D/A output voltage, no
                % modifications
                obj.RP.SetTagVal('chan1Scaler', obj.channel1Scale);
                obj.RP.SetTagVal('chan2Scaler', obj.channel2Scale);
            end
            
            % zero tag the relevant buffers
            obj.RP.ZeroTag('audioChannel1');
            obj.RP.ZeroTag('audioChannel2');
            obj.RP.ZeroTag('triggerIdx');
            obj.RP.ZeroTag('triggerVals');
            obj.RP.ZeroTag('buttonPressValue');
            obj.RP.ZeroTag('buttonPressSample');

            % now attempt to actually run the circuit
            obj.RP.Run;
            if obj.RP.GetStatus ~= 7
                obj.RP.close();
                error('TDT connection error. Try rebooting the TDT.');
            end
            
            % display some status information to the user
            fprintf('Channel 1, float [-1.0, 1.0] --> [-%2.4f, %2.4f] V\n', ...
                 obj.channel1Scale, obj.channel1Scale);
            fprintf('Channel 2, float [-1.0, 1.0] --> [-%2.4f, %2.4f] V\n', ...
                 obj.channel2Scale, obj.channel2Scale);
             
            % get and store sample rate once the circut is loaded and
            % running
            obj.sampleRate = obj.RP.GetSFreq();
            
            % validate the xorVal
            buttonBoxOK = obj.configure_button_box();
            if ~buttonBoxOK
               btState = warning('backtrace');
               warning('backtrace', 'off');
               warning('Button box might be incorrectly configured.') 
               warning('backtrace', btState.state);
            end

            obj.status = sprintf('No stimulus loaded.');
            if useRZ6
                obj.deviceType = 'RZ6';
            else
                obj.deviceType = 'RP2.1';
            end
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
        % should be specified in seconds. Trigger values and index values
        % should be non-negative.
        %
        % last updated: 2016-02-08, LV, lennyv_at_bu_dot_edu
            
            %%%%%%%%%%%%%%%%%%%%
            % input validation %
            %%%%%%%%%%%%%%%%%%%%

            if nargin < 3
                triggerInfo = [];
            else
                if iscolumn(triggerInfo)
                    triggerInfo = triggerInfo';
                end
                if size(triggerInfo, 2) == 1
                    triggerInfo(:, 2) = 1:size(triggerInfo, 1);
                    warning(['triggerInfo should be specified as ',...
                             '[idx, val], array, and the values '...
                             'should all be positive...only one column '...
                             'detected, so using those as event index and '...
                             'using sequential event values'])
                end
                
                if (size(triggerInfo, 2) ~= 2) || ...
                    (length(size(triggerInfo)) ~= 2) || ...
                    any(triggerInfo(:) < 0)
                    error(['triggerInfo must be specified as ',...
                           '[idx, val], array, and the values '...
                           'should all be positive.'])
                end
            end
            
            if any(abs(audioData) > 1)
                error('All audio data must be scaled between -1.0 and 1.0.')
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
            
            if (strcmpi(obj.paradigmType, 'playback_1channel') || ...
               strcmpi(obj.paradigmType, 'playback_2channel'))
                % convert down to single precision floating point, 
                % since that's what the TDT natively uses for DAC
                if ~isa(audioData, 'single')
                    audioData = single(audioData);
                end
            else
                % otherwise scale up by 2**15 for integer transfer
                audioData = audioData .* (2^15);
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
            
            if size(triggerInfo, 1) > 2290
                error(['Circuit can only support a maximum of 2290 ', ...
                       'trigger values. Reduce the number of ', ...
                       'triggers specified in triggerInfo.'])
            end

            %%%%%%%%%%%%%%%%%%%%%
            % write data to TDT %
            %%%%%%%%%%%%%%%%%%%%%
           
            % hack - the WriteTagVEX methods don't like single value inputs
            % also correct for 1 sample difference in index but add 1 back to
            % account for one sample zero padding at beginning...so +1 -1
            % cancel out
            triggerIdx = [triggerIdx; - 1];
            triggerVals = [triggerVals; 0];
            
            % reset buffer indexing and zeroTag everything
            obj.reset_buffers(true)
            % size +2/+1 are intentional on next lines
            %obj.RP.SetTagVal('stimSize', size(audioData, 1)+2);
            obj.stimSize = size(audioData, 1) + 1;

            % note: 0 padding below appears to eliminate the clicking noise
            % when buffers are written to or accessed - LV 2016-02-08
            if ~strcmpi(obj.paradigmType, 'playback_2channel_16bit')
                fprintf('Writing to channel 1 buffer...\n')
                curStatus = obj.RP.WriteTagVEX('audioChannel1', 0, 'F32',...
                                                [0; audioData(:, 1)]);
                if ~curStatus
                    error('Error writing to audioChannel1 buffer.')
                end

                if obj.nChans == 2
                    fprintf('Writing to channel 2 buffer...\n')
                    curStatus = obj.RP.WriteTagVEX('audioChannel2', 0, ...
                                                   'F32', ...
                                                   [0; audioData(:, 2)]);
                    if ~curStatus
                        error('Error writing to audioChannel2 buffer.')
                    end
                end
            else
                fprintf('Writing 2 channels to audio buffer...\n')
                curStatus = obj.RP.WriteTagVEX('audioChannel1', 0, 'I16',...
                                                [[0; 0], audioData']);
                if ~curStatus
                    error('Error writing to audioChannel1 buffer.')
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
        % last updated: 2016-04-21, LV, lennyv_at_bu_dot_edu

            if obj.stimSize == 0
                error('No stimulus loaded.')
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
        % last updated: 2015-04-03, LV, lennyv_at_bu_dot_edu
            
            stat = obj.RP.SetTagVal('stopSample', 0);
            if ~stat
                error('Error setting stop sample.')
            end
            pause(0.02);
            currentSample = obj.get_current_sample();
            obj.status = sprintf('stopped at buffer index %d', currentSample);
        end
       
        
        function play_blocking(obj, stopAfter, debugMode)
        % tdt.play_blocking(stopAfter, debugMode)
        %
        % Plays the contents of the audio buffers on the TDT and holds up
        % Matlab execution while doing so.
        %
        % Inputs:
        % --------------------------------------------------------------------
        % stopAfter - the sample number at which playback should cease. If not
        % specified, playback will continue until the end of the stimulus is 
        % reached. Note: buffer positions will automatically reset to 0
        % after the end fo the stimulus is reached.
        %
        % debugMode - if set, display the buffer sample numbers on screen
        % while Matlab is being held up
        %
        % version added: 1.1
        % last updated: 2016-04-22, LV, lennyv_at_bu_dot_edu

            if obj.stimSize == 0
                error('No stimulus loaded.')
            end

            if nargin < 2 || isempty(stopAfter)
                stopAfter = [];
                stopAfterSample = obj.stimSize;
            else
                stopAfterSample = stopAfter;
            end
            
            if nargin < 3 || isempty(debugMode)
                debugMode = false;
            end
            
            if stopAfterSample < obj.get_current_sample()
                error(['Buffer index already passed desired stop point. ' ...
                       'Did you mean to rewind the buffer first?'])
            end
            
            stat = obj.RP.SetTagVal('stopSample', stopAfterSample);
            if ~stat
                error('Error setting stop sample.')
            end
            
            currentSample = obj.get_current_sample();
            fprintf('Playing stimulus in blocking mode...\n');
            fprintf('    Started playing at sample %d\n', currentSample);
            obj.RP.SoftTrg(1);
            pause(0.01)
            currentSample = obj.RP.GetTagVal('chan1BufIdx');
            try
                if debugMode
                    fprintf(1, '    Current sample: %08d\n', ....
                            currentSample);
                end
                while currentSample > 0 && currentSample <= stopAfterSample
                    currentSample = obj.RP.GetTagVal('chan1BufIdx');
                    if debugMode && currentSample
                        fprintf(1, '\b\b\b\b\b\b\b\b\b%08d\n', ...
                                currentSample);
                    end
                    pause(0.1);
                end
                fprintf(1, '    Sample %d reached\n', stopAfterSample);
            catch ME
                obj.pause()
                fprintf(['\n' obj.status]);
                throw(ME);
            end
            currentSample = obj.get_current_sample();
            obj.status = sprintf('stopped at buffer index %d', currentSample);
            fprintf('...done.\n')
            if currentSample >= obj.stimSize
                fprintf(1, ['All samples played; ', ...
                    'rewind or reset the buffers.\n']);
            end
        end

        
        function rewind(obj)
        % tdt.rewind()
        %
        % Rewinds the audio buffer without clearing it, and clears the button
        % pres buffers. Useful when new audio data does not need to be loaded 
        % into the TDT (to minimize the communication time between PC and TDT).
        %
        % last updated: 2015-04-21, LV, lennyv_at_bu_dot_edu
        
            obj.reset_buffers(false);
            currentSample = obj.get_current_sample();
            obj.get_button_presses;
            obj.status = sprintf('stopped at buffer index %d', currentSample);
        end
       
        
        function reset(obj)
        % tdt.reset()
        %
        % Rewinds the buffer and sets all values in the buffers to 0. 
        %
        % last updated: 2015-03-11, LV, lennyv_at_bu_dot_edu

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
        % last updated: 2015-03-11, LV, lennyv_at_bu_dot_edu
        
            statusVal = obj.RP.SetTagVal('arbitraryEvent', eventVal);
            if ~statusVal
                error('Event could not be written.')
            end
            pause(0.01);
            obj.RP.SoftTrg(4);
        end
        
        
        function [pressVals, pressSamples] = get_button_presses(obj)
           % [pressVals, pressSamples] = tdt.get_button_presses()
           %
           % retrieves button press values and the sample that each press
           % occurred relative to the start of playback
           %
           % note: will only store 2000 button presses without requiring the
           % reset() function to be called due to the size limitation on the 
           % button information storage buffers. Until reset() or rewind() is 
           % called, subsequent button presses will not be recorded.
           % Whenever rewind() is called, however, the button presses
           % recorded up until that point will be stored in the TDT object.
           % 
           % Note that pressVals are expressed as 2^(realButtonNumber-1), i.e., 
           % button1 = 1, button2 = 2, button3 = 4, button4 = 8.
           %
           % pressSamples output is expressed in *samples* relative to the
           % stimulus. There is a constant 3-4 sample delay when using the
           % "WordIn" component of the TDT, which for practical purposes makes
           % little difference in the accuracy of reaction times calculated
           % from pressSamples.
           %
           % Will return NaN for both outputs if no buttons have been pressed.
           %
           % version added: 1.3
           % last updated: 2016-04-22 LV, lennyv_at_bu_dot_edu

           nPress = obj.RP.GetTagVal('nPress');

           if nPress > 0
               pressVals = obj.RP.ReadTagVEX('buttonPressVals',0, ...
                                             nPress, 'I32', 'F64', 1);
               % convert back to Matlab-style indexing here with +1...but
               % +1 may not be necessary anymore because of leading 0 that was
               % appened to the stimulus? 
               pressSamples = obj.RP.ReadTagVEX('buttonPressSamples', 0,...
                                                nPress, 'I32', 'F64', 1);
           else
               pressVals = NaN;
               pressSamples = NaN;
           end

        end
        
        
        function [currentSample1, trigBufSample1] = ...
                get_current_sample(obj, checks)
        % [audioIdx, triggerIdx] = tdt.get_current_sample(checks)
        %
        % Gets the current buffer position for the audio stimuli (output 1) and
        % for triggers (output 2).
        %
        % An error is raised if the audio buffers or the trigger buffers become
        % misaligned.
        %
        % last updated: 2016-04-21, LV, lennyv_at_bu_dot_edu

            if nargin == 1
                checks = true;
            end
        
            currentSample1= obj.RP.GetTagVal('chan1BufIdx');
            if checks && strcmpi(obj.paradigmType, 'playback_2channel')
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
        function buttonBoxOK = configure_button_box(obj)    
        % try auto-setting the button box xor value; if incorrectly
        % configured, the button box buffer will be receiving meaningless
        % input
        % version added: 1.6
        % last modified: 2016-04-21 LV, lennyv_at_bu_dot_edu
            buttonBoxOK = false;
            if isempty(obj.xorVal)
                buttonPresses = obj.get_button_presses();
                if length(buttonPresses) == 1 && isnan(buttonPresses(1))
                    fprintf(1, 'xorVal autodetected successfully! (0)\n');
                    buttonBoxOK = true;
                elseif all(buttonPresses == buttonPresses(1))
                    tryVal = buttonPresses(1);
                    obj.RP.SetTagVal('xorVal', tryVal);
                    obj.reset();
                    pause(0.5);
                    buttonPresses = obj.get_button_presses();
                    if length(buttonPresses) == 1 && isnan(buttonPresses(1))
                        fprintf(1, 'xorVal autodetected successfully! ');
                        fprintf(1, '(%d)\n', tryVal);
                        buttonBoxOK = true;
                        obj.xorVal = tryVal;
                    else
                        fprintf(1, 'WARNING: ');
                        fprintf(1, 'xorVal was set to %d, ', tryVal);
                        fprintf(1, 'but that does not seem to be working.\n');
                        fprintf(1, 'Button box may not function correctly.\n');
                    end
                else
                    fprintf(1, 'WARNING: ');
                    fprintf(1, 'Could not auto-configure xorVal.\n');
                    fprintf(1, 'Button box may not function correctly.\n');
                end
            else
                obj.RP.SetTagVal('xorVal', obj.xorVal);
                pause(0.01);
                obj.reset()
                pause(0.5);
                buttonPresses = obj.get_button_presses();
                if length(buttonPresses) ~= 1 || ~isnan(buttonPresses(1));
                    fprintf(1, 'WARNING: ');
                    fprintf(1, ['xorVal is likely incorrect; '...
                                'try setting it to [] to autoconfigure.\n']);
                    fprintf(1, 'Button box might not function correctly.\n');
                else
                    buttonBoxOK = true;
                end
            end
        end
        
        
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
        % 0. Note: will always reset the button-press buffers.
        %
        % last updated: 2015-10-25, LV, lennyv_at_bu_dot_edu           
            
            obj.RP.SoftTrg(2);
            pause(0.01);
            
            if clearBuffer
                obj.RP.ZeroTag('audioChannel1');
                obj.RP.ZeroTag('audioChannel2');
                obj.RP.ZeroTag('triggerIdx');
                obj.RP.ZeroTag('triggerVals');
                obj.RP.SetTagVal('stopSample', 0);
                obj.stimSize = 0;
            end
            
            % always clear button press buffers
            obj.RP.ZeroTag('buttonPressValue');
            obj.RP.ZeroTag('buttonPressSample');
            
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
        % last updated: 2015-03-11, LV, lennyv_at_bu_dot_edu

            obj.reset_buffers(true);
            obj.RP.Halt;
            pause(0.01);
            obj.RP.ClearCOF;
            close(obj.f1);
            close(obj.hiddenFigure);
            obj.status = sprintf('Not connected.');
        end
        
        
        function rateTag = parse_inputs(obj, paradigmType, ...
                                        requestedSampleRate, scaling, varargin)
        % helper function to help setup the constructor
        % version added: 1.6
        % last modified: 2016-04-21 LV, lennyv_at_bu_dot_edu
        % support old and new style parsing and defaults
            if any(cellfun(@ischar, varargin)) || isempty(varargin) % kwargs
                p = inputParser();
                addRequired(p, 'paradigmType', @ischar);
                addRequired(p, 'requestedSampleRate', @isnumeric);
                addRequired(p, 'scaling', @isnumeric);
                addOptional(p, 'triggerDuration', 0.005, @isnumeric);
                addOptional(p, 'buttonHoldDuration', 0.2, @isnumeric);
                addOptional(p, 'xorVal', [], @isnumeric);
                addOptional(p, 'figNum', 9999, @isnumeric);
                
                p.parse(paradigmType, requestedSampleRate, scaling, varargin{:})
    
                inputs = p.Results;
            else % all numeric, like old style
                btState = warning('backtrace');
                warning('backtrace', 'off');
                warning('Deprecation warning!');
                warning('Use keyword argument pairs in future releases (see help).') 
                warning('backtrace', btState.state);

                inputs.paradigmType = paradigmType;
                inputs.requestedSampleRate = requestedSampleRate;
                inputs.scaling = scaling;
                if ~isempty(varargin)
                    inputs.triggerDuration = varargin{1};
                else
                    inputs.triggerDuration = [];
                end
                if isempty(inputs.triggerDuration)
                    inputs.triggerDuration = [];
                end

                if length(varargin) > 1
                    inputs.buttonHoldDuration = varargin{2};
                else
                    inputs.buttonHoldDuration = [];
                end
                if isempty(inputs.buttonHoldDuration)
                    inputs.buttonHoldDuration = 0.2;
                end
                
                if length(varargin) > 3
                    inputs.xorVal = varargin{3};
                else
                    inputs.xorVal = [];
                end
                if isempty(inputs.xorVal)
                    inputs.xorVal = [];
                end
                
                if length(varargin) > 4
                    inputs.figNum = varargin{4};
                else
                    inputs.figNum = [];
                end
                if isempty(inputs.figNum)
                    inputs.figNum = 9999; 
                end
            end
            
            %% input checks
            
            % by default, set scaling equal on both channels (2 channel), or
            % use monaural playback
            if length(inputs.scaling) < 2
                if strcmpi(inputs.paradigmType, 'playback_2channel')
                    inputs.scaling(2) = inputs.scaling(1);
                else
                    inputs.scaling(2) = 0;
                end
            end
            
            %%% voltage scaling
            if any(abs(inputs.scaling) > 10)
                error('Scaling must be specified and be within +/-10 V.')
            end
            
            if 97 == inputs.requestedSampleRate
                rateTag = 4;
            elseif 48 == inputs.requestedSampleRate 
                rateTag = 3;
            elseif 24 == inputs.requestedSampleRate
                rateTag = 2;
            elseif 12 == inputs.requestedSampleRate 
                rateTag = 1;
            else
                error('invalid sample rate specified (must be 97, 48, 24, 12)')
            end
            
            if inputs.buttonHoldDuration <= 0
                error('buttonHoldDuration must be positive')
            end
            
            if inputs.triggerDuration <= 0
                error('triggerDuration must be positive')
            end
            
            % scaling factors
            obj.channel1Scale = single(inputs.scaling(1));
            obj.channel2Scale = single(inputs.scaling(2));
            
            obj.paradigmType = inputs.paradigmType;
            obj.buttonHoldDuration = inputs.buttonHoldDuration;
            obj.triggerDuration = inputs.triggerDuration;
            obj.xorVal = inputs.xorVal;
            
            % Start ActiveX controls and hides the figure at the start of
            % each block
            obj.f1 = figure(inputs.figNum);
            set(obj.f1,'Position', [5 5 30 30], 'Visible', 'off');
            obj.RP = actxcontrol('RPco.x', [5 5 30 30], obj.f1);
            % open up a new figure and hide it, so that the first plot command
            % doesn't screw things up
            obj.hiddenFigure = figure(inputs.figNum + 1);
            set(gcf, 'Visible', 'off');

        end
    end
end
