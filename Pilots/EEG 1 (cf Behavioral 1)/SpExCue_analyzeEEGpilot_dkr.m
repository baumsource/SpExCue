 % This script includes methods to plot multiple curves from a given
% experimental session, i.e. from the attendBack and/or attendFig sessions.
%
%  NOTE: this analysis script requires that the raw BDF file has already
%  been processed AND that specific filename/filepath conventions were
%  followed.

%**************************************************************************
% Preliminary stuff to ensure that the analysis is setup to be executed.
%**************************************************************************
% Ensure that necessary files from 'eeglab' are available
if ~exist('pop_loadset', 'file')
    eeglab
end
clear; clc; close force all;


% Primary path definitions required for this analysis
dirBDF = '/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Experiments/Pilot/Results/RB/';  % directory with the .bdf/.set files
dirLOCS = '/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Tools/EEG/';       % directory with the eeglab locations file
dirSupp = '/Users/rbaumgartner/Documents/ARI/Projects/SpExCue/Experiments/MATLAB_general/';  % directory with required analysis functions, e.g. 'eeglabTrig2dkrTrigID.m' 



%**************************************************************************
% User specified parameters for analysis
%**************************************************************************
% subjectID = {'001', '002'};
subjectID = {'RB'};

% Determine which .set file to be used for the analysis...
%   if '', run on data without ICA components removed
%   if '_ICAclean', run on data with ICA components removed 
setFileLabel = '_ICAraw'; 
% Other parameters that influence the proper selection of the .set file
fs = 100;  % ASSUME the sampling rate... note that there is a subsequent check to confirm that this value is valid
filtFreqLbl = '25';  % STRING used to identify the filename from the 'preprocessed scripts'. This value is the upper frequency of the bandpass filter in my preprocessing scripts.


% If multiple curves are desired per experimental session, use a cell array
% for 'subsetConditions_' such that each cell entry consists of a double
% array with a list of experimental conditions to be combined into a single
% curve. ALSO, multiple curves will require that 'subsetResponses_' is a
% cell array containing appropriate strings for each curve to be generated.
% If only a single curve is desired, 'subsetConditions_' can simply be a
% double array and 'subsetResponses_' ans simply be a single string.
%
% NOTE: for the "difference plots" the order of the 'subsetConditions'
% arrays are CRUCIAL.  They must be specified as:
%                             {[44 54], [42 52], [41 51], [40 50]} 
% for appropriate plots.


% % --------
subsetResponses_detectFig = {'all', 'all'};  % VALID INPUTS: 'all', 'correct', 'wrong'
subsetConditions_detectFig = {[45], [25]}; % 10 (fig. std.), 20 (back. std.), 30 (fig. deviant), 40 (back. deviant)
subsetConditionsString_detectFig = {'BackDEV', 'BackSTD'};  % Used only for the 'laymanPlotLabelsFLAG' plots

subsetResponses_detectBack = {'all', 'all'};  % VALID INPUTS: 'all', 'correct', 'wrong'
subsetConditions_detectBack = {[35], [15]}; % 10 (fig. std.), 20 (back. std.), 30 (fig. deviant), 40 (back. deviant)
subsetConditionsString_detectBack = {'FigDEV', 'FigSTD'};  % Used only for the 'laymanPlotLabelsFLAG' plots




nChans = 33;  % Total number of EEG channels in the data set
nChansLocsFilename = 'biosemi_eloc.locs';  % filename of the eeglab locations file
% nChans = 32;  % Total number of EEG channels in the data set
% nChansLocsFilename = 'biosemi_eloc.locs';  % filename of the eeglab locations file

epochStart = -0.5;  % seconds before the trigger event
epochEnd = 1; % seconds after the trigger event
baselineCorrectInterval = [-100 0];

threshEpoch = 70; %(Inyong used 70)  % magnitude in microvolts for which epochs should be excluded.  IFF empty, don't threshold.

alphaLim = [8 12];  % For plot where alpha frequencies are removed... range ([low high]) of frequencies in Hz to band-reject.

savePlotsFLAG = 0;
laymanPlotLabelsFLAG = 1;  % specify if plots should have simple legend labels for the general audience
plotIndividualDataFLAG = 1;  % specify if plots for each individual subject should be generated
plotPopulationDataFLAG = 1;  % specify if plots for the *population* grand average should be generated
chanNum = [32];%[5,8,23,26,32];  % channel number(s) to be used in the plots

% Settings for for time domain auditory ERPs
plotTimeDomainFLAG = 1;  % plots of raw epochs
plotTimeDomainAlphaModFLAG = 0; % plots of epochs where alpha frequencies are removed
plotTimeDomainDiffPlotsFLAG = 1; % difference plots, i.e. switch minus no-switch ERPs
gridFreq = 1/2;  % This is presentation rate (in Hz) for the chords of the osullivant stimulus.  
ylimInd = 7;  % Maximum range (+/-) for the plots for INDIVIDUAL data
ylimPop = 4;  % Maximum range (+/-) for the plots for POPULATION data

% Settings for scalp plots for *difference* waveform 
plotTopoDiffFLAG = 0;  % specify if the topoplots should be generated
topoPlotTimes = [-50 150 400];  % Times in MILLISECONDS of the epoch to generate the topoplot 

% Settings for Event-Related Spectral Pertubation (ERSP) plots
plotERSPplotsFLAG = 0;  % specify if the topoplots should be generated
cycles = [1 0.5]; % The number of cycles in the window size
freqRng = [2 55];
padratio = 8;
winsize = 256;
alpha = nan; %0.15; % If NaN, don't plot with computed significance
timesout = 200;  % default 200
plotFreqRng = [2 20];

% Settings for ERP Image plot
plotERPimageFLAG = 0;  % specify if the ERP image plot should be generated
smooth = 10;   % Smoothing parameter (number of trials). {Default: 5} -- erpimage() equivalent: 'avewidth' 
decimate = 1;  % Decimate parameter (i.e. ratio of trials_in/trials_out) -- erpaimge() equivalent: 'decimate' {Default: 0}



%---------------------------------------------------
% Plot settings
%---------------------------------------------------
imageFileType = 'PDF'; % '-depsc'; '-dpng'; 'PDF'; % specify what type of file the images should be saved as
imgWidth = 10; % Parameters for PDF images
imgHeight = 6; % Parameters for PDF images

% An array of RGB values for the different curves of a given session
col = [...
    0 0 1;...  % blue
    1 0 0;...  % red
    0 1 0;...  % green
    1 0 1;...  % magenta
    0 1 1;...  % cyan  
    1 0.8398 0;...          % gold
    0.6016 0.8008 0.1953;... %YellowGreen
    0.6445 0.1641 0.1641;...% brown
    0.5781 0 0.8242;... % darkViolet
    1 0.6250 0.4766;... % lightSalmon
    ];
% col=hsv(10);  % Alternate color scheme

ftSz_axis = 10;  % Specify font size across plots;
ftSz_lgnd = 6;
lineWid = 1.5; % Specify line widths across plots
lineWidGrid = 0.5;  % Line width of the grid (if used)
lgndLocation = 'NorthWest';

anLabels = {'Fp1 (Ch1)','AF3 (Ch2)', 'F7 (Ch3)', 'F3 (Ch4)','FC1 (Ch5)', ... 
    'FC5 (Ch6)', 'T7 (Ch7)', 'C3 (Ch8)', 'CP1 (Ch9)','CP5 (Ch10)', ...
    'P7 (Ch11)', 'P3 (Ch12)', 'Pz (Ch13)', 'PO3 (Ch14)','O1 (Ch15)', ...
    'Oz (Ch16)', 'O2 (Ch17)', 'PO4 (Ch18)', 'P4 (Ch19)', 'P8 (Ch20)',...
    'CP6 (Ch21)', 'CP2 (Ch22)', 'C4 (Ch23)', 'T8 (Ch24)','FC6 (Ch25)',...
    'FC2 (Ch26)',  'F4 (Ch27)', 'F8 (Ch28)', 'AF4 (Ch29)', 'Fp2 (Ch30)',...
    'Fz (Ch31)',  'Cz (Ch32)','VEOG (Ch33)'};
     

% chanLabels = {'Fp1 (Ch1)', 'AF7 (Ch2)', 'AF3 (Ch3)', 'F1 (Ch4)', 'F3 (Ch5)',... 
%     'F5 (Ch6)', 'F7 (Ch7)', 'FT7 (Ch8)', 'FC5 (Ch9)', 'FC3 (Ch10)', 'FC1 (Ch11)',...
%     'C1 (Ch12)', 'C3 (Ch13)', 'C5 (Ch14)', 'T7 (Ch15)', 'TP7 (Ch16)',...
%     'CP5 (Ch17)', 'CP3 (Ch18)', 'CP1 (Ch19)', 'P1 (Ch20)', 'P3 (Ch21)',...
%     'P5 (Ch22)', 'P7 (Ch23)', 'P9 (Ch24)', 'PO7 (Ch25)', 'PO3 (Ch26)',...
%     'O1 (Ch27)', 'Iz (Ch28)', 'Oz (Ch29)', 'POz (Ch30)', 'Pz (Ch31)', 'CPz (Ch32)',...
%     'Fpz (Ch33)', 'Fp2 (Ch34)', 'AF8 (Ch35)', 'AF4 (Ch36)', 'AFz (Ch37)',...
%     'Fz (Ch38)', 'F2 (Ch39)', 'F4 (Ch40)', 'F6 (Ch41)', 'F8 (Ch42)', 'FT8 (Ch43)',...
%     'FC6 (Ch44)', 'FC4 (Ch45)', 'FC2 (Ch46)', 'FCz (Ch47)', 'Cz (Ch48)',...
%     'C2 (Ch49)', 'C4 (Ch50)', 'C6 (Ch51)', 'T8 (Ch52)', 'TP8 (Ch53)',...
%     'CP6 (Ch54)', 'CP4 (Ch55)', 'CP2 (Ch56)', 'P2 (Ch57)', 'P4 (Ch58)',...
%     'P6 (Ch59)', 'P8 (Ch60)', 'P10 (Ch61)', 'PO8 (Ch62)', 'PO4 (Ch63)', 'O2 (Ch64)',... 
%     'VEOG (Ch65)'};


%--------------------------------------------------------------------------
%  Calculate useful values and run some checks on the input parameters
%--------------------------------------------------------------------------

% Force the input parameters to be of type CELL, regardless of how the user
% specified the parameters
if ~iscell(subsetResponses_detectFig)
    subsetResponses_detectFig = {subsetResponses_detectFig};
    subsetConditions_detectFig = {subsetConditions_detectFig};
    subsetResponses_detectBack = {subsetResponses_detectBack};
    subsetConditions_detectBack = {subsetConditions_detectBack};
end

% Calculate how many curves per experimental session should be generated.
nCurves = length(subsetResponses_detectFig); 

% Calculate the number of samples based on the duration of the epoch
nSamps = length(epochStart:1/fs:(epochEnd-1/fs));

% CHECK appropriate dimensions for the 'Responses_' and 'Conditions_'
% input parameters.
tmp = [...
    length(subsetResponses_detectFig),...
    length(subsetConditions_detectFig),...
    length(subsetResponses_detectBack),...
    length(subsetConditions_detectBack)];
    
if length(unique(tmp)) > 1
    error('DKR: dimensions for the ''Responses_'' and ''Conditions_'' user specified parameters MUST AGREE.')
end

% Compute the indices for the onset of each chord in the osullivan stimului
% to be used in the plots.
if gridFreq ~= 0
    gridIndices = (-100:1:100)*(1/gridFreq);
    gridIndices = gridIndices(gridIndices >= epochStart & gridIndices <= epochEnd);
end

% Need this path to access additional analysis functions, e.g.
% 'eeglabTrig2dkrTrigID.m' 
addpath(dirSupp) 


%--------------------------------------------------------------------------
% INITIALIZE output arrays 
%--------------------------------------------------------------------------
% Arrays to store the *average* epoch for each subjectID, i.e. the grand
% average epoch for each subject is saved.
populationMeanDATA_M1 = zeros(nChans, nSamps, length(subjectID), nCurves);  
populationMeanDATA_M0 = zeros(nChans, nSamps, length(subjectID), nCurves);
% populationMeanDATAnoAlpha_fig = zeros(nChans, nSamps, length(subjectID), nCurves);   
% populationMeanDATAnoAlpha_back = zeros(nChans, nSamps, length(subjectID), nCurves);

% Arrays were *individual* epochs are concatenated for all subjectID
populationDATA_M1 = cell(1,nCurves); 
populationDATA_M0 = cell(1,nCurves);
% populationDATAnoAlpha_fig = cell(1,nCurves);
% populationDATAnoAlpha_back = cell(1,nCurves);

% Arrays to provide an ID for the concatenated individual epochs
populationSUBJID_M1 = cell(1,nCurves); 
populationSUBJID_M0 = cell(1,nCurves);

% % Arrays to store the ERSP
% popERSPfig_ersp = cell(length(subjectID), nCurves);
% popERSPfig_itc = cell(length(subjectID), nCurves);
% popERSPfig_powbase = cell(length(subjectID), nCurves);
% popERSPfig_times = cell(length(subjectID), nCurves);
% popERSPfig_freqs = cell(length(subjectID), nCurves);
% 
% popERSPback_ersp = cell(length(subjectID), nCurves);
% popERSPback_itc = cell(length(subjectID), nCurves);
% popERSPback_powbase = cell(length(subjectID), nCurves);
% popERSPback_times = cell(length(subjectID), nCurves);
% popERSPback_freqs = cell(length(subjectID), nCurves);


%**************************************************************************
% Loop through the list of 'subjectID' values to be analyzed. 
%**************************************************************************
hwait = waitbar(0,'Please wait...'); ctr = 1;
for n = 1:length(subjectID)
    
    %----------------------------------------------------------------------
    % load existing datasets saved from file
    %----------------------------------------------------------------------
    % Construct the filename and filepath for a given subject's data
    filename = ['EEGpilot_subj',subjectID{n},'_resamp100_filt25',setFileLabel,'.set'];%['exp2subj', , '_detectFig_resamp', num2str(fs), 'filt', filtFreqLbl, setFileLabel, '.set']; 
    filepath = dirBDF;
    
    % Load the actual data
    EEG = pop_loadset('filename', filename, 'filepath', filepath); 
    EEG = eeg_checkset(EEG);  % verify consistency of the fields of an EEG dataset
    EEG = pop_chanedit(EEG, 'load',{[dirLOCS, nChansLocsFilename] 'filetype' 'autodetect'});  % Set channel locations

    clear filename filepath

    % HACK check to confirm that the sampling rate used to calculate the
    % value for nSamps is consistent with the actual data.
    if fs ~= EEG.srate
        error('DKR: Initial sampling rate is not consistent with actual sampling rate of the data')
    end
    
    %----------------------------------------------------------------------
    % Convert trigger values to decimal values using ONLY the lower 8-bits
    % of binary information.
    %----------------------------------------------------------------------
    eventList = eeglabTrig2dkrTrigID([EEG.event.type]); 
    for i = 1:length(eventList)
        EEG.event(i).type = eventList(i);
    end
    
    % Archive the EEG data prior to further manipulation of the trigger
    % values for subjectID{n}
    EEGRAW = EEG;
    
    % Initialize arrays to save the number of epochs used for each curve
    % for each individual.
    numIndEpochs = nan(nCurves,1);
    
    
    %----------------------------------------------------------------------
    % Create masks indicating triggers corresponding to 'correct' (203) or
    % 'wrong' (204) trials for both the detectFig and detectBack sessions.
    %   NOTE: Must shift mask values by 1 and 2 bins to get the indices for
    %   the trigger corresponding to the actual deviant, not just to the
    %   correct/wrong response (or to the end of the sound stimulus).   
    %----------------------------------------------------------------------
    modValTargM1 = round(mod(eventList,100)/10) == 30; % The deviant has a 5 in the ones position for the trigger integers
    modValTargM0 = round(mod(eventList,100)/10) == 10;  % The deviant has a 5 in the ones position for the trigger integers

%     tmp = eventListFig == 203;          % for 'detectFIG' session
%     tmp1 = circshift(tmp,-1);
%     tmp2 = circshift(tmp,-2);
%     sv1 = tmp1 & modValTargFig;
%     sv2 = tmp2 & modValTargFig;
%     corrIndicesFig = sv1 | sv2;
%     
%     tmp = eventListFig == 204;          
%     tmp1 = circshift(tmp,-1);
%     tmp2 = circshift(tmp,-2);
%     sv1 = tmp1 & modValTargFig;
%     sv2 = tmp2 & modValTargFig;
%     wrongIndicesFig = sv1 | sv2;
    

    %**********************************************************************
    % Loop through the data for both sessions depending on how many curves
    % were requested per the input parameters.
    %**********************************************************************
    for j = 1:nCurves
        
        % Select the data for the user specified curve to be generated.
%         subsetResponses_detectFigTMP = subsetResponses_detectFig{j};
%         subsetConditions_detectFigTMP = subsetConditions_detectFig{j};
%         subsetResponses_detectBackTMP = subsetResponses_detectBack{j};
%         subsetConditions_detectBackTMP = subsetConditions_detectBack{j};
        
        
        %******************************************************************
        % Create an array of trigger indicies that should be used for
        % further processing.
        %-----------------------------------------------------
        % User specifies subset of responses: 'all', 'correct', 'wrong'
        % User specifies subset of conditions: an array of integers
        %                                corresponding to test condition ID
        % Do this separately for 'detectFig' and 'detectBack' so that
        % target and non-target conditions can be compared between
        % experimental sessions.
        %******************************************************************

        % For 'detectFIG'... 
%         switch lower(subsetResponses_detectFigTMP)  
%             case 'all'
%                 respIndicesFig = corrIndicesFig + wrongIndicesFig;
%             case 'correct'
%                 respIndicesFig = corrIndicesFig;
%             case 'wrong'
%                 respIndicesFig = wrongIndicesFig;
%             otherwise
%                 error('DKR: Invalid string for ''subsetResponses_detectFig''');
%         end
%         
%         condIndicesFig = zeros(length(eventListFig), 1);
%         for i = 1:length(subsetConditions_detectFigTMP)
%             tmp = eventListFig == subsetConditions_detectFigTMP(i);
%             condIndicesFig = condIndicesFig + tmp;
%         end


        %------------------------------------------------------------------
        % Logical comparison to mark only the desired indices based on the
        % user specified subset parameters: 'Responses_' and 'Conditions_' 
        %------------------------------------------------------------------
%         desiredIndicesFig = respIndicesFig & condIndicesFig;
%         desiredIndicesBack = respIndicesBack & condIndicesBack;
    desiredIndicesM1 = modValTargM1;
    desiredIndicesM0 = modValTargM0;
%         % FIX BECAUSE EVEN A 'NO RESPONSE' CAN BE CORRECT
%         desiredIndicesFig = condIndicesFig;
%         desiredIndicesBack = condIndicesBack;

        %------------------------------------------------------------------
        %  Re-label the trigger values of the "desired" trials, i.e. trials
        %  for a certain test condition.
        %------------------------------------------------------------------
        tmp = find(desiredIndicesM1 ~= 0);  % 'tmp' should be ONLY indicies of the desired triggers
        for i = 1:length(tmp)
            EEG.event(tmp(i)).type = 3;  % reset desired trigger to 1
        end

        tmp = find(desiredIndicesM0 ~= 0); % 'tmp' should be ONLY indicies of the desired triggers
        for i = 1:length(tmp)
            EEG.event(tmp(i)).type = 1;  % reset desired trigger to 1
        end

        % Define CELL ARRAY of trigger value to be epoched with 'pop_epoch'
        desiredTrigsM1_cell = {3};
        desiredTrigsM0_cell = {1};
        

        %------------------------------------------------------------------
        % Epoch the data for each condition...  handle the selection of
        % triggers associated with a given test condition
        %------------------------------------------------------------------
        tmpEEGM1 = pop_epoch(EEG, desiredTrigsM1_cell, [epochStart epochEnd], 'newname', ['subsetTrigs'], 'epochinfo', 'yes');
        tmpEEGM1 = pop_rmbase(tmpEEGM1, baselineCorrectInterval);  % baseline correction
        tmpEEGM1 = eeg_checkset(tmpEEGM1);

        tmpEEGM0 = pop_epoch(EEG, desiredTrigsM0_cell, [epochStart epochEnd], 'newname', ['subsetTrigs'], 'epochinfo', 'yes');
        tmpEEGM0 = pop_rmbase(tmpEEGM0, baselineCorrectInterval);  % baseline correction
        tmpEEGM0 = eeg_checkset(tmpEEGM0);
        
        
        %******************************************************************
        % Select only the epochs that meet the thresholding specification
        %******************************************************************
        if ~isempty(threshEpoch)
            % Find epochs where the amplitude exceeds the thresholding value
            threshMaskM1 = squeeze(max(max(abs(tmpEEGM1.data),[],2)) < threshEpoch);
            threshMaskM0 = squeeze(max(max(abs(tmpEEGM0.data),[],2)) < threshEpoch);
            
            %-----------------------------------------------------
            %  Modify the pre-existing 'desiredIndices' to include
            %  only the trials that also meet the threshold criteria.
            %-----------------------------------------------------
            for i = 1:length(threshMaskM1)
                tmpInd = tmpEEGM1.event(i).urevent; % obtain the *absolute* trigger event index
                desiredIndicesM1(tmpInd) = desiredIndicesM1(tmpInd) && threshMaskM1(i);
            end
            
            for i = 1:length(threshMaskM0)
                tmpInd = tmpEEGM0.event(i).urevent; % obtain the *absolute* trigger event index
                desiredIndicesM0(tmpInd) = desiredIndicesM0(tmpInd) && threshMaskM0(i);
            end
            
            %-----------------------------------------------------
            %  Re-label the trigger values of the "desired" trials
            %  to account for the threshold criteria.
            %-----------------------------------------------------
            tmp = find(desiredIndicesM1 ~= 0);  % 'tmp' should be ONLY indicies of the desired triggers
            for i = 1:length(tmp)
                EEG.event(tmp(i)).type = 3;  % reset desired trigger to 2
            end

            tmp = find(desiredIndicesM0 ~= 0); % 'tmp' should be ONLY indicies of the desired triggers
            for i = 1:length(tmp)
                EEG.event(tmp(i)).type = 1;  % reset desired trigger to 2
            end

            % Define CELL ARRAY of trigger value to be epoched with 'pop_epoch'
            threshedTrigsM1_cell = {3};
            threshedTrigsM0_cell = {1};
            
            %-----------------------------------------------------
            % RE-epoch the data to select the epochs that meet the
            % thresholding criteria.
            %-----------------------------------------------------
            EEG_M1 = pop_epoch(EEG, threshedTrigsM1_cell, [epochStart epochEnd], 'newname', ['subsetTrigs'], 'epochinfo', 'yes');
            EEG_M1 = pop_rmbase(EEG_M1, baselineCorrectInterval);  % baseline correction
            EEG_M1 = eeg_checkset(EEG_M1);
            
            EEG_M0 = pop_epoch(EEG, threshedTrigsM0_cell, [epochStart epochEnd], 'newname', ['subsetTrigs'], 'epochinfo', 'yes');
            EEG_M0 = pop_rmbase(EEG_M0, baselineCorrectInterval);  % baseline correction
            EEG_M0 = eeg_checkset(EEG_M0);
            
        else
            % Thresholding wasn't applied, thus, use all trials that meet
            % the 'conditions' and 'responses' criteria.
            EEG_M1 = tmpEEGM1;
            EEG_M0 = tmpEEGM0;
        end
        


        %%%--------------------------------------
        %%% Try brickwall filter of alpha by setting all bands of alpha to zero
%         EEGNoAlpha = EEG;
%         EEGNoAlpha.data = removeAlpha(EEG.data, EEG.srate, alphaLim);
%         EEGNoAlpha = pop_rmbase(EEGNoAlpha, baselineCorrectInterval);  % baseline correction
%         
%         EEGbackNoAlpha = EEGback;        
%         EEGbackNoAlpha.data = removeAlpha(EEGback.data, EEGback.srate, alphaLim);        
%         EEGbackNoAlpha = pop_rmbase(EEGbackNoAlpha, baselineCorrectInterval);  % baseline correction

        
        %******************************************************************
        % Generate the MEAN epoch for the individual and save 
        %******************************************************************
        populationMeanDATA_M1(:,:,n,j) = mean(EEG_M1.data,3);
        populationMeanDATA_M0(:,:,n,j) = mean(EEG_M0.data,3);

%         populationMeanDATAnoAlpha_fig(:,:,n,j) = mean(EEGNoAlpha.data,3);
%         populationMeanDATAnoAlpha_back(:,:,n,j) = mean(EEGbackNoAlpha.data,3);
        
                
        %******************************************************************
        % Concatenate the individual epochs for each subject into one large
        % array.   
        %******************************************************************
        populationDATA_M1{j} = cat(3, populationDATA_M1{j}, EEG_M1.data);
        populationDATA_M0{j} = cat(3, populationDATA_M0{j}, EEG_M0.data);

%         populationDATAnoAlpha_fig{j} = cat(3, populationDATAnoAlpha_fig{j}, EEGNoAlpha.data);
%         populationDATAnoAlpha_back{j} = cat(3, populationDATAnoAlpha_back{j}, EEGbackNoAlpha.data);
        
        % Generate an array that provides a way to index which epochs
        % correspond to a given subjectID value.
        tmp = ones(size(EEG_M1.data,3), 1) * str2double(subjectID{n});   
        populationSUBJID_M1{j} = [populationSUBJID_M1{j}; tmp];
        numIndEpochsM1(j) = length(tmp);   % save the number of epochs used for each cuve of a given subject
        tmp = ones(size(EEG_M0.data,3), 1) * str2double(subjectID{n}); 
        populationSUBJID_M0{j} = [populationSUBJID_M0{j}; tmp];
        numIndEpochsM0(j) = length(tmp);   % save the number of epochs used for each cuve of a given subject
        
        %******************************************************************
        % Calculate the ERSP... plotting of individual data will occur here
        % if the flags are set appropriately.
        %******************************************************************
        if plotERSPplotsFLAG ~= 0
            if length(chanNum) == 1
                mainAnalysisPlotInd_erspITC

                % Save the data for the population average
                popERSPfig_ersp{n,j} = erspFig;
                popERSPfig_itc{n,j} = itcFig;
                popERSPfig_powbase{n,j} = powbaseFig;
                popERSPfig_times{n,j} = timesFig;
                popERSPfig_freqs{n,j} = freqsFig;

                popERSPback_ersp{n,j} = erspBack;
                popERSPback_itc{n,j} = itcBack;
                popERSPback_powbase{n,j} = powbaseBack;
                popERSPback_times{n,j} = timesBack;
                popERSPback_freqs{n,j} = freqsBack;
            else
                warning('DKR: Generation of the ERSP plot is not yet implemented for cases when the ''chanNum'' input parameter has more than one channel specified.')
            end
        end
        
        %%%-----------------------------------------
        %%% Plot the ERP image plot... helps to visualize temporal
        %%% coherence across the individual trials  
        if plotERPimageFLAG ~= 0
            if length(chanNum) == 1
                mainAnalysisPlotInd_erpimg
            else
                warning('DKR: Generation of the ERP image plot is not yet implemented for cases when the ''chanNum'' input parameter has more than one channel specified.')
            end
        end
            
       
        %******************************************************************
        % Reload the raw data into EEG and EEGback structures so that a
        % new curve can be generated for subjectID{n}
        %******************************************************************
        EEG = EEGRAW;
        
        % Waitbar to help the user know how long the process will take
        waitbar(ctr/(nCurves*length(subjectID)), hwait)
        ctr = ctr+1;
        
    end
    

    %**********************************************************************
    % Plotting of INDIVIDUAL data... see ERSP and ERP plots within the
    % nCurves for-loop for additional individual data plots.
    %**********************************************************************
    if plotIndividualDataFLAG ~= 0
        fs1 = EEG.srate;
        t1 = epochStart:1/fs1:(epochEnd-1/fs1);
        
        
       
        %%%-----------------------------------------
        %%% Time domain plot (raw data)
        if plotTimeDomainFLAG ~= 0
            mainAnalysisPlotInd_timedomRaw
        end
        
        %%%-----------------------------------------
        %%% Time domain plot  (for data without alpha frequencies)
        if plotTimeDomainAlphaModFLAG ~= 0
            mainAnalysisPlotInd_timedomAlphamod
        end
        
        %%%-----------------------------------------
        %%% DIFFERENCE WAVEFORM - Time domain plot (raw data) 
        if plotTimeDomainDiffPlotsFLAG ~= 0 && ((length(subsetResponses_detectFig) == 4) && (length(subsetResponses_detectBack) == 4))
            mainAnalysisPlotInd_timedomDiffraw
        end
        
        %%%-----------------------------------------
        %%% Scalp map DIFFERENCE plot
        if plotTopoDiffFLAG ~= 0 && ((length(subsetResponses_detectFig) == 4) && (length(subsetResponses_detectBack) == 4))
            mainAnalysisPlotInd_timedomScalpDiffraw
        end 
        
        if get(gcf, 'Number') > 25
            close all;
        end
    end
    
end




%**************************************************************************
% Begin plotting of POPULATION data
%**************************************************************************
if plotPopulationDataFLAG ~= 0
    fs1 = EEG.srate;
    t1 = epochStart:1/fs1:(epochEnd-1/fs1);

    %----------------------------------------------------------------------
    % Taking mean across all INDIVIDUAL epochs across all subjects
    %----------------------------------------------------------------------
    %%%-----------------------------------------
    %%% Time domain plot (raw data)
    if plotTimeDomainFLAG ~= 0
        mainAnalysisPlotPop_timedomRawIndepochs
    end
    
    %%%-----------------------------------------
    %%% Time domain plot  (for data without alpha frequencies)
    if plotTimeDomainAlphaModFLAG ~= 0
        mainAnalysisPlotPop_timedomAlphamodIndepochs
    end
    
    
    
    %----------------------------------------------------------------------
    % Taking mean across the MEAN epochs for each subject
    %----------------------------------------------------------------------
    %%%-----------------------------------------
    %%% Time domain plot (raw data)
    if plotTimeDomainFLAG ~= 0
        mainAnalysisPlotPop_timedomRawMeanepochs
    end

    %%%-----------------------------------------
    %%% Time domain plot  (for data without alpha frequencies)
    if plotTimeDomainAlphaModFLAG ~= 0
        mainAnalysisPlotPop_timedomAlphamodMeanepochs
    end
    
    %%%-----------------------------------------
    %%% DIFFERENCE WAVEFORM - Time domain plot (raw data) 
    if plotTimeDomainDiffPlotsFLAG ~= 0 && ((length(subsetResponses_detectFig) == 4) && (length(subsetResponses_detectBack) == 4))
        mainAnalysisPlotPop_timedomDiffrawMeanepochs
    end
    
    %%%-----------------------------------------
    %%% Scalp map DIFFERENCE plot
    if plotTopoDiffFLAG ~= 0 && ((length(subsetResponses_detectFig) == 4) && (length(subsetResponses_detectBack) == 4))
        mainAnalysisPlotPop_timedomScalpDiffrawMeanepochs
    end
    
    %%%-----------------------------------------
    %%% ERSP plot and ITC plot
    if plotERSPplotsFLAG ~= 0
        if length(chanNum) == 1   
            caxis([-1 1])  % adjust color scale mapping for image plots
            mainAnalysisPlotPop_ersp
            caxis([0 1])
            mainAnalysisPlotPop_itc
            caxis('auto')

            % If all four conditions were simulated, then the difference ERSP
            % plot can be generated
            if (length(subsetResponses_detectFig) == 4) && (length(subsetResponses_detectBack) == 4)
                mainAnalysisPlotPop_erspDiff
            end
        else
            warning('DKR: Generation of the ERSP and ERP image plot are not yet implemented for cases when the ''chanNum'' input parameter has more than one channel specified.')
        end
    end

end


% Close the waitbar
close(hwait)






%**************************************************************************
% Other plotting ideas
%**************************************************************************
% chanNum = 32;
% cycles = [1 0.5]; % The number of cycles in the window size
% epochlim = [-500 1000];  % in milliseconds
% freqRng = [4 20];
% padratio = 8;
% winsize = 512; 256;
% alpha = 0.1;
% timesout = 400;
% 
% 
% %%% Plot the event-related spectral perturbation (power spectrum)
% %-------------------------COMPARE 'attendFig' vs. 'attendBack' conditions
% figure;
% [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
%     newtimef({EEG.data(chanNum,:,:) EEGback.data(chanNum,:,:)},...
%     EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate,...
%     cycles,...
%     'freqs', freqRng,...
%     'padratio', padratio,...
%     'winsize', winsize,...
%     ...'alpha', alpha,...
%     'timesout', timesout);
% 
% 
% 
% %-------------------------
% % %  I think these settings are close to what Brigi and I came up with for
% % %    the 'attend FIGURE' data
% % figure; pop_newtimef(EEG, 1, 32, [-500  1000], [1 0.5] , 'baseline',[0],...
% %     'alpha',0.15, 'padratio', 8, 'winsize', 256, 'freqs', [4 20],...
% %     'plotphase', 'off', 'padratio', 8, 'basenorm', 'on', 'trialbase', 'full',...
% %     'timesout', 400);

