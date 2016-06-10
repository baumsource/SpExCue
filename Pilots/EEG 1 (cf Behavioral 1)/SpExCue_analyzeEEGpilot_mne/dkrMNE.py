# -*- coding: utf-8 -*-

def preProc(dirString, subjString, dateString, expNum, freqVals, trigVals, outputDir=None, refChans=['EXG1', 'EXG2'], eogChan = ['A1'], rejectThresh=None, proj = False):
    """Preprocess a given participant's data and save data *for each trigger
       value* to file.

    Parameters
    ----------
    dirString :     STRING of the general directory of the experimental data,
                    i.e. one folder above the individual subject data.
    subjString :    STRING of the subject number, e.g. '001'
    dateString :    STRING of the date code used for the individual subjects'
                    data, e.g. '01042016'
    expNum :        INTEGER (or STRING) of the experiment number.  NOTE that 
                    '1' and '2' are currently the only valid values
    freqVals :      LIST of *two* frequency values in Hz that define the upper 
                    and lower cutoff filter frequencies. The values MUST be in
                    the proper order, i.e. [lowFreq highFreq]
    trigVals :      LIST of MNE trigger values, i.e. the 24-bit values, to be 
                    extracted from the data set.  It is possible to extract 
                    multiple conditions by providing a list. Furthermore, it is
                    possible to extract multiple trigger values per condition 
                    by creating a nested list, e.g. 
                        trigVals = [[3880, 3890], [3881, 3891], [3882, 3892]] 
    outputDir :     (optional) STRING of the directory in which the output data
                    should be saved. Default is the current working directory.
    refChans :      (optional) LIST of reference chanels to be used
    eogChan :       (optional) LIST of EOG channel to be used for eye blink 
                    artifact removal
    rejectThresh :  (optional) Threshold (in VOLTS, e.g. 120e-6) for rejection 
                    peak-to-peak threshold.  If 'None', thresholding is not
                    applied. 
    proj :          (optional) Apply the Signal Space Projection to remove the
                    eyeblink artifact
    

    Returns
    -------
    nothing... the data arrays are saved as .mat files. The output dictionary 
    is organized as follows:
    
    data :      the data array with dimesions
                    [numEpochs x 1 (i.e, numChannels) x numSamples]
    trigs :     full list of trigger values used for the given condition
    fs :        sampling rate of data
    tmin :      start of each epoch (###ms before the trigger)
    tmax :      end of each epoch (###ms after the trigger)    
    fLow :      LOWER cutoff frequency for bandpass filter applied to data
    fHigh :     UPPER cutoff frequency for bandpass filter applied to data
    baseline :  time (in seconds) to begin the baseline correction window 
    badChan :   unused external channels that were removed from the analysis
 
    """
    
    import mne  # Import MNE package for EEG analysis 
    from scipy import io  # used for saving (see io.savemat())
    from anlffr.helper import biosemi2mne as bs  # used for importing data with Hari's function
    import anlffr.preproc  # required for eye blink projection calculation
    from mne.preprocessing.ssp import compute_proj_epochs  # required for eye blink projection calculation

    #********************************************
    # Preprocess settings
    #********************************************
    fLow = freqVals[0]  # LOWER cutoff frequency for bandpass filter applied to data
    fHigh = freqVals[1]  # UPPER cutoff frequency for bandpass filter applied to data
    tmin = -0.5  # start of each epoch (###ms before the trigger)
    tmax = 1.0  # end of each epoch (###ms after the trigger)    
    baselineVal = -0.100 # Can be either None (default) or some NEGATIVE value corresponding to the time *prior* to t = 0 to begin the baseline correction window    
    # proj = True  # Indicate if the eyeblink projection should be applied    

    
    # Specify the unused external that should be removed from the analysis based on experiment number
    if expNum == 1:
        badChan = ['EXG4', 'EXG5', 'EXG6', 'EXG7', 'EXG8']  
    elif expNum == 2:
        badChan = ['EXG5', 'EXG6', 'EXG7', 'EXG8']  # Channels that should be removed from the analysis
    else:
        badChan = []
    
    # Define the output directory as the current working directory if none was specified
    if outputDir == None:
        import os
        outputDir = os.getcwd()
        

    #********************************************
    # Begin running preprocessing 
    #********************************************
#    directory = dirString + '\\subj' + subjString + '\\'
#    bdf_fnameFig = 'dkr_exp' + str(expNum) + '_subj' + subjString + '_' + dateString + '_detectFig.bdf'
#    bdf_fnameBack = 'dkr_exp' + str(expNum) + '_subj'+ subjString + '_' + dateString + '_detectBack.bdf'    
    
    directory = dirString + '//subj' + subjString + '//'
    bdf_fnameFig = 'dkr_exp' + str(expNum) + '_subj' + subjString + '_' + dateString + '_detectFig.bdf'
    bdf_fnameBack = 'dkr_exp' + str(expNum) + '_subj'+ subjString + '_' + dateString + '_detectBack.bdf'
    outDataString = '//exp' + str(expNum) + '_subj' + subjString + '_epoch'    # String used for the filename of the saved output file

#    ## Read in the BDF file (MNE methods)
#    rawFig = mne.io.read_raw_edf(directory+bdf_fnameFig, preload = True)    # NOTE:   directory+bdf_fname  will concatenate these two strings
#    rawBack = mne.io.read_raw_edf(directory+bdf_fnameBack, preload = True)    # NOTE:   directory+bdf_fname  will concatenate these two strings
#    rawFig_ref, _ = mne.io.set_eeg_reference(rawFig, ref_channels=refChans, copy=False) ## (re)-Reference the data with the first two external channels
#    rawBack_ref, _ = mne.io.set_eeg_reference(rawBack, ref_channels=refChans, copy=False)  ## (re)-Reference the data with the first two external channels    
        
    ## Read in the BDF file (Hari's method which automatically adds the Biosemi montage)
    rawFig_ref, evesFig = bs.importbdf(directory+bdf_fnameFig, refchans=refChans)
    rawBack_ref, evesBack = bs.importbdf(directory+bdf_fnameBack, refchans=refChans) 
    
    ## Remove the unused external channels
    rawFig_ref.drop_channels(badChan)
    rawBack_ref.drop_channels(badChan)
    
    ## Use a generic name for the dataset to be used in subsequent processing
    dataFig = rawFig_ref
    dataBack = rawBack_ref
    
    ## Obtain list of events as an integer array
    fs = dataFig.info['sfreq']
    trigDur = 0.001  # this was defined in my matlab experiment scripts, i.e. TDT definitions
    min_duration = trigDur - 1/float(fs)
    eventsFig = mne.find_events(dataFig, min_duration=min_duration)  
    eventsBack = mne.find_events(dataBack, min_duration=min_duration)
    

    #**************************************
    # Filter the data
    #**************************************
    # To use filter, it is necessary that the data to be processed is preloaded.
    # This can be done when loading the data by setting the input parameter, 
    # 'preload', to True.
    # e.g.,    raw = mne.io.read_raw_edf(dataDir+bdf_fname, preload = True)  
    dataFig.filter(fLow, fHigh, l_trans_bandwidth=fLow/2)
    dataBack.filter(fLow, fHigh,  l_trans_bandwidth=fLow/2) 
    
    
    #**************************************
    # Find eye blinks and create a signal-space projection
    #**************************************    
    eogFiltLow = 0.5  # highpass cutoff frequency in Hz for the signal to be used when finding the eye blinks
    eogFiltHigh = 10  # lowpass cutoff frequency in Hz for the signal to be used when finding the eye blinks
    eogThresh = 100e-6 # Threshold in VOLTS for finding eye blinks
    eventBlinkID = 998 # Event ID for indicating the eye blinks
    blinkWinPre = -0.20  # time in SECONDS prior to the estimated center of the eye blink to include in the SSP calculation
    blinkWinPost = 0.20  # time in SECONDS subsequent to the estimated center of the eye blink to include in the SSP calculation 
    
    blinksFig = anlffr.preproc.find_blinks(dataFig, event_id=eventBlinkID, ch_name=eogChan, thresh=eogThresh, l_freq=eogFiltLow, h_freq=eogFiltHigh)  # use Hari's code to find the eye blinks
    blinksBack = anlffr.preproc.find_blinks(dataBack, event_id=eventBlinkID, ch_name=eogChan, thresh=eogThresh, l_freq=eogFiltLow, h_freq=eogFiltHigh)  # use Hari's code to find the eye blinks
    tmpFig = mne.Epochs(dataFig, blinksFig, eventBlinkID, blinkWinPre, blinkWinPost, baseline=(None,0))  # select epochs around the eye blink
    tmpBack = mne.Epochs(dataBack, blinksBack, eventBlinkID, blinkWinPre, blinkWinPost, baseline=(None,0))  # select epochs around the eye blink
    blinksFigProj = mne.preprocessing.ssp.compute_proj_epochs(tmpFig)  # calculate the signal-space projection of the eye blinks
    blinksBackProj = mne.preprocessing.ssp.compute_proj_epochs(tmpBack)  # calculate the signal-space projection of the eye blinks
    
    ## Add projection(s) to the MNE-raw object BUT DO NOT APPLY YET... 
    ## the projection is applied during the epoching (assuming mne.Epochs() is
    ## called with proj=True)!
    dataFig.add_proj(blinksFigProj)    
    dataBack.add_proj(blinksBackProj)    
    
#    ##DEBUG
#    import time
#    print('************  DEBUG begin')
#    print(dataFig.info['projs'])    
#    print(dataBack.info['projs'])    
#    print(len(dataFig.ch_names))   
#    print(len(dataBack.ch_names))
#    print('************  DEBUG end')    
#    h = input('press 1, then <Enter>')
#    time.sleep(10)
#    h = input('press 1, then <Enter>')
    
    #------------------------------------------------------------------------------
    # PLOT eye blink projections  (CODE FROM LE - see April 5, 2016 email)
    #------------------------------------------------------------------------------
    plotProj = False
    savePlot = True    
    
    from mne.viz.topomap import _prepare_topo_plot as ptopo
    import pylab as pl  # used in my script for plotting
    
    picksFig, pos2dFig, merge_grads, ch_names, ch_type = ptopo(dataFig,'eeg',layout=None)   
    picksBack, pos2dBack, merge_grads, ch_names, ch_type = ptopo(dataBack,'eeg',layout=None)   
                                       
    proj0Fig = blinksFigProj[0]['data']['data']
    proj1Fig = blinksFigProj[1]['data']['data']
    proj0Back = blinksBackProj[0]['data']['data']
    proj1Back = blinksBackProj[1]['data']['data']
    
    if plotProj == True:
        pl.figure()
        pl.subplot(2,2,1)
        tp,cn = mne.viz.topomap.plot_topomap(proj0Fig.squeeze(), pos2dFig[:proj0Fig.shape[1],])
        pl.title('Subj' + subjString + ' | Proj 1, AttFig')
        pl.subplot(2,2,2)
        tp,cn = mne.viz.topomap.plot_topomap(proj1Fig.squeeze(), pos2dFig[:proj1Fig.shape[1],])
        pl.title('Subj' + subjString + ' | Proj 2, AttFig')
        pl.subplot(2,2,3)
        tp,cn = mne.viz.topomap.plot_topomap(proj0Back.squeeze(), pos2dBack[:proj0Back.shape[1],])
        pl.title('Subj' + subjString + ' | Proj 1, AttBack')
        pl.subplot(2,2,4)
        tp,cn = mne.viz.topomap.plot_topomap(proj1Back.squeeze(), pos2dBack[:proj1Back.shape[1],])
        pl.title('Subj' + subjString + ' | Proj 2, AttBack')

        if savePlot == True:
            pl.savefig(outputDir+'//subj' + subjString +'_projectionsPlot.pdf', format='pdf')   
    
        
        
    # Generate and save the dictionary that contains the data for the eyeblink 
    # projections for BOTH the figure and background data
    dictEyeproj = dict(proj0Fig = proj0Fig,   
                        proj1Fig = proj1Fig, 
                        pos2dFig = pos2dFig, 
                        proj0Back = proj0Back,  
                        proj1Back = proj1Back,  
                        pos2dBack = pos2dBack, 
                        subjString = subjString
                        )   
    io.savemat(outputDir + '//exp' + str(expNum) + '_subj' + subjString + '_eyeProj', dictEyeproj)    
    
    #**************************************
    # Extract epochs
    #**************************************
    baseline = (baselineVal, 0)  # define the baseline period... means from the first instant to t = 0
        
    for trg in range(len(trigVals)):
        if rejectThresh == None:        
            epochsFig = mne.Epochs(dataFig, eventsFig, trigVals[trg], tmin, tmax, proj=proj, baseline=baseline, preload=True, reject = None)  # Read epochs
            epochsBack = mne.Epochs(dataBack, eventsBack, trigVals[trg], tmin, tmax, proj=proj, baseline=baseline, preload=True, reject = None)  # Read epochs
            tmpThresh = []  # Because 'None' is not a valid data type to be saved to a Matlab file, a different value must be used
        else:
            epochsFig = mne.Epochs(dataFig, eventsFig, trigVals[trg], tmin, tmax, proj=proj, baseline=baseline, preload=True, reject = dict(eeg = rejectThresh))  # Read epochs
            epochsBack = mne.Epochs(dataBack, eventsBack, trigVals[trg], tmin, tmax, proj=proj, baseline=baseline, preload=True, reject = dict(eeg = rejectThresh))  # Read epochs
            tmpThresh = rejectThresh            
            
        dictFig = dict(data = epochsFig._data,   # the data array with dimesion [numEpochs x numChannels x numSamples]
                        trigs = trigVals[trg], # full list of trigger values used for the given condition
                        fs = epochsFig.info['sfreq'],  # sampling rate of data
                        tmin = tmin,    # start of each epoch (###ms before the trigger)
                        tmax = tmax,    # end of each epoch (###ms after the trigger)    
                        fLow = fLow,    # LOWER cutoff frequency for bandpass filter applied to data
                        fHigh = fHigh,  # UPPER cutoff frequency for bandpass filter applied to data
                        baseline = baselineVal, # time (in seconds) to begin the baseline correction window 
                        badChan = badChan,  # Channels that should be removed from the analysis
                        rejectThresh = tmpThresh  # Reject epochs based on peak-to-peak amplitude
                        )
        
        dictBack = dict(data = epochsBack._data,   # the data array with dimesion [numEpochs x 1 (i.e, numChannels) x numSamples]
                        trigs = trigVals[trg], # full list of trigger values used for the given condition
                        fs = epochsBack.info['sfreq'],  # sampling rate of data
                        tmin = tmin,    # start of each epoch (###ms before the trigger)
                        tmax = tmax,    # end of each epoch (###ms after the trigger)    
                        fLow = fLow,    # LOWER cutoff frequency for bandpass filter applied to data
                        fHigh = fHigh,  # UPPER cutoff frequency for bandpass filter applied to data
                        baseline = baselineVal, # time (in seconds) to begin the baseline correction window 
                        badChan = badChan,  # Channels that should be removed from the analysis
                        rejectThresh = tmpThresh  # Reject epochs based on peak-to-peak amplitude
                        )        
        
        # Save the dictionary that contains the data file.  NOTE that only the 
        # FIRST trigger value for a given condition is used in the filename
        io.savemat(outputDir+outDataString+str(trigVals[trg][0])+'_figMNE', dictFig)
        io.savemat(outputDir+outDataString+str(trigVals[trg][0])+'_backMNE', dictBack)
        
        # Save the MNE-epochs object so that there is more flexibility in subsequent processing
        epochsFig.save(outputDir+outDataString+str(trigVals[trg][0])+'_figMNE-epo.fif')
        epochsBack.save(outputDir+outDataString+str(trigVals[trg][0])+'_backMNE-epo.fif')
 

##*****************************************************************************
##*****************************************************************************
def getTrigVals(dirString, subjString, dateString, expNum):
    """Print to the console the MNE trigger values with the corresponding
       trigger values that I used in the actual experiment

    Parameters
    ----------
    dirString :     STRING of the general directory of the experimental data,
                    i.e. one folder above the individual subject data.
    subjString :    STRING of the subject number, e.g. '001'
    dateString :    STRING of the date code used for the individual subjects'
                    data, e.g. '01042016'
    expNum :        INTEGER (or STRING) of the experiment number.  NOTE that 
                    '1' and '2' are currently the only valid values

    Returns
    -------
    nothing... the information is printed to the console window
 
    """
    
    import mne  # Import MNE package for EEG analysis
    import numpy as np  # used in my script for creating the time vector
    
    directory = dirString + '//subj' + subjString + '//'
    bdf_fnameFig = 'dkr_exp' + str(expNum) + '_subj' + subjString + '_' + dateString + '_detectFig.bdf'
    bdf_fnameBack = 'dkr_exp' + str(expNum) + '_subj'+ subjString + '_' + dateString + '_detectBack.bdf'
    

    # Read in the BDF file
    rawFig = mne.io.read_raw_edf(directory+bdf_fnameFig, preload = False)    # NOTE:   directory+bdf_fname  will concatenate these two strings
    rawBack = mne.io.read_raw_edf(directory+bdf_fnameBack, preload = False)    # NOTE:   directory+bdf_fname  will concatenate these two strings

    # Obtain list of events as an integer array
    eventsFig = mne.find_events(rawFig)  
    eventsBack = mne.find_events(rawBack)
    
    
    ## Use only the lower 8 bits of events to determine the actual trigger 
    ## values that I used in the experiment and print these updated trigger ID 
    ## values to the screen 
    print('*******************************')        
    print('Trigger values for "Fig" data')
    print('*******************************')        
    un = np.unique(eventsFig[:,2])
    for i in un:
        tmp = bin(i)
        tmp = tmp[(len(tmp)-8):len(tmp)]
        tmpDec = int(tmp,2)
        print(i)    
        print(tmpDec)
        print('-----------')
    
    print('*******************************')        
    print('Trigger values for "Back" data')
    print('*******************************') 
    un = np.unique(eventsBack[:,2])
    for i in un:
        tmp = bin(i)
        tmp = tmp[(len(tmp)-8):len(tmp)]
        tmpDec = int(tmp,2)
        print(i)    
        print(tmpDec)
        print('-----------')





#*****************************************************************************
#*****************************************************************************
def grandAvgEpoch(inputDir, subjList, triggerNum, expNum, expCond, channel, plotMe=0, outputDir=None, dirString=None):   
    """Access the output data from dkrMNE.proProc() and compute a grand-average
        epoch across the input 'subjList' for a *single* epoch.

    Parameters
    ----------
    inputDir :      STRING of the directory where the dkrMNE.preProc() data are
                    located.  If specified as 'None', use the optional input
                    'dirString'.
    subjList :      LIST of STRING of the subject numbers, e.g. ['001', '012']
    triggerNum :    INTEGER of event number to be used in the grand average.
                    If multiple events were used for a given condition, NOTE 
                    that it is essential that the *first* trigger event is 
                    specified here, i.e. the trigger number used in the
                    filename of the data output from dkrMNE.preProc()
    expNum :        INTEGER (or STRING) of the experiment number.  NOTE that 
                    '1' and '2' are currently the only valid values
    expCond :       STRING of the experiment condition, i.e. 'fig' or 'back'. 
    channel :       LIST of integers) for the BioSemi channel number(s) to be 
                    used in the analysis.  The *mean* of these channels will be
                    used for the grandAvgEpoch() calculation
    plotMe :        (optional) INTEGER flag. If not zero, plot the grand 
                    average epoch for the given triggerNumber across all 
                    subjects in 'subjList'
    outputDir :     (optional) STRING of the directory in which the output data
                    should be saved. Default is the current working directory.
    dirString :     (optional) STRING.  This can be used iff 'inputDir' is
                    specified as 'None'. If so, the output from 
                    dkrMNE.preProc() is saved within individual folders labeled
                    for each subject and 'dirString' is the directory one level
                    above the individual subject directories.
    

    Returns
    -------
    dictOut :  DICT of output data and parameters.  The output dictionary is 
               organized as follows:
        data ->     The data array with dimesions
                    [numEpochs x 1 (i.e, numChannels) x numSamples]
        subjList -> List of subjects that compose the grand average epoch
        channel->   The channels used to generate the grand average epoch.
                    Note that these are the ACTUAL biosemi channel number, 
                    i.e. NOT the python array index values.
        trigs ->    Full list of trigger values used for the given condition
        fs ->       Sampling rate of data
        tmin ->     Start of each epoch (###ms before the trigger)
        tmax ->     End of each epoch (###ms after the trigger)    
        fLow ->     LOWER cutoff frequency for bandpass filter applied to data
        fHigh ->    UPPER cutoff frequency for bandpass filter applied to data
        baseline -> Time (in seconds) to begin the baseline correction window 
        badChan ->  Unused external channels that were removed from the
                    analysis
 
    """
    from scipy import io  # used for saving (see io.savemat())
    import numpy as np
        
    # Define the output directory as the current working directory if none was specified
    if outputDir == None:
        import os
        outputDir = os.getcwd()
        
        
    chanN = [x-1 for x in channel]  # *index* of the channel number, i.e. channel number minus one    
    
    # Create a string to be used for the filename of the output data
    outDataString = '//exp' + str(expNum) + '_GrandAvg_epoch' + str(triggerNum)  + '_Ch' + str(channel[0]) + '_' + expCond.lower()  # String used for the filename of the saved output file    
    
    
    
    # Loop through all subjects to obtain a mean epoch for each subject 
    for i in range(len(subjList)):    
        
        # Handle the optional input directory string
        if inputDir != None:
            directory = inputDir
        else:
            directory = dirString + '\\subj' + subjList[i] + '\\'
        
        mat_fname = 'exp' + str(expNum) + '_subj' + subjList[i] + '_epoch' + str(triggerNum) + '_' + expCond.lower() + 'MNE.mat'
        
           
        # Load the data for a given subject
        tmp = io.loadmat(directory+mat_fname)
        tmpDat = tmp['data']
        
        # Initialize the grand average output array on the first run through the for-loop
        if i == 0:
            grandAvg = [0] * len(tmpDat[0,0,:])
        
        
        # B/C all channels were used in the epoching, select the desired channels now
        tmpDat = tmpDat[:,chanN,:]  
        
        # Take the mean across the selected channels
        tmpDat = np.mean(tmpDat,1)
        
        # Take the mean across all epochs
        tmpMean = np.mean(tmpDat,0)
        
        # Add the given subject's mean epoch to the main grand average output array
        grandAvg = grandAvg + tmpMean
        
        # DEBUG
        print(str(subjList[i]) + ' -> numEpochs ' + str(len(tmpDat[:,0])))
    
    
    # Create the actual grand average epoch across all subjects
    grandAvg = grandAvg/len(subjList)
    
    
    dictOut = dict(data = grandAvg,   # the data array with dimesion [numSamples]
                    subjList = subjList, # list of subjects that compose the grand average epoch
                    channel = channel,  # The ACTUAL biosemi channel number, i.e. NOT the python array index value
                    trigs = tmp['trigs'], # full list of trigger values used for the given condition
                    fs = tmp['fs'],        # sampling rate of data
                    tmin = tmp['tmin'],    # start of each epoch (###ms before the trigger)
                    tmax = tmp['tmax'],    # end of each epoch (###ms after the trigger)    
                    fLow = tmp['fLow'],    # LOWER cutoff frequency for bandpass filter applied to data
                    fHigh = tmp['fHigh'],  # UPPER cutoff frequency for bandpass filter applied to data
                    baseline = tmp['baseline'], # time (in seconds) to begin the baseline correction window 
                    badChan = tmp['badChan']  # Channels that should be removed from the analysis
                    ) 
                 
                 
    if plotMe != 0:                
        t = np.arange(tmp['tmin'],tmp['tmax']+(1/tmp['fs']),(1/float(tmp['fs'])))  # ****  Because np.arrange does not include the final data point, it is necessary to manually add an additional sample (1/fs) to the time array
        import pylab as pl  # used in my script for plotting
        pl.figure()
        pl.plot(t, grandAvg)
        pl.xlabel('Time (s)')
        pl.ylabel('Amplitude (microvolts)')
        pl.title('Grand average')
        pl.ticklabel_format(style='sci', scilimits = (-3, 1), axis='y')  
    

    io.savemat(outputDir+outDataString, dictOut) 
    return(dictOut)
    

#*****************************************************************************
#*****************************************************************************
def grandAvgPlot(inputDir, subjList, trigList, expNum, channel, outputDir=None, savePlot=False):
    """Access the output data from dkrMNE.preProc() and compute an average
        epoch across the input 'subjList' for the list of conditions specified
        by 'trigList'.  Note that a grand average plot for *each subject* will 
        also be generated.
        

    Parameters
    ----------
    inputDir :      STRING of the directory where the .mat data from 
                    dkrMNE.preProc() are located.
    subjList :      LIST of strings for the subjectIDs to be used in the 
                    desired grand average plot.  
    trigList :      LIST of strings or integers OR a numpy array of integers...
                    of event number to be used in the grand average. If 
                    multiple events were used for a given condition, NOTE that 
                    it is essential that the *first* trigger event is specified 
                    here, i.e. the trigger number used in the filename of the 
                    data output from dkrMNE.preProc()
    expNum :        INTEGER (or STRING) of the experiment number.  NOTE that 
                    '1' and '2' are currently the only valid values
    channel :       LIST of integer(s) for the BioSemi channel number(s) to be 
                    used in the analysis.  The *mean* of these channels will be 
                    used for the grandAvgPlot() calculation
    outputDir :     (optional) STRING of the directory in which the output data
                    should be saved. Default is the current working directory. 
    savePlot :      (optional) BOOLEAN to turn on/off the saving of all plots
                    as PDFs.  Note that savePlot can be True/False or !0/0

    Returns
    -------
    dictOut :  Nothing.... a PDF of the plot is saved to the specified output 
                directory
     
    """    

    from scipy import io  # used for saving (see io.savemat())
    import numpy as np
    import pylab as pl  # needed for plotting    
    
    # Define the output directory as the current working directory if none was specified
    if outputDir == None:
        import os
        outputDir = os.getcwd()
    
    colorOrder = ['b', 'r', 'g', 'm', 'c']  # order of colors to be used for the different trigerList values
    chanN = [x-1 for x in channel]  # *index* of the channel number, i.e. channel number minus one   
    
    for j in range(len(subjList)):      
        pl.figure()  # create a new figure window
        
        #-----------------------------------------------------
        # Plot the FIGURE data    
        #-----------------------------------------------------
        for i in range(len(trigList)):  
            
            # Load the data for a given subject   
            tmp = io.loadmat(inputDir+'//exp' + str(expNum) + '_subj' + subjList[j] + '_epoch' + str(trigList[i]) + '_figMNE.mat')
            tmpDat = tmp['data']
            
            # Initialize the grand average output array on the first run through the for-loop
            if (i == 0) & (j == 0):
                grandAvgFig = np.zeros((len(trigList), len(subjList), len(tmpDat[0,0,:])))
            
            # B/C all channels were used in the epoching, select the desired channels now
            tmpDat = tmpDat[:,chanN,:]  
            # Take the mean across the selected channels
            tmpDat = np.mean(tmpDat,1)
            # Take the mean across all epochs
            tmpMean = np.mean(tmpDat,0)
            # Insert the given subject's mean epoch to the main grand average output array
            grandAvgFig[i,j,:] = tmpMean
    
            # get the sampling frequency and the vector of time points for plotting
            if (i == 0) & (j == 0):
                fs = tmp['fs']  # get the sampling frequency
                t = np.arange(tmp['tmin'], tmp['tmax']+(1/fs),(1/float(fs)))  # ****  Because np.arrange does not include the final data point, it is necessary to manually add an additional sample (1/fs) to the time array
            
            # Plot the trace across all subjects for the given trigger number
            pl.plot(t, tmpMean, '-'+colorOrder[i])        
            if (i == 0) & (j == 0):
                pl.hold(True)
    
        #-----------------------------------------------------
        # Plot the BACKGROUND data    
        #-----------------------------------------------------    
        for i in range(len(trigList)):  
            tmp = io.loadmat(inputDir+'//exp' + str(expNum) + '_subj' + subjList[j] + '_epoch' + str(trigList[i]) + '_backMNE.mat')
            tmpDat = tmp['data']
            
            # Initialize the grand average output array on the first run through the for-loop
            if (i == 0) & (j == 0):
                grandAvgBack = np.zeros((len(trigList), len(subjList), len(tmpDat[0,0,:])))
            
            # B/C all channels were used in the epoching, select the desired channels now
            tmpDat = tmpDat[:,chanN,:]  
            # Take the mean across the selected channels
            tmpDat = np.mean(tmpDat,1)
            # Take the mean across all epochs
            tmpMean = np.mean(tmpDat,0)
            # Insert the given subject's mean epoch to the main grand average output array
            grandAvgBack[i,j,:] = tmpMean
    
            # get the sampling frequency and the vector of time points for plotting
            if (i == 0) & (j == 0):
                fs = tmp['fs']  # get the sampling frequency
                t = np.arange(tmp['tmin'], tmp['tmax']+(1/fs),(1/float(fs)))  # ****  Because np.arrange does not include the final data point, it is necessary to manually add an additional sample (1/fs) to the time array
            
            # Plot the trace across all subjects for the given trigger number
            pl.plot(t, tmpMean, '--'+colorOrder[i])        
            if (i == 0) & (j == 0):
                pl.hold(True)
        
        #----------------------------------------------------
        # Final stuff for individual subject plot    
        #----------------------------------------------------
        pl.xlabel('Time (s)')
        pl.ylabel('Amplitude (volts)')
        pl.title('Grand average for Subj' + subjList[j])
        pl.ticklabel_format(style='sci', scilimits = (-3, 1), axis='y')  
        
        if savePlot == True:
            pl.savefig(outputDir+'//subj' + subjList[j]+'_GAplot.pdf', format='pdf')    
    
    
    
    #------------------------------------------------------------------------------
    # Now plot the grandAverage across the subjList
    #------------------------------------------------------------------------------
    
    pl.figure()  # create a new figure window
    GAfig = np.mean(grandAvgFig,1)   # Take mean across the subjList
    GAback = np.mean(grandAvgBack,1)   # Take mean across the subjList
    
    for i in range(len(trigList)):  
        tmp = np.squeeze(GAfig[i,:])
        pl.plot(t, tmp, '-'+colorOrder[i])  
        pl.hold(True)
        
    for i in range(len(trigList)):  
        tmp = np.squeeze(GAback[i,:])
        pl.plot(t, tmp, '--'+colorOrder[i])
    
        
    pl.xlabel('Time (s)')
    pl.ylabel('Amplitude (volts)')
    pl.title('Grand average')
    pl.ticklabel_format(style='sci', scilimits = (-3, 1), axis='y')  
    
    if savePlot == True:
        pl.savefig(outputDir+'//subjCOMBO_GAplot.pdf', format='pdf')       
    

#*****************************************************************************
#*****************************************************************************
def grandAvgPlot_Legacy(inputDir, trigList, expNum, channel, outputDir=None, savePlot=False):
#def grandAvgPlot(inputDir, trigList, expNum, expCond, channel, outputDir=None): 
    """Access the output data from dkrMNE.grandAvgEpoch() and  generate a
        plot with the grand-average epochs across the 'trigList'.  
        
    SUGGEST to use dkrMNE.grandAvgPlot()

    Parameters
    ----------
    inputDir :      STRING of the directory where the data from grandAvgEpoch()
                    are located.  
    trigList :      numpy array of INTEGERS.... of event number to be used in 
                    the grand average. If multiple events were used for a given 
                    condition, NOTE that it is essential that the *first* 
                    trigger event is specified here, i.e. the trigger number 
                    used in the filename of the data output from 
                    dkrMNE.preProc()
    expNum :        INTEGER (or STRING) of the experiment number.  NOTE that 
                    '1' and '2' are currently the only valid values
    channel :       INTEGER for the BioSemi channel number used for the 
                    *filename* in dkrMNE.grandAvgEpoch().
    outputDir :     (optional) STRING of the directory in which the output data
                    should be saved. Default is the current working directory.
    savePlot :      (optional) BOOLEAN to turn on/off the saving of all plots
                    as PDFs.  Note that savePlot can be True/False or !0/0
                    
    

    Returns
    -------
    dictOut :  Nothing.... a PDF of the plot is saved to the specified output 
                directory
     
    """    

    from scipy import io  # used for saving (see io.savemat())
    import numpy as np
    import pylab as pl  # needed for plotting    
    
    colorOrder = ['b', 'r', 'g', 'm']
    pl.figure()  # create a new figure window
    
    
    #-----------------------------------------------------
    # Plot the FIGURE data    
    #-----------------------------------------------------
    for i in range(len(trigList)):  
#        tmp = io.loadmat(inputDir+'//exp' + str(expNum) + '_GrandAvg_epoch' + trigList[i] + '_fig')
        tmp = io.loadmat(inputDir+'//exp' + str(expNum) + '_GrandAvg_epoch' + trigList[i] + '_Ch' + str(channel) + '_fig')
        dat = tmp['data']
        
        if i == 0:
            fs = tmp['fs']  # get the sampling frequency
            t = np.arange(tmp['tmin'], tmp['tmax']+(1/fs),(1/float(fs)))  # ****  Because np.arrange does not include the final data point, it is necessary to manually add an additional sample (1/fs) to the time array
        
                
        pl.plot(t, dat[0,:], '-'+colorOrder[i])
        pl.hold(True)
        
        
    #-----------------------------------------------------
    # Plot the BACKGROUND data    
    #-----------------------------------------------------    
    for i in range(len(trigList)):  
#        tmp = io.loadmat(inputDir+'//exp' + str(expNum) + '_GrandAvg_epoch' + trigList[i] + '_back')
        tmp = io.loadmat(inputDir+'//exp' + str(expNum) + '_GrandAvg_epoch' + trigList[i] + '_Ch' + str(channel) + '_back')
        dat = tmp['data']
        
        if i == 0:
            fs = tmp['fs']  # get the sampling frequency
            t = np.arange(tmp['tmin'], tmp['tmax']+(1/fs),(1/float(fs)))  # ****  Because np.arrange does not include the final data point, it is necessary to manually add an additional sample (1/fs) to the time array
        
                
        pl.plot(t, dat[0,:], '--'+colorOrder[i])
        pl.hold(True)
        
    
    #----------------------------------------------------
    # Final stuff for plot    
    #----------------------------------------------------
    pl.xlabel('Time (s)')
    pl.ylabel('Amplitude (microvolts)')
    pl.title('Grand average')
    pl.ticklabel_format(style='sci', scilimits = (-3, 1), axis='y')  
    
    if savePlot == True:
        if outputDir == None:
            pl.savefig('comboGAplot.pdf', format='pdf')    
        else:
            pl.savefig(outputDir+'//comboGAplot.pdf', format='pdf')    


     
    
        
#*****************************************************************************
#*****************************************************************************
def grandAvgITC(inputDir, subjList, triggerNum, expNum, expCond, channel, freqs, plotMe=0, outputDir=None, dirString=None):   
    """Access the output data from dkrMNE.proProc() and compute the 
        grand-average inter-trial coherence (ITC) across the input 'subjList' 
        for a *single* epoch.

    Parameters
    ----------
    inputDir :      STRING of the directory where the dkrMNE.preProc() data are
                    located.  If specified as 'None', use the optional input
                    'dirString'.
    subjList :      LIST of STRING of the subject numbers, e.g. ['001', '012']
    triggerNum :    INTEGER of event number to be used in the grand average.
                    If multiple events were used for a given condition, NOTE 
                    that it is essential that the *first* trigger event is 
                    specified here, i.e. the trigger number used in the
                    filename of the data output from dkrMNE.preProc()
    expNum :        INTEGER (or STRING) of the experiment number.  NOTE that 
                    '1' and '2' are currently the only valid values
    expCond :       STRING of the experiment condition, i.e. 'fig' or 'back'. 
    channel :       INTEGER for the BioSemi channel number to be used in the 
                    analysis.  NOTE that only a *single channel* can be 
                    specified here.
    freqs :         INTEGER LIST containing the frequencies on which the ITC
                    should be calculated.
    plotMe :        (optional) INTEGER flag. If not zero, plot the grand 
                    average epoch
    outputDir :     (optional) STRING of the directory in which the output data
                    should be saved. Default is the current working directory.
    dirString :     (optional) STRING.  This can be used iff 'inputDir' is
                    specified as 'None'. If so, the output from 
                    dkrMNE.preProc() is saved within individual folders labeled
                    for each subject and 'dirString' is the directory one level
                    above the individual subject directories.
    

    Returns
    -------
    dictOut :  DICT of output data and parameters.  The output dictionary is 
               organized as follows:
        data ->     The data array with dimesions
                    [numEpochs x 1 (i.e, numChannels) x numSamples]
        subjList -> List of subjects that compose the grand average epoch
        channel->   The channels used to generate the grand average epoch.
                    Note that these are the ACTUAL biosemi channel number, 
                    i.e. NOT the python array index values.
        trigs ->    Full list of trigger values used for the given condition
        fs ->       Sampling rate of data
        tmin ->     Start of each epoch (###ms before the trigger)
        tmax ->     End of each epoch (###ms after the trigger)    
        fLow ->     LOWER cutoff frequency for bandpass filter applied to data
        fHigh ->    UPPER cutoff frequency for bandpass filter applied to data
        baseline -> Time (in seconds) to begin the baseline correction window 
        badChan ->  Unused external channels that were removed from the
                    analysis
 
    """
    
    from scipy import io  # used for saving (see io.savemat())
#    import numpy as np
    import mne.time_frequency.tfr as tfr  # required for the ITC calculation
        
    # Define the output directory as the current working directory if none was specified
    if outputDir == None:
        import os
        outputDir = os.getcwd()
        
    
    # Create a string to be used for the filename of the output data
    outDataString = '//exp' + str(expNum) + '_GrandAvgITC_epoch' + str(triggerNum) + '_' + expCond.lower() + '_ch' + str(channel)  # String used for the filename of the saved output file    
    
    
    
    # Loop through all subjects to obtain a mean epoch for each subject 
    for i in range(len(subjList)):    
        
        # Handle the optional input directory string
        if inputDir != None:
            directory = inputDir
        else:
            directory = dirString + '//subj' + subjList[i] + '//'
        
        mat_fname = 'exp' + str(expNum) + '_subj' + subjList[i] + '_epoch' + str(triggerNum) + '_' + expCond.lower() + 'MNE.mat'
        
           
        # Load the data for a given subject
        tmp = io.loadmat(directory+mat_fname)
        tmpDat = tmp['data']
        

        #---------------------------------------------------------------
        # Run intertrial coherence analysis (as suggested by Le and only 
        # requires the data array not the entire Epochs object)
        #---------------------------------------------------------------
        fs = tmp['fs']  # get the sampling frequency         
        power, itc = tfr._induced_power_mtm(tmpDat, fs, freqs)
        
        # Initialize the grand average output array on the first run through the for-loop
        if i == 0:
            gaITC = itc
            gaPower = power
        else:
            gaITC = gaITC + itc
            gaPower = gaPower + power

        # DEBUG
        print(subjList[i])
        
        ##----------------------------------------------
        ## Save the individual subject ITC data to file        
        ##----------------------------------------------
        outDataStringSubj = '//exp' + str(expNum) + '_subj' + subjList[i] + '_epoch' + str(triggerNum) + '_' + expCond.lower() + '_ch' + str(channel-1) + 'ITC'  # String used for the filename of the saved output file            
        dictOutSubj = dict(dataITC = itc,   # the data array with dimesion [numChannels, len(freqs), numSamples]
                    dataPower = power,   # the data array with dimesion [numChannels, len(freqs), numSamples]
                    freqs = freqs, # list of frequencies used for the ITC calculation
                    subjList = subjList[i], # list of subjects that compose the grand average epoch
                    channel = channel,  # The ACTUAL biosemi channel number, i.e. NOT the python array index value
                    trigs = triggerNum, # trigger values used for the given condition ITC plot
                    fs = tmp['fs'],        # sampling rate of data
                    tmin = tmp['tmin'],    # start of each epoch (###ms before the trigger)
                    tmax = tmp['tmax'],    # end of each epoch (###ms after the trigger)    
                    fLow = tmp['fLow'],    # LOWER cutoff frequency for bandpass filter applied to data
                    fHigh = tmp['fHigh'],  # UPPER cutoff frequency for bandpass filter applied to data
                    baseline = tmp['baseline'], # time (in seconds) to begin the baseline correction window 
                    badChan = tmp['badChan']  # Channels that should be removed from the analysis
                    )    
        io.savemat(outputDir+outDataStringSubj, dictOutSubj)  
    
    
    
    
    
    # Create the actual grand average epoch across all subjects
    gaITC = gaITC/len(subjList)
    gaPower = gaPower/len(subjList)
    
    
    dictOut = dict(dataITC = gaITC,   # the data array with dimesion [numChannels, len(freqs), numSamples]
                    dataPower = gaPower,   # the data array with dimesion [numChannels, len(freqs), numSamples]
                    freqs = freqs, # list of frequencies used for the ITC calculation
                    subjList = subjList, # list of subjects that compose the grand average epoch
                    channel = channel,  # The ACTUAL biosemi channel number, i.e. NOT the python array index value
                    trigs = triggerNum, # trigger values used for the given condition ITC plot
                    fs = tmp['fs'],        # sampling rate of data
                    tmin = tmp['tmin'],    # start of each epoch (###ms before the trigger)
                    tmax = tmp['tmax'],    # end of each epoch (###ms after the trigger)    
                    fLow = tmp['fLow'],    # LOWER cutoff frequency for bandpass filter applied to data
                    fHigh = tmp['fHigh'],  # UPPER cutoff frequency for bandpass filter applied to data
                    baseline = tmp['baseline'], # time (in seconds) to begin the baseline correction window 
                    badChan = tmp['badChan']  # Channels that should be removed from the analysis
                    )      
    
    if plotMe != 0:
        tmin = float(tmp['tmin'])
        tmax = float(tmp['tmax'])
        fmin = float(min(freqs))
        fmax = float(max(freqs))
        
        import pylab as pl  # used in my script for plotting
        pl.figure()
        pl.imshow(gaITC[channel-1,:,:].squeeze(),aspect='auto',origin=(tmin,fmin), extent=[tmin, tmax, fmin, fmax])
        pl.xlabel('Time (s)')
        pl.ylabel('Frequency (Hz)')
        pl.title('Grand average ITC - Ch. ' + str(channel) + ' Trig. ' + str(triggerNum) + ' (' + expCond.upper() + ')')


        pl.figure()
        pl.imshow(gaPower[channel-1,:,:].squeeze(),aspect='auto',origin=(tmin,fmin), extent=[tmin, tmax, fmin, fmax]) 
        pl.xlabel('Time (s)')
        pl.ylabel('Frequency (Hz)')
        pl.title('Grand average induced power - Ch. ' + str(channel) + ' Trig. ' + str(triggerNum)+ ' (' + expCond.upper() + ')')
    
    
    io.savemat(outputDir+outDataString, dictOut)  
    return(dictOut)
    
    
    
#*****************************************************************************
#*****************************************************************************
def eyeprojPlot(inputDir, subjList, expNum, outputDir=None, savePlot=False):
    """Access the output data from dkrMNE.preProc() and plot the two primary
    Signal-Space Projections (SSP) of the eye blinks.
        

    Parameters
    ----------
    inputDir :      STRING of the directory where the .mat data from 
                    dkrMNE.preProc() are located.
    subjList :      LIST of strings for the subjectIDs to be used in the 
                    desired grand average plot.  
    expNum :        INTEGER (or STRING) of the experiment number.  NOTE that 
                    '1' and '2' are currently the only valid values
    outputDir :     (optional) STRING of the directory in which the output data
                    should be saved. Default is the current working directory. 
    savePlot :      (optional) BOOLEAN to turn on/off the saving of all plots
                    as PDFs.  Note that savePlot can be True/False or !0/0

    Returns
    -------
    Nothing....  just plots the eye blink projections in a single figure
     
    """   
    import mne
#    from mne.viz.topomap import _prepare_topo_plot as ptopo
#    import mne.viz.topomap.plot_topomap
    import pylab as pl  # used in my script for plotting
    from scipy import io  # used for saving (see io.savemat())
    
    
    # Define the output directory as the current working directory if none was specified
    if outputDir == None:
        import os
        outputDir = os.getcwd()
    
    #------------------------------------------------------------------------------
    # PLOT eye blink projections  (CODE FROM LE - see April 5, 2016 email)
    #------------------------------------------------------------------------------
    for j in range(len(subjList)):      
        pl.figure()  # create a new figure window
          
        # Load the data for a given subject   
        tmp = io.loadmat(outputDir + '//exp' + str(expNum) + '_subj' + subjList[j] + '_eyeProj.mat')    
        
        proj0Fig = tmp['proj0Fig']
        proj1Fig = tmp['proj1Fig']
        proj0Back= tmp['proj0Back']
        proj1Back= tmp['proj1Back']
        pos2dFig = tmp['pos2dFig']
        pos2dBack = tmp['pos2dBack']
        subjString = tmp['subjString']
        subjString  = subjString[0]
            
#        picksFig, pos2dFig, merge_grads, ch_names, ch_type = ptopo(dataFig,'eeg',layout=None)   
#        picksBack, pos2dBack, merge_grads, ch_names, ch_type = ptopo(dataBack,'eeg',layout=None)   
                                           
#        proj0Fig = blinksFigProj[0]['data']['data']
#        proj1Fig = blinksFigProj[1]['data']['data']
#        proj0Back = blinksBackProj[0]['data']['data']
#        proj1Back = blinksBackProj[1]['data']['data']
        
        pl.figure()
        pl.subplot(2,2,1)
        tp,cn = mne.viz.topomap.plot_topomap(proj0Fig.squeeze(), pos2dFig[:proj0Fig.shape[1],])
        pl.title('Subj' + subjString + ' | Proj 1, AttFig')
        pl.subplot(2,2,2)
        tp,cn = mne.viz.topomap.plot_topomap(proj1Fig.squeeze(), pos2dFig[:proj1Fig.shape[1],])
        pl.title('Subj' + subjString + ' | Proj 2, AttFig')
        pl.subplot(2,2,3)
        tp,cn = mne.viz.topomap.plot_topomap(proj0Back.squeeze(), pos2dBack[:proj0Back.shape[1],])
        pl.title('Subj' + subjString + ' | Proj 1, AttBack')
        pl.subplot(2,2,4)
        tp,cn = mne.viz.topomap.plot_topomap(proj1Back.squeeze(), pos2dBack[:proj1Back.shape[1],])
        pl.title('Subj' + subjString + ' | Proj 2, AttBack')
    
        if savePlot == True:
            pl.savefig(outputDir+'//subj' + subjString +'_projectionsPlot.pdf', format='pdf')  
