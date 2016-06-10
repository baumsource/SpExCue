# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 11:04:18 2016

@author: dkreed
"""
import sys
import dkrMNE
#dirString = '//projectnb//anl1//darrin//Exp1//'
#outputDir = '//projectnb//anl1//darrin//Exp1//_MNEpreProc//'
#outputDir = '//projectnb//anl1//darrin//Exp1//_MNEpreProc//2016April25_mne_eyeblinksRemoved//'
#outputDir = '//projectnb//anl1//darrin//Exp1//_MNEpreProc//2016April25_mne_eyeblinksNotRemoved//'
dirString = '//Users//rbaumgartner//Documents//ARI//ARIcloud//SpExCue//Experiments//Pilots//EEG\ 1\ \(cf\ Behavioral\ 1\)//'
outputDir = '//Users//rbaumgartner//Documents//ARI//ARIcloud//SpExCue//Experiments//Pilots//EEG\ 1\ \(cf\ Behavioral\ 1\)//SpExCue_analyzeEEGpilot_mne//'
expNum = 1
freqVals = [0.5, 25]
trigVals = [[3880, 3890], [3881, 3891], [3882, 3892], [3884, 3894]] # for exp1 
#trigVals = [3855, 3865, 3875, 3885] # for exp2
#channel = [32] #[5, 26, 31, 32]  # For exp1, Cz = 32; For exp2, Cz = 48
#i = int(sys.argv[1])
#dkrMNE.getTrigVals(dirString, subjString, dateString, expNum)
rejectThresh = None #120e-6
proj = False  # Apply eye blink SSP projection (True or False)


subjID_array = ['001', '010', '012', '013', '014', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024']
subjDate_array = ['01042016', # 001
            '01192016', # 010
            '12162015', # 012
            '01112016', # 013
            '01152016', # 014
            '01142016', # 015
            '01142016', # 016
            '01192016', # 017
            '01202016', # 018
            '01202016', # 019
            '02012016', # 020
            '01282016', # 021
            '01252016', # 022
            '02012016', # 023
            '01272016', # 024
            ]


##------------------
##  Run preprocessing on subject list.... note that 'i' is passed in through the SCC
##------------------
#dkrMNE.preProc(dirString, subjID_array[i], subjDate_array[i], expNum, freqVals, trigVals, outputDir, refChans=['EXG1', 'EXG2'], eogChan = ['A1'], rejectThresh=rejectThresh, proj=proj)
     
for i in range(len(subjID_array)):
    dkrMNE.preProc(dirString, subjID_array[i], subjDate_array[i], expNum, freqVals, trigVals, outputDir)
     

