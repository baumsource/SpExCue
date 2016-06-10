# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 11:15:48 2016

@author: dkreed
"""


import dkrMNE
import numpy as np
inputDir = 'C://Users//dkreed//Documents//EEGdata//___TEST DATA//Exp1//_MNEpreProc//2016April25_mne_eyeblinksNotRemoved//'
outputDir = 'C://Users//dkreed//Documents//EEGdata//___TEST DATA//Exp1//_MNEpreProc//2016April25_mne_eyeblinksNotRemoved//'
expNum = 1
savePlot = True

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



dkrMNE.eyeprojPlot(inputDir, subjID_array, expNum, outputDir=outputDir, savePlot=savePlot)
