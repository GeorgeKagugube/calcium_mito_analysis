#!/Users/gwk/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 10:00:00 2023"

Python script for analysing Fura-2 AM calcium imaging data from a dataframe
This script uses the package to perform analysis on Fura-2 AM calcium imaging data.
It reads data from a CSV file, processes it, and generates plots to visualize the results.
The results are saved to a CSV file and the plots are displayed using matplotlib.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from scipy.integrate import simpson
from scipy.optimize import curve_fit, OptimizeWarning
import seaborn as sns
import warnings
import os
import sys
from pathlib import Path

## Import own written functions needed in the analysis of the data
output_dir = Path('/Users/gwk/Desktop/Bioinformatics/calcium_mito_analysis/results')
data_files1 = Path('/Users/gwk/Desktop/Bioinformatics/calcium_mito_analysis/data/Platting_1/csv_files/Day1Recording_M3_Hom_Fura_2Ratio.csv')
path_to_rhod = Path('/Users/gwk/Desktop/Bioinformatics/calcium_mito_analysis/data/Platting_1/csv_files/Day1Recording_M3_Hom_rhod123.csv')

from utils.io_utils import load_data, insertTime, cleanVarNames
from utils.plotting import plot_fura2_rhodamine
from utils.fura_analysis import analyze_single_trace, analyze_all_cells

# Set plotting style
sns.set_theme(style='whitegrid', palette='muted')
# Suppress warnings from scipy.optimize
warnings.simplefilter('ignore', OptimizeWarning)
# Set default figure size
plt.rcParams['figure.figsize'] = (10, 6)    
# Set default font size
plt.rcParams['font.size'] = 12  
# Set default line width
plt.rcParams['lines.linewidth'] = 2

## Load the data to be analysed here
hom1 = load_data(data_files1)
rhod1 = load_data(path_to_rhod)

# Check the header of the files here
print (hom1.head())
print (rhod1.head())

## Insert time in ms at the first colum here
hom1 = cleanVarNames(insertTime(hom1)) 
rhod1 = cleanVarNames(insertTime(rhod1))

## Check the tope five rows of the data frame 
print (hom1.head())
print (rhod1.head())

## Visualise a few variables here 
time, trace1, trace2 =  hom1['Time'], hom1['Mean Intensity '], hom1['Mean Intensity .1']
rhod_1, rhod_2 = rhod1['490/10(466.1) Mean Intensity '], rhod1['490/10(466.1) Mean Intensity .1']

# Call the plotting function here
plot_fura2_rhodamine(time, trace1, rhod_1)

## Analyse the calcuim trace signal here 
stim_regions = [(0,200),(250,1200)]
analyze_all_cells(hom1, stim_regions = stim_regions, plot_each = True)
