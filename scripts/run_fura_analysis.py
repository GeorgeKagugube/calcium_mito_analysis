#!/Users/gwk/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 10:00:00 2023"
""
""
"Starter Python notebook for analysing Fura-2 AM calcium imaging data from a dataframe
This script uses the FuraAnalysis class to perform analysis on Fura-2 AM calcium imaging data.
It reads data from a CSV file, processes it, and generates plots to visualize the results.
The results are saved to a CSV file and the plots are displayed using matplotlib.
"""
# Starter Python notebook for analysing Fura-2 AM calcium imaging data from a dataframe
# Import necessary libraries
import os
import ../utils.fura_analysis as fa
import ../utils.io_
# Ensure that the 'results' directory exists
if not os.path.exists('results'):
    os.makedirs('results')  
# Import the FuraAnalysis and FuraAnalysisConfig classes from the fura_analysis module
# Adjust the path to the fura_analysis module as necessary
# This assumes that the fura_analysis module is in a directory named 'utils' relative to this script
# If the module is in a different location, adjust the import statement accordingly
# Import the FuraAnalysis and FuraAnalysisConfig classes from the fura_analysis module
# Import necessary libraries

import sys
# Add the path to the fura_analysis module
sys.path.append('/Users/gwk/Desktop/Bioinformatics/calcium_mito_analysis/utils')  # Adjust this path
# Import the FuraAnalysis and FuraAnalysisConfig classes from the fura_analysis module

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from scipy.integrate import simpson
from scipy.optimize import curve_fit, OptimizeWarning
import seaborn as sns
import warnings


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