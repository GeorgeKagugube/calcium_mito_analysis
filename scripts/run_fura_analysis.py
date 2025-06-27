#!/Users/gwk/anaconda3/bin/python3
"""
Created on Thu Oct 19 10:00:00 2023"

Python script for analysing Fura-2 AM calcium imaging data from a dataframe
This script uses the package to perform analysis on Fura-2 AM calcium imaging data.
It reads data from a CSV file, processes it, and generates plots to visualize the results.
The results are saved to a CSV file and the plots are displayed using matplotlib.
"""

import pandas as pd
import numpy as np
import warnings
import os
import sys
from pathlib import Path
from glob import glob

## Import own written functions needed in the analysis of the data
output_dir = Path('/Users/gwk/Desktop/Bioinformatics/Calcium_analysis/data')
file_path = Path('/Users/gwk/Desktop/Bioinformatics/calcium_mito_analysis/data/Platting_2/csv_file')

fura_data, rhodamine_data = [], []

for item in list(file_path.glob('*.csv')):
    if 'Fura' in item.stem or 'fura' in item.stem:
        fura_data.append(item)
    else:
        rhodamine_data.append(item)

from utils.io_utils import load_data, insertTime, cleanVarNames
from utils.plotting import plot_fura2_rhodamine
#from utils.fura_analysis import analyze_single_trace, analyze_all_cells
from utils.fura2_smoothed_scaled_analysis import batch_process_csv, analyze_trace, plot_trace_scaling

## Load the data to be analysed here
hom1 = load_data(sorted(fura_data)[4])
rhod1 = load_data(sorted(rhodamine_data)[4])

## Check the header of the files here
#print (hom1.head())
#print (rhod1.head())

## Insert time in ms at the first colum here
hom1 = cleanVarNames(insertTime(hom1)) 
rhod1 = cleanVarNames(insertTime(rhod1))

## for some datasets, the coverslip moved during imaging,
## In order to remove these artifacts, data in the ranges where this happened were set to the baseline

## Check the tope five rows of the data frame 
#print (hom1)
#print (rhod1.head())

## Select the half of the curv to analyse at this point
#df = hom1[(hom1["Time"] >= 400) & (hom1["Time"] <= 800)]
#df2 = rhod1[(rhod1["Time"] >= 400) & (rhod1["Time"] <= 820)]
#hom, rhod = hom1, rhod1
#hom1.loc[(hom1["Time"] >= 0) & (hom1["Time"] <= 280), hom1.columns[1:]] = 0.6
#hom1.loc[(hom1["Time"] >= 430) & (hom1["Time"] <= 600), hom1.columns[1:]] = 0.67
#print (df.head())
#print (df.shape)
#print (df2.shape)

# Unfiltered 
df = hom1

## Visualise a few variables here 
time, trace1, trace2 =  df['Time'], df['Mean Intensity '], df['Mean Intensity .1']
rhod_1, rhod_2 = df['Mean Intensity '], df['Mean Intensity .1']

# Call the plotting function here
plot_fura2_rhodamine(time, trace1, rhod_1)


## Analyse the calcuim trace signal here 
summary_results = batch_process_csv(df, plot_scaling = True, 
                                auto_window=True,
                               window_size=30,  # seconds
                               step_size=5, output_dir = output_dir)

'''
fdata = pd.DataFrame(summary_results)
#fdata.drop(['SNR','Decay_Tau'], axis=1, inplace=True)
#fdata.dropna(inplace=True)
print (fdata.head(20))

## Export the dataframe for downstrem analysis in r.
## Ask the user to enter a file name here 
#fname = input('Enter the name that you wish to save the file as here: ')
#fdata.to_csv(f'{output_dir}/{fname}.csv')
'''