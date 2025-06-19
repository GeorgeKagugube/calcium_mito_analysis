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
output_dir = Path('/Users/gwk/Desktop/Bioinformatics/calcium_mito_analysis/results')
file_path = Path('/Users/gwk/Desktop/Bioinformatics/calcium_mito_analysis/data/Platting_1/csv_files')

fura_data, rhodamine_data = [], []

for item in list(file_path.glob('*.csv')):
    if 'Fura' in item.stem or 'fura' in item.stem:
        fura_data.append(item)
    else:
        rhodamine_data.append(item)

from utils.io_utils import load_data, insertTime, cleanVarNames
from utils.plotting import plot_fura2_rhodamine
from utils.fura_analysis import analyze_single_trace, analyze_all_cells

## Load the data to be analysed here
hom1 = load_data(sorted(fura_data)[1])
rhod1 = load_data(sorted(rhodamine_data)[1])

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


#hom, rhod = hom1, rhod1
#hom.loc[(hom["Time"] >= 0) & (hom["Time"] <= 280), hom.columns[1:]] = 0.6
#hom.loc[(hom["Time"] >= 430) & (hom["Time"] <= 600), hom.columns[1:]] = 0.67

## Visualise a few variables here 
time, trace1, trace2 =  hom1['Time'], hom1['Mean Intensity '], hom1['Mean Intensity .1']
rhod_1, rhod_2 = rhod1['490/10(466.1) Mean Intensity '], rhod1['490/10(466.1) Mean Intensity .1']

# Call the plotting function here
plot_fura2_rhodamine(time, trace1, rhod_1)


## Analyse the calcuim trace signal here 
stim_regions = [(350,850)]#,(600,800)]#,(550,720),(750, 840)]
summary_results = analyze_all_cells(hom1, stim_regions = stim_regions, plot_each = False)
fdata = pd.DataFrame(summary_results)
fdata.drop(['SNR','Decay_Tau'], axis=1, inplace=True)
fdata.dropna(inplace=True)

print(fdata.head(50))

## Export the dataframe for downstrem analysis in r.
## Ask the user to enter a file name here 
fname = input('Enter the name that you wish to save the file as here: ')
fdata.to_csv(f'{output_dir}/{fname}.csv')
