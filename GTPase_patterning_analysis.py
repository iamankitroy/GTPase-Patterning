#!/Users/roy/mambaforge/bin/python

import sys
import pickle
import numpy as np
import pandas as pd
from skimage import img_as_float, img_as_ubyte, io

# Import images
def readImg(filename):
    img = img_as_float(io.imread(filename))
    return img

# Import background data
def getBackground(background_file):
    with open(background_file, "rb") as fh:
        ch_background = pickle.load(fh)
    return ch_background

# Calculate patterning index
def calcPatterningRatio(img, segmented_img, channel, background):

    total_t = len(segmented_img)            # total time in timeseries

    # store pattern data
    patternData = {
        'frame' : np.zeros(total_t),
        'avgI_in' : np.zeros(total_t),
        'avgI_out' : np.zeros(total_t),
    }

    # calculate pattern statistics
    for t in range(total_t):
        # current frame
        patternData['frame'][t] = t
        # average instensity inside segmented regions
        patternData['avgI_in'][t] = np.sum(img[t]*segmented_img[t]) / np.sum(segmented_img[t])
        # average instensity outside segmented regions
        patternData['avgI_out'][t] = np.sum(img[t]*~segmented_img[t]) / np.sum(~segmented_img[t])

    # create data frame from pattern data
    patternData = pd.DataFrame.from_dict(patternData)
    
    # caculate Patterning Index
    patternData['PI'] = patternData['avgI_in'] / patternData['avgI_out']

    # caculate pattern statistics with background subtraction
    patternData['avgI_in_BS'] = patternData['avgI_in'] - background[channel]
    patternData['avgI_out_BS'] = patternData['avgI_out'] - background[channel]
    patternData['PI_BS'] = patternData['avgI_in_BS'] / patternData['avgI_out_BS']
    
    # set frame as index
    patternData = patternData.set_index('frame')
    # rename columns with channel prefix
    col_names = {col: f"{channel}_{col}" for col in patternData.columns}
    patternData.rename(columns=col_names, inplace=True)

    return patternData

# Normalize activity channel for enrichment calculations
def preEnrichment_normalization(data, active, total):
    return data[active]/data[total]

# Calculate enrichment
def calc_enrichment(data, avgI_in, avgI_out):
    return data[avgI_in]/data[avgI_out]

# Write output files
def dataOUT(data, filename):
    outname = '-'.join(filename.split('_C=')[:-1])
    outname = f'{outname}-analyzed.csv'

    data.to_csv(outname, index=True)

# Main function
def main():
    segmented_file = sys.argv[1]            # segmentation file
    lipid_file = sys.argv[2]                # lipid channel file
    total_GTPase_file = sys.argv[3]         # total GTPase channel file
    active_GTPase_file = sys.argv[4]        # active GTPase channel file
    background_file = sys.argv[5]           # background data file

    segmentation_map = np.bool_(readImg(segmented_file))    # segmentation map
    lipid_img = readImg(lipid_file)                         # lipid channel image data
    total_GTPase_img = readImg(total_GTPase_file)           # total GTPase channel image data
    active_GTPase_img = readImg(active_GTPase_file)         # active GTPase channel image data
    ch_background = getBackground(background_file)          # background data

    # Pattern statistics for the lipid channel
    lipid_pattern = calcPatterningRatio(
        lipid_img,
        segmentation_map,
        "Lipid",
        ch_background
        )
    # Pattern statistics for the total GTPase channel
    total_GTPase_pattern = calcPatterningRatio(
        total_GTPase_img,
        segmentation_map,
        "Total",
        ch_background
        )
    # Pattern statistics for the active GTPase channel
    active_GTPase_pattern = calcPatterningRatio(
        active_GTPase_img,
        segmentation_map,
        "Active",
        ch_background
        )

    # Combine pattern statistics from all channels
    analyzed_data = pd.concat(
        [lipid_pattern, total_GTPase_pattern, active_GTPase_pattern],
        axis=1
    )

    # Normalize active GTPase channel for enrichment calculation
    analyzed_data['Active_Norm_in'] = preEnrichment_normalization(analyzed_data, "Active_avgI_in", "Total_avgI_in")
    analyzed_data['Active_Norm_out'] = preEnrichment_normalization(analyzed_data, "Active_avgI_out", "Total_avgI_out")
    analyzed_data['Active_Norm_in_BS'] = preEnrichment_normalization(analyzed_data, "Active_avgI_in_BS", "Total_avgI_in_BS")
    analyzed_data['Active_Norm_out_BS'] = preEnrichment_normalization(analyzed_data, "Active_avgI_out_BS", "Total_avgI_out_BS")

    # Calculate enrichment
    analyzed_data['Active_Enrichment_in'] = calc_enrichment(analyzed_data, "Active_Norm_in", "Active_Norm_out")
    analyzed_data['Active_Enrichment_in_BS'] = calc_enrichment(analyzed_data, "Active_Norm_in_BS", "Active_Norm_out_BS")

    # Write output files
    dataOUT(analyzed_data, segmented_file)

main()

# Ankit Roy
# 9th January, 2023
# 28th April, 2023          >>      Changed deprecated bool8 function to bool_
# 20th June, 2023           >>      Revamped entire script to be more Pythonic.
#                           >>      Added calculations for all metrics with background subtraction.
#                           >>      Uses background data from a pre-generated pickle file.
#                           >>      Added comments.