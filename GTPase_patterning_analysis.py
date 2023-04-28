#!/Users/roy/mambaforge/bin/python

import sys
import numpy as np
import pandas as pd
from skimage import img_as_float, img_as_ubyte, io

def readImg(filename):
    img = img_as_float(io.imread(filename))
    return img

def calcPatterningRatio(img, segmented_img):
    ratio = [(np.sum(img[t]*segmented_img[t])/np.sum(segmented_img[t]),         # average intensity in patterned regions
            np.sum(img[t]*~segmented_img[t]/np.sum(~segmented_img[t])),         # average intensity outside patterened regions
            (np.sum(img[t]*segmented_img[t])/np.sum(segmented_img[t]))/(np.sum(img[t]*~segmented_img[t]/np.sum(~segmented_img[t]))))        # patterning ratio
            for t in range(len(segmented_img))]
    return ratio

def dataOUT(data, filename):
    outname = '-'.join(filename.split('_C=')[:-1])
    outname = f'{outname}-analyzed.csv'

    data.to_csv(outname, index=False)

def main():
    segmented_file = sys.argv[1]
    lipid_file = sys.argv[2]
    total_GTPase_file = sys.argv[3]
    active_GTPase_file = sys.argv[4]

    segmentation_map = np.bool_(readImg(segmented_file))
    lipid_img = readImg(lipid_file)
    total_GTPase_img = readImg(total_GTPase_file)
    active_GTPase_img = readImg(active_GTPase_file)

    lipid_pattern = calcPatterningRatio(lipid_img, segmentation_map)
    total_GTPase_pattern = calcPatterningRatio(total_GTPase_img, segmentation_map)
    active_GTPase_pattern = calcPatterningRatio(active_GTPase_img, segmentation_map)

    analyzed_data = pd.DataFrame({
        'frame' : [n for n in range(len(segmentation_map))],
        'lipid_avgI_in' : [n[0] for n in lipid_pattern],
        'totalGTPase_avgI_in' : [n[0] for n in total_GTPase_pattern],
        'activeGTPase_avgI_in' : [n[0] for n in active_GTPase_pattern],
        'lipid_avgI_out' : [n[1] for n in lipid_pattern],
        'totalGTPase_avgI_out' : [n[1] for n in total_GTPase_pattern],
        'activeGTPase_avgI_out' : [n[1] for n in active_GTPase_pattern],
        'lipid_PI' : [n[2] for n in lipid_pattern],
        'totalGTPase_PI' : [n[2] for n in total_GTPase_pattern],
        'activeGTPase_PI' : [n[2] for n in active_GTPase_pattern],
    })

    analyzed_data['activeGTPase_norm_in'] = analyzed_data['activeGTPase_avgI_in']/analyzed_data['totalGTPase_avgI_in']
    analyzed_data['activeGTPase_norm_out'] = analyzed_data['activeGTPase_avgI_out']/analyzed_data['totalGTPase_avgI_out']
    analyzed_data['activeGTPase_enrichment_in'] = analyzed_data['activeGTPase_norm_in']/analyzed_data['activeGTPase_norm_out']

    dataOUT(analyzed_data, segmented_file)

main()

# Ankit Roy
# 9th January, 2023
# 28th April, 2023          >>      Changed deprecated bool8 function to bool_