#!/Users/roy/mambaforge/bin/python

import os
import sys
import numpy as np
import pickle
import glob
from skimage import img_as_float, io

# Read image file
def readImg(filename):
    img = img_as_float(io.imread(filename))
    return img

# Calculate mean intensity for background in a given channel
def mean_channelBackground(background_dir, channel=""):

    path_regex = os.path.join(background_dir, f"*{channel}.tif")
    imgs = np.array([readImg(file) for file in glob.glob(path_regex)])

    return imgs.mean()

# Compute backgrounds for all channels
def compute_channelBackground(background_dir, channels):

    ch_backgrounds = {ch: mean_channelBackground(background_dir, ch) for ch in channels}
    return ch_backgrounds

# Pickle output file
def dataOUT(backgrounds, filename="background_dict.p"):
    with open(filename, "wb") as fh:
        pickle.dump(backgrounds, fh)

# Main function
def main():
    background_dir = sys.argv[1]                # path to images for background calculation
    channels = ["Lipid", "Total", "Active"]     # channel names

    # store background values for all channels
    ch_backgrounds = compute_channelBackground(background_dir, channels)
    
    # pickle dictionary with background values
    dataOUT(ch_backgrounds)
    
# Run main
main()

# Ankit Roy
# 19th June, 2023