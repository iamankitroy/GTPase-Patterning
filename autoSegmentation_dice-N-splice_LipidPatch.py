#!/Users/roy/mambaforge/bin/python

import sys
import numpy as np
from skimage import io, img_as_float, img_as_ubyte, filters
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.segmentation import random_walker
from scipy import ndimage as nd
import cv2
import matplotlib.pyplot as plt
import pandas as pd

def getImage(filename):
    img = io.imread(filename)
    return img

def denoiseImg(img):
    sigma_est = np.mean(estimate_sigma(img, channel_axis=None))
    denoise = denoise_nl_means(img, h=1.15 * sigma_est, fast_mode=False, patch_size=5, patch_distance=6)
    return denoise

def runCLAHE(img, size=4):
    img = img_as_ubyte(img)
    clahe = cv2.createCLAHE(clipLimit = np.max(img), tileGridSize = (size,size))
    equalized = clahe.apply(img)
    return equalized

def segment_RandomWalker(img, upper_percentile=70, lower_percentile=30):

    markers = np.zeros(img.shape, dtype=np.uint)
    markers[(img >= np.percentile(img, upper_percentile))] = 2
    markers[(img <= np.percentile(img, lower_percentile))] = 1

    labels = random_walker(img, markers, beta=10, mode='bf')
    
    return labels

def cleanUp(img, kernel_size=(9,9), rounds=5):
    for run in range(rounds):
        img = nd.binary_closing(img, np.ones((5,5)))
        img = nd.binary_opening(img, np.ones(kernel_size))
        
    return img

def dilate(img, rounds=1):
    for run in range(rounds):
        img = nd.binary_dilation(img)
        img = nd.binary_closing(img, np.ones((5,5)))
        img = nd.binary_opening(img, np.ones((5,5)))
        
    return img

def imageSave(filename, suffix, img):
    img = np.bool8(img)

    outname = '.'.join(filename.split('.')[:-1])
    outname = f'{outname}_{suffix}.tif'

    io.imsave(outname, img)

def main():
    filename = sys.argv[1]
    img = getImage(filename)

    slices = np.linspace(0, img.shape[1], 3, dtype=int)

    segmented_img = np.zeros(img.shape)

    for t in range(len(img)):
        spliced_img = np.zeros(img[t].shape, dtype=bool)

        for x_index in range(len(slices)-1):
            for y_index in range(len(slices)-1):
                
                denoise = denoiseImg(img[t,slices[x_index]:slices[x_index+1]:,slices[y_index]:slices[y_index+1]:])
                equalized = runCLAHE(denoise, 2)
                upper_thr, lower_thr = (50, 20)
                labels = segment_RandomWalker(equalized, upper_thr, lower_thr)
                segment = labels == 2
                spliced_img[slices[x_index]:slices[x_index+1], slices[y_index]:slices[y_index+1]] = segment

        cleaned_seg = cleanUp(spliced_img)
        # dilated_seg = dilate(cleaned_seg, 2)
        segmented_img[t] = cleaned_seg
        print(f"{t+1:>5d} of {len(img)}")

    imageSave(filename, 'segmented', segmented_img)

main()

# Ankit Roy
# 7th June, 2024