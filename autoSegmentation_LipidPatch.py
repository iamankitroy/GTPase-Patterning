#!/Users/roy/mambaforge/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import cv2
from scipy import ndimage as nd
from skimage import img_as_float, img_as_ubyte, io
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.filters import threshold_otsu
from skimage.segmentation import random_walker

def readImg(filename):
    img = img_as_float(io.imread(filename))
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

def segment_RandomWalker(img, modifier=0.2):
    threshold = threshold_otsu(img)

    markers = np.zeros(img.shape, dtype=np.uint)
    markers[(img >= (threshold * (1+modifier)))] = 2
#    markers[(img >= (threshold * (1)))] = 2
#    markers[(img < (threshold * (1-modifier)))] = 1
    markers[(img < (threshold * (1)))] = 1

    labels = random_walker(img, markers, beta=10, mode='bf')
    
    return labels

def cleanUp(img, kernel_size=(9,9), rounds=5):
    for run in range(rounds):
        img = nd.binary_closing(img, np.ones((3,3)))
        img = nd.binary_opening(img, np.ones(kernel_size))
        
    return img

def dilate(img, rounds=1):
    for run in range(rounds):
        img = nd.binary_dilation(img)
        
    return img

def imageSave(filename, suffix, img):
    img = np.bool8(img)

    outname = '.'.join(filename.split('.')[:-1])
    outname = f'{outname}_{suffix}.tif'

    io.imsave(outname, img)

def main():
    filename = sys.argv[1]
    img = readImg(filename)
    segmented_img = np.zeros(np.shape(img))

    for t in range(len(img)):
        denoise = denoiseImg(img[t])
        equalized = runCLAHE(denoise, 2)
        labels = segment_RandomWalker(equalized, 0.2)
        segment = labels == 2
        cleaned_seg = cleanUp(segment)
        dilated_seg = dilate(cleaned_seg, 3)
        segmented_img[t] = dilated_seg
        print(f"{t+1:>5d} of {len(img)}")

    imageSave(filename, 'segmented', segmented_img)

main()

# Ankit Roy
# 9th January, 2023
