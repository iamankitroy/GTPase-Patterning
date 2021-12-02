#!/Users/roy/anaconda3/bin/python

import numpy as np
import pandas as pd
import sys


#--- Get input data
def dataIN(filename):
	data = pd.read_csv(filename)
	return data


#--- Get number of bins
def get_numBins(data):
	x_bins = int((max(data["POSITION_X"]) // binsize) + 1)				# bins in x-axis
	y_bins = int((max(data["POSITION_X"]) // binsize) + 1)				# bins in y-axis

	return (x_bins, y_bins)


#--- Bin coordinates of a single frame
def bin_SingleFrame(data):

	# store binned frame data
	binnedFrame = [[[] for yi in range(y_bins)] for xi in range(x_bins)]

	# place coordinates in appropriate bins
	for index, row in data.iterrows():
		x_coord = row["POSITION_X"]						# spot x-coordinate
		y_coord = row["POSITION_Y"]						# spot y-coordinate

		xi = int(x_coord // binsize)					# x bin index
		yi = int(y_coord // binsize)					# y bin index

		binnedFrame[xi][yi].append((x_coord, y_coord))	# store spot in bin

	return binnedFrame


#--- Bin spots from all frames
def SpotBinning(data):
	start_frame = int(min(data["FRAME"]))				# start frame
	end_frame = int(max(data["FRAME"]))					# end frame

	binnedData = {}										# stores binned data from all frames

	# bin spots from every frame
	for frame in range(start_frame, end_frame + 1):
		frameData = data[data["FRAME"] == frame]		# single frame data
		binnedFrame = bin_SingleFrame(frameData)		# binned data from single frame
		binnedData[frame] = binnedFrame					# store single frame binned data

	return binnedData


#--- Get neighbouring bins
def get_Neighbours(xi, yi):

	x_search = [xi+t for t in [1, 0, -1] if (xi+t >= 0) and (xi+t < x_bins)]			# x index of neighbouring bins
	y_search = [yi+t for t in [1, 0, -1] if (yi+t >= 0) and (yi+t < y_bins)]			# y index of neighbouring bins

	# x and y indices of neighbouring bins
	neighbours = [(x, y) for x in x_search for y in y_search if not((x == xi) and (y == yi))]

	return neighbours


#--- Calculate Eucledian distances
def calc_distance(x1, y1, x2, y2):

	d = round(np.sqrt((x1 - x2)**2 + (y1 - y2)**2), 3)		# distance
	return d


#--- Calculate minimum pairwise distance between spots in a single frame
def calc_minPairwiseDistances(cellCoords, neighbourCoords):

	dists = {}

	# empty cells
	if len(cellCoords) == 0:
		return []

	# Calculate distances between points inside the central cell and in neighbouring cells
	for spot1_index in range(len(cellCoords)):
		if spot1_index not in dists.keys():
			dists[spot1_index] = []						# initialise distance list for new spot1 inside central cell

		# distances in the central cell
		for spot2_index in range(spot1_index + 1, len(cellCoords)):
			if spot2_index not in dists.keys():
				dists[spot2_index] = []					# initialise distance list for new spot2 inside central cell

			x1 = cellCoords[spot1_index][0]				# central cell spot1 x coordinate
			y1 = cellCoords[spot1_index][1]				# central cell spot1 y coordinate
			x2 = cellCoords[spot2_index][0]				# central cell spot2 x coordinate
			y2 = cellCoords[spot2_index][1]				# central cell spot2 y coordinate

			d = calc_distance(x1, y1, x2, y2)			# Eucledian distance between central cell spot1 and spot2

			dists[spot1_index].append(d)				# add distance for spot1
			dists[spot2_index].append(d)				# add distance for spot2

		# distances with spots in the neighbouring cells
		for spot2_index in range(len(neighbourCoords)):	
			x1 = cellCoords[spot1_index][0]				# central cell spot1 x coordinate
			y1 = cellCoords[spot1_index][1]				# central cell spot1 y coordinate
			x2 = neighbourCoords[spot2_index][0]		# neighbouring cell spot2 x coordinate
			y2 = neighbourCoords[spot2_index][1]		# neighbouring cell spot2 y coordinate

			d = calc_distance(x1, y1, x2, y2)			# Eucledian distance between central cell spot1 and neighbouring cell spot2

			dists[spot1_index].append(d)				# add distance for spot1

	# store minimum pairwise distances
	min_dists = [min(dists[key]) for key in dists if len(dists[key]) != 0]		# minimum distances for every spot

	return min_dists


#--- Get distances for spots in binned data
def get_Distances(binnedData):

	pairwiseDists = {}					# store pairwise distances

	# Get pairwise distances from spots in a binned search space
	for frame in binnedData:
		binnedFrame = binnedData[frame]						# single frame with binned spots
		pairwiseDists[frame] = []							# initiate empty list for frame

		# get search space for every bin
		# search space includes central bin and the adjacent bins
		for xi in range(len(binnedFrame)):
			for yi in range(len(binnedFrame[xi])):
				cellCoords = binnedFrame[xi][yi]			# central cell coordinates
				neighbours = get_Neighbours(xi, yi)			# search space
				neighbourCoords = []						# stores all coordinates from search space

				# store all coordinates from search space
				for bin_coord in neighbours:
					neighbourCoords.extend(binnedFrame[bin_coord[0]][bin_coord[1]])		# coordinates for neighbouring spots

				# calculate and store minimum pairwise distances for all spots in central cell
				pairwiseDists[frame].extend(calc_minPairwiseDistances(cellCoords, neighbourCoords))

	return pairwiseDists


#--- Write data
def dataOUT(filename, distData):

	outname = filename[:-4]
	outname = f"{outname}_dist.csv"			# output file name

	frames = []								# stores all frames
	dists = []								# stores all distances

	# split data into frames and distance lists
	for frame in distData:
		for dist in distData[frame]:
			frames.append(frame)			# add frame
			dists.append(dist)				# add distance

	# data frame with columns to store frames (FRAME) and distance data (DIST)
	distData = pd.DataFrame({"FRAME": frames, "DIST": dists})

	# write output file
	distData.to_csv(outname, index=False, float_format="%.3f")


#--- Main function
def main():
	filename = sys.argv[1]					# filename

	global x_bins, y_bins, binsize

	binsize = 5								# binsize in µm

	data = dataIN(filename)					# spot data
	x_bins, y_bins = get_numBins(data)		# number of bins in x axis and y axis

	binnedData = SpotBinning(data)			# binned spot data

	pairwiseDists = get_Distances(binnedData)		# minimum pairwise distances for all spots

	dataOUT(filename, pairwiseDists)		# write data out

main()

# Ankit Roy
# 1st December, 2021
