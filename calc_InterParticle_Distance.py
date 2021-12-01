#!/Users/roy/anaconda3/bin/python

import numpy as np
import pandas as pd
import sys


def dataIN(filename):
	data = pd.read_csv(filename)
	return data


def get_numBins(data):
	x_bins = int((max(data["POSITION_X"]) // binsize) + 1)
	y_bins = int((max(data["POSITION_X"]) // binsize) + 1)

	return (x_bins, y_bins)


def bin_SingleFrame(data):

	binnedFrame = [[[] for yi in range(y_bins)] for xi in range(x_bins)]

	for index, row in data.iterrows():
		x_coord = row["POSITION_X"]
		y_coord = row["POSITION_Y"]

		xi = int(x_coord // binsize)
		yi = int(y_coord // binsize)

		binnedFrame[xi][yi].append((x_coord, y_coord))

	return binnedFrame


def SpotBinning(data):
	start_frame = int(min(data["FRAME"]))
	end_frame = int(max(data["FRAME"]))

	binnedData = {}

	for frame in range(start_frame, end_frame + 1):
		frameData = data[data["FRAME"] == frame]
		binnedFrame = bin_SingleFrame(frameData)
		binnedData[frame] = binnedFrame

	return binnedData


def get_Neighbours(xi, yi):

	x_search = [xi+t for t in [1, 0, -1] if (xi+t >= 0) and (xi+t < x_bins)]
	y_search = [yi+t for t in [1, 0, -1] if (yi+t >= 0) and (yi+t < y_bins)]

	neighbours = [(x, y) for x in x_search for y in y_search if not((x == xi) and (y == yi))]

	return neighbours


def calc_distance(x1, y1, x2, y2):

	d = round(np.sqrt((x1 - x2)**2 + (y1 - y2)**2), 3)
	return d


def calc_minPairwiseDistances(cellCoords, neighbourCoords):

	dists = {}

	if len(cellCoords) == 0:
		return []

	# Calculate distances between points inside the central cell and in neighbouring cells
	for spot1_index in range(len(cellCoords)):
		if spot1_index not in dists.keys():
			dists[spot1_index] = []

		# distances in the central cell
		for spot2_index in range(spot1_index + 1, len(cellCoords)):
			if spot2_index not in dists.keys():
				dists[spot2_index] = []

			x1 = cellCoords[spot1_index][0]
			y1 = cellCoords[spot1_index][1]
			x2 = cellCoords[spot2_index][0]
			y2 = cellCoords[spot2_index][1]

			d = calc_distance(x1, y1, x2, y2)

			dists[spot1_index].append(d)
			dists[spot2_index].append(d)

		# distances with spots in the neighbouring cells
		for spot2_index in range(len(neighbourCoords)):	
			x1 = cellCoords[spot1_index][0]
			y1 = cellCoords[spot1_index][1]
			x2 = neighbourCoords[spot2_index][0]
			y2 = neighbourCoords[spot2_index][1]

			d = calc_distance(x1, y1, x2, y2)

			dists[spot1_index].append(d)

	# store minimum pairwise distances
	min_dists = [min(dists[key]) for key in dists if len(dists[key]) != 0]

	return min_dists


def get_Distances(binnedData):

	pairwiseDists = {}

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
					neighbourCoords.extend(binnedFrame[bin_coord[0]][bin_coord[1]])

				pairwiseDists[frame].extend(calc_minPairwiseDistances(cellCoords, neighbourCoords))

	return pairwiseDists


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

	distData = pd.DataFrame({"FRAME": frames, "DIST": dists})

	# write output file
	distData.to_csv(outname, index=False, float_format="%.3f")


def main():
	filename = sys.argv[1]

	global x_bins, y_bins, binsize

	binsize = 5				# in µm

	data = dataIN(filename)
	x_bins, y_bins = get_numBins(data)

	binnedData = SpotBinning(data)

	pairwiseDists = get_Distances(binnedData)

	dataOUT(filename, pairwiseDists)

main()

# Ankit Roy
# 1st December, 2021
