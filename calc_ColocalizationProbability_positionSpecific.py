#!/Users/roy/anaconda3/bin/python

import sys
import pandas as pd

#--- Get data
def dataIN(filename):
	data = pd.read_csv(filename, comment="#",
			dtype = {
				"TRACK_ID" : str,
				"PSEUDO_TRACK_ID" : str,
				"COLOCALIZED_SPOT" : bool,
				"COLOCALIZATION_ID": str,
				"CHANNEL" : str
				})
	return data

#--- Get single channel data
def get_singleChannel(data, channel="GTPase"):
	data = data[data["CHANNEL"] == channel]
	return data

#--- Normalize track start
def normalize_trackStart(data):

	# stores data for all track starting positions
	all_trackStarts = {}

	# group tracks according to pseudo track id
	groupings = data.groupby(["PSEUDO_TRACK_ID"])
	
	# find and store track start for all tracks
	for gid, group in groupings:
		trackStart = min(group["FRAME"])		# track start frame
		all_trackStarts[gid] = trackStart		# store track start

	data["NORM_FRAME"] = data.apply(lambda x: x["FRAME"] - all_trackStarts[x["PSEUDO_TRACK_ID"]], axis=1)

	return data

#--- Get position specific colocalization probabilities
def get_posProbs(data):

	# stores all position specific colocalization probabilities and observation counts
	# stores data as a tuple with the following indices
	# 0. Position specific colocalization probability
	# 1. Observation count from which probability was calculated
	posProbs = {}

	# group data according to the normalized frame
	groupings = data.groupby(["NORM_FRAME"])

	# calculate the colocalization probability at every position and store along with the number of observations at that position
	for gid, group in groupings:
		nobs = len(group)										# number of observations at specific normalized frame
		prob = round(group["COLOCALIZED_SPOT"].sum()/nobs, 4)	# colocalization probability
#		print(f"{gid},{prob},{nobs}")

		posProbs[gid] = (prob, nobs)		# store number of observations and colocalization probability

	return posProbs

#--- Generate plotabble data
def gen_plotOut(posProbs, time_resolution):

	frames = sorted(list(posProbs.keys()))			# normalized frames
	probs = [posProbs[key][0] for key in frames]	# position specific colocalization probabilities
	nobs = [posProbs[key][1] for key in frames]		# number of observations

	# plottable data frame
	posProbs_df = pd.DataFrame({
		"FRAME": frames,
		"PROBS": probs,
		"NOBS": nobs
	})

	# calculate lifetime
	posProbs_df["LIFETIME"] = posProbs_df["FRAME"] * time_resolution

	return posProbs_df

#--- Generate plot file
def gen_plotFile(posProbs, filename):

	outname = filename[:-4]
	outname = f"{outname}_posProbPlot.csv"		# output file name

	# write plot file
	posProbs.to_csv(outname, index=False, float_format="%.4f")

#--- Main function
def main():
	filename = sys.argv[1]				# file name

	time_resolution = 0.022				# s

	data = dataIN(filename)				# load single molecule data
	pd.set_option('display.max_columns', None)

	data = get_singleChannel(data)		# get single channel

	data = normalize_trackStart(data)	# normalize track start positions

	posProbs = get_posProbs(data)		# positional colocalization probability and observation counts

	posProbs = gen_plotOut(posProbs, time_resolution)	# convert to plottable data frame

	gen_plotFile(posProbs, filename)	# generate plot file

#--- Run main
main()

# Ankit Roy
# 2nd February, 2022
# 10th February, 2022	--> Explicitly states the data types for certain columns of input file
# 9th February, 2024	--> Now calculates lifetime from time resolution