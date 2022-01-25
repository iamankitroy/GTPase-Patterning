#!/Users/roy/anaconda3/bin/python

import pandas as pd
import sys
 
__author__ = "Ankit Roy"
__copyright__ = "Copyright 2022, Bieling Lab, Max Planck Institute of Molecular Physiology"
__license__ = "GPL"
__maintainer__ = "Ankit Roy"
__status__ = "Development"

# Get data
def dataIN(filename):
	data = pd.read_csv(filename, comment='#')
	return data

# Get single channel data
# Default: GTPase channel
def singleChannel(data):
	data = data[data["CHANNEL"] == channel]
	return data

# Get recruitment frames
def get_recruitmentFrames(data):
	max_rec_frame = min(data["FRAME"]) + frame_threshold			# max recruitment frame number
	subset_data = data.loc[(data["FRAME"] < max_rec_frame), ]		# recruitment frame data

	recruitment_events = sum(subset_data["COLOCALIZED_SPOT"])		# recruitment events
	total_frames = len(subset_data.index)							# total number of recruitment frames

	return (recruitment_events, total_frames, round(recruitment_events/frame_threshold, 2))

# Get extraction frames
def get_extractionFrames(data):
	min_ext_frame = max(data["FRAME"]) - frame_threshold			# min extraction frame number
	subset_data = data.loc[(data["FRAME"] > min_ext_frame), ]		# extraction frame data

	extraction_events = sum(subset_data["COLOCALIZED_SPOT"])		# extraction events
	total_frames = len(subset_data.index)							# total number of extraction frames

	return (extraction_events, total_frames, round(extraction_events/frame_threshold, 2))

# Get internal frames
def get_internalFrames(data):
	min_int_frame = min(data["FRAME"]) + frame_threshold			# min internal frame number
	max_int_frame = max(data["FRAME"]) - frame_threshold			# max internal frame number
	subset_data = data.loc[(data["FRAME"] >= min_int_frame) & (data["FRAME"] <= max_int_frame), ]	# internal frame data

	internal_events = sum(subset_data["COLOCALIZED_SPOT"])			# internal colocalization events
	total_frames = len(subset_data.index)							# total number of internal frames

	# Return NA if no internal frames are found 
	if total_frames != 0:
		return (internal_events, total_frames, round(internal_events/total_frames, 2))
	else:
		return (internal_events, total_frames, 'NA')

# Get recruitment frames
# Default: First 3 frames
def classifyFrames(data):

	# group by pseudo track ids
	groupings = data.groupby(["PSEUDO_TRACK_ID"])

	# recruitment, extraction and internal colocalization statistics
	# for each pseudo track id stores a tuple with 3 values:
	# 	1. recruitment/extraction/internal colocalization event count
	# 	2. total frames qualifying class criteria
	#	3. event probability
	recruitmentStats = {}
	extractionStats = {}
	internalStats = {}

	# store event probabilities alone
	recruitmentProbs = []
	extractionProbs = []
	internalProbs = []

	# 
	for gid, group in groupings:
		min_frame = min(group["FRAME"])				# first frame
		max_frame = max(group["FRAME"])				# last frame
		track_length = max_frame - min_frame + 1	# track length

		# Skip smaller tracks
		if track_length < min_track_length:
			continue

		recruitmentStats[gid] = get_recruitmentFrames(group)		# recruitment statistics
		extractionStats[gid] = get_extractionFrames(group)			# extraction statistics
		internalStats[gid] = get_internalFrames(group)				# internal colocalization statistics

		recruitmentProbs.append(recruitmentStats[gid][2])			# store recruitment probabilities
		extractionProbs.append(extractionStats[gid][2])				# store extraction probabilities

		# store internal colocalization probability if present
		if internalStats[gid][2] != 'NA':
			internalProbs.append(internalStats[gid][2])

	return recruitmentProbs, extractionProbs, internalProbs


# Main function
def main():

	global frame_threshold
	global min_track_length
	global channel

	frame_threshold = 3				# frames to consider for recruitment and extraction
	min_track_length = 5			# minimum track length to consider for analysis
	channel = "GTPase"				# channel for analysis

	pd.set_option('display.max_columns', None)
	filename = sys.argv[1]			# input colocalization subset file name
	data = dataIN(filename)			# colocalization subset data

	data = singleChannel(data)		# GTPase channel data

	# classify frames into recruitment, extraction or internal and calculate probabilities
	recruitmentProbs, extractionProbs, internalProbs = classifyFrames(data)


# Run main function
main()

# Ankit Roy
# 21st January, 2022
# 25th January, 2022		>>		Classifies tracks into recruitment, extraction and internal frames and returns event probabilities
