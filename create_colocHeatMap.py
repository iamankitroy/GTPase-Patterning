#!/Users/roy/anaconda3/bin/python

import sys
import pandas as pd

#--- Get data
def dataIN(filename):
	data = pd.read_csv(filename,
		comment="#",
		dtype = {
 				"TRACK_ID" : str,
 				"PSEUDO_TRACK_ID" : str,
 				"COLOCALIZED_SPOT" : bool,
 				"COLOCALIZATION_ID": str,
 				"CHANNEL" : str
 				})

	return data

#--- Get start and end frame for all tracks
def get_frame_limits(data):

	# group by tracks from different channels
	groupings = data.groupby(["PSEUDO_TRACK_ID", "CHANNEL"])

	all_track_starts = {}							# store all starting frames
	all_track_ends = {}								# store all ending frames

	# get start and end frame for every group
	for gid, group in groupings:
		track_start = min(group["FRAME"])			# track starting frame
		track_end = max(group["FRAME"])				# track ending frame

		all_track_starts[gid] = track_start			# store track starting
		all_track_ends[gid] = track_end				# store track ending

	return all_track_starts, all_track_ends

#--- Normalize frame start and end
def normalize_frames(data, all_track_starts, all_track_ends):

	# normalize start position
	data["NORM_START"] = data.apply(lambda x: x["FRAME"] - all_track_starts[(x["PSEUDO_TRACK_ID"], x["CHANNEL"])], axis=1)
	# normalize end position
	data["NORM_END"] = data.apply(lambda x: x["FRAME"] - all_track_ends[(x["PSEUDO_TRACK_ID"], x["CHANNEL"])], axis=1)
	# track duration
	data["TRACK_LENGTH"] = data.apply(lambda x: all_track_ends[(x["PSEUDO_TRACK_ID"], x["CHANNEL"])] - all_track_starts[(x["PSEUDO_TRACK_ID"], x["CHANNEL"])] + 1, axis=1)
	return data

#--- Generate plotting data
def gen_plotData(data):

	data["COLOCALIZATION_STATUS"] = "None"
	data["COLOCALIZATION_STATUS"] = data.apply(lambda x: "Internal" if x["COLOCALIZED_SPOT"] else x["COLOCALIZATION_STATUS"], axis=1)
	data["COLOCALIZATION_STATUS"] = data.apply(lambda x: "Recruitment" if x["ANNOTATION_SPOT"] == "Recruitment" else x["COLOCALIZATION_STATUS"], axis=1)
	data["COLOCALIZATION_STATUS"] = data.apply(lambda x: "Extraction" if x["ANNOTATION_SPOT"] == "Extraction" else x["COLOCALIZATION_STATUS"], axis=1)

	plotdata = data[["PSEUDO_TRACK_ID", "CHANNEL", "FRAME", "NORM_START", "NORM_END", "COLOCALIZATION_STATUS", "TRACK_LENGTH", "ANNOTATION_TRACK"]]

	return plotdata

#--- Write plot data
def writeOUT(plotdata, filename):

	# output file name
	outname = filename[:-4]
	outname = f"{outname}_heatPlotData.csv"

	# write output file
	plotdata.to_csv(outname, index=False)

#--- Main function
def main():
	filename = sys.argv[1]				# input file name
	data = dataIN(filename)				# colocalization data

	all_track_starts, all_track_ends = get_frame_limits(data)				# store track start and end positions
	data = normalize_frames(data, all_track_starts, all_track_ends)			# normalize frames to start and end positions
	plotdata = gen_plotData(data)											# generate plotable data
	writeOUT(plotdata, filename)											# write plot file
	

#--- Run main function
main()

# Ankit Roy
# 15th February, 2022
