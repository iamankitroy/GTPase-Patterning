#!/Users/roy/anaconda3/bin/python3

"""
Identifies colocalized spots in dual channel single molecule data.

Uses spot statistics obtained by tracking single molecule spots with TrackMate (All Spots statistics file).
Takes in two input spot statistics files for the GTPase channel and the GDI channel respectively.
Colocalization events are defined as instances where spots from different channels appear within a threshold distance.
The output file contains a field ("COLOCALIZED_SPOT") which indicates whether or not a spot was colocalized.
The output file has data combined from both the channels and uses the "CHANNEL" field to descrimate which channel the spots belong to.
These fields indicate the colocalization status of a spot and the spot's channel respectively.
Additional fields generated by this script and their brief description:
--> PSEUDO_TRACK_ID: contains <TRACK_ID> if present or <Label> otherwise
--> COLOCALIZED_SPOT: <True> if spot was colocalized <False> otherwise
--> COLOCALIZATION_ID: <GTPase-PSEUDO_TRACK_ID>-<GDI-PSEUDO_TRACK_ID>
--> CHANNEL: <GTPase> if spot belongs to GTPase channel <GDI> otherwise
"""

import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp

__author__ = "Ankit Roy"
__copyright__ = "Copyright 2021, Bieling Lab, Max Planck Institute of Molecular Physiology"
__license__ = "GPL"
__maintainer__ = "Ankit Roy"
__status__ = "Development"

#--- Fetch arguments
def get_args():
	parser = argparse.ArgumentParser()

	# Required arguments group
	required_args = parser.add_argument_group(title = "Required arguments")

	# GTPase channel file name
	required_args.add_argument("-gp", "--gtpase",
								help = "GTPase channel spot statistics file",
								required = True)

	# GDI channel file name
	required_args.add_argument("-gd", "--gdi",
								help = "GDI channel spot statistics file",
								required = True)

	# Distance cutoff
	parser.add_argument("-d", "--dist",
								help = "(default = 0.5 µm) Colocalization distance cutoff",
								default = 0.5,
								type = float)

	# Field of view cutoff
	parser.add_argument("-fov", "--field",
								help = "(default = 1.00) Field of view",
								default = 1.00,
								type = float)

	# Pixel size
	parser.add_argument("-ps", "--pixel_size",
								help = "(default = 0.178 µm) Pixel size",
								default = 0.178,
								type = float)

	# Image size
	parser.add_argument("-is", "--image_size",
								help = "(default = 512 px) Image size",
								default = 512,
								type = int)

	# First frame
	parser.add_argument("--first_frame",
								help = "First frame",
								type = int)

	# Last frame
	parser.add_argument("--last_frame",
								help = "Last frame",
								type = int)

	# Analyse control cases
	parser.add_argument("--control",
								help = "(default = False) If set to true only analyses tracks that appear in the first few frames. The number of such frames can be set with --control_frame_limit.",
								choices = ['True', 'False'],
								default = 'False')

	# Control frame limit
	parser.add_argument("--control_frame_limit",
								help = "(default = 3) Set the number of initial frames to use for analysis of control data.",
								default = 3,
								type = int)

	# Minimum track length for GTPases
	parser.add_argument("--gtpase_track_min_length",
								help = "(default = 5) Minimum track length of GTPases.",
								default = 5,
								type = int)

	# Output file name
	parser.add_argument("--outfile",
								help = "(default = Colocalization.csv) Output file name",
								default = "Colocalization.csv")

	args = parser.parse_args()

	return args

#--- Get data
def dataIN(filename):
	data = pd.read_csv(filename)
	return data

#--- Tracks for control cases
def get_control_tracks(data):
	# control tracks
	control_tracks = set(data[data["FRAME"] <= args.control_frame_limit]["PSEUDO_TRACK_ID"])

	# Keep only control tracks
	data = data[data["PSEUDO_TRACK_ID"].isin(control_tracks)]

	return data

#--- Eliminate tracks that appear before first frame
def eliminate_preexisting_tracks(data):
	# pre-existing tracks
	preexisting_tracks = set(data[data["FRAME"] < args.first_frame]["TRACK_ID"])
	preexisting_tracks = [n for n in preexisting_tracks if n != "None"]
	
	# eliminate tracks
	for track_id in preexisting_tracks:
		data = data[data["TRACK_ID"] != track_id]
	
	return data

#--- Remove spot data
def remove_spots(data):
	# removing spots
	data = data[data["TRACK_ID"] != "None"]
	return data

#--- Remove short tracks
def remove_short_tracks(data):
	# group by TRACK_ID and remove smaller tracks
	groupings = data.groupby(["TRACK_ID"])

	# long tracks
	long_track_ids = []

	# get TRACK_ID of long tracks
	for gid, group in groupings:
		min_frame = min(group["FRAME"])				# first frame
		max_frame = max(group["FRAME"])				# last frame
		track_length = max_frame - min_frame + 1	# frame length

		# keep track of long tracks
		if track_length >= args.gtpase_track_min_length:
			long_track_ids.append(gid)

	# keep only long tracks
	data = data[data["TRACK_ID"].isin(long_track_ids)]
	
	return data


#--- Keep specified frames
def filter_frames(data):

	# limit analysis to initial frames for control cases
	if args.control == "True":
		data_filtered = get_control_tracks(data)

	else:
		# filter frames according to user specifed first and last frame
		if (args.first_frame >= 0) and (args.last_frame >= 0):
			# eliminate tracks that originate before first frame
			data = eliminate_preexisting_tracks(data)
			# filter by start and end frame
			data_filtered = data[(data["FRAME"] >= args.first_frame) & (data["FRAME"] < args.last_frame)]

		# keep all frames after first
		elif args.first_frame >= 0:
			# eliminate tracks that originate before first frame
			data = eliminate_preexisting_tracks(data)
			# filter by start frame
			data_filtered = data[data["FRAME"] >= args.first_frame]

		# keep all frames up till the last
		else:
			data_filtered = data[data["FRAME"] < args.last_frame]

	return data_filtered

#--- Exclude spots outside field of view
def filter_fov(data):

	# field of view threshold
	threshold = args.image_size * args.pixel_size * args.field

	# outlier pseudo ids
	data_outliers = data.loc[(data["POSITION_X"] >= threshold) | (data["POSITION_Y"] >= threshold), "PSEUDO_TRACK_ID"]
	# outlier filter
	data_outliers_filter = ~data["PSEUDO_TRACK_ID"].isin(data_outliers)
	# filtering outliers
	data_fov = data[data_outliers_filter]

	return data_fov

#--- Calculate distance
def calc_dist(x1, y1, x2, y2):
	d = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)
	return round(d, 2)

#--- Get colocalizations for a single frame
def get_coloc_single(gtpase_data_frame, gdi_data_frame, frame, dist):

	# progress status
	print("\r# Frame: {}".format(frame), end="", flush=True)

	# colocalized GTPase spots
	gtpase_coloc = gtpase_data_frame
	gtpase_coloc["COLOCALIZED_SPOT"] = False
	gtpase_coloc["COLOCALIZATION_ID"] = np.nan
	
	# colocalized GDI spots
	gdi_coloc = gdi_data_frame
	gdi_coloc["COLOCALIZED_SPOT"] = False
	gdi_coloc["COLOCALIZATION_ID"] = np.nan

	# Calculate colocalizations in a single frame
	for gtpase_index in range(len(gtpase_data_frame)):
		for gdi_index in range(len(gdi_data_frame)):
			# GTPase spot coordinates
			x1 = gtpase_data_frame.iloc[gtpase_index,]["POSITION_X"]
			y1 = gtpase_data_frame.iloc[gtpase_index,]["POSITION_Y"]
			id1 = gtpase_data_frame.iloc[gtpase_index,]["PSEUDO_TRACK_ID"]

			# GDI spot coordinates
			x2 = gdi_data_frame.iloc[gdi_index,]["POSITION_X"]
			y2 = gdi_data_frame.iloc[gdi_index,]["POSITION_Y"]
			id2 = gdi_data_frame.iloc[gdi_index,]["PSEUDO_TRACK_ID"]

			# spot distance
			d = calc_dist(x1, y1, x2, y2)

			# store colocalized spots
			if d <= dist:
				# create colocalization id
				coloc_id = '{}-{}'.format(id1, id2)

				# Colocalized GTPase spot
				gtpase_coloc.loc[(gtpase_coloc["POSITION_X"] == x1) & (gtpase_coloc["POSITION_Y"] == y1), "COLOCALIZED_SPOT"] = True
				# Assign colocalization id to GTPase spot
				gtpase_coloc.loc[(gtpase_coloc["POSITION_X"] == x1) & (gtpase_coloc["POSITION_Y"] == y1), "COLOCALIZATION_ID"] = coloc_id

				# Colocalized GDI spot
				gdi_coloc.loc[(gdi_coloc["POSITION_X"] == x2) & (gdi_coloc["POSITION_Y"] == y2), "COLOCALIZED_SPOT"] = True
				# Assign colocalization id to GDI spot
				gdi_coloc.loc[(gdi_coloc["POSITION_X"] == x2) & (gdi_coloc["POSITION_Y"] == y2), "COLOCALIZATION_ID"] = coloc_id

	return [gtpase_coloc, gdi_coloc]

#--- Add pseudo track IDs
def add_PsedoTrackID(data):
	data = data.copy()
	data["PSEUDO_TRACK_ID"] = data.apply(lambda x: x["TRACK_ID"] if x["TRACK_ID"] != "None" else x["Label"], axis=1)
	return data

#--- Get all colocalizations
def get_coloc(gtpase_data, gdi_data):

	# total number of frames
	total_frames = min(max(gtpase_data["FRAME"]), max(gdi_data["FRAME"]))

	# processes for parallelization
	processes = []

	# Calculate colocalization for every pair of spots per frame
	for frame in range(0, total_frames):

		# fetch spots from single frame
		gtpase_data_frame = gtpase_data[gtpase_data["FRAME"] == frame]
		gdi_data_frame = gdi_data[gdi_data["FRAME"] == frame]
		
		# Skip frames that do not have any spots
		if (len(gtpase_data_frame) == 0) or (len(gdi_data_frame) == 0):
			continue

		# add sub-process
		processes.append((gtpase_data_frame, gdi_data_frame, frame, args.dist))

	# start multiprocessing
	pool = mp.Pool(mp.cpu_count())
	results = pool.starmap(get_coloc_single, processes)
	pool.close

	# combine colocalization results from single frames
	gtpase_coloc = pd.concat([gp[0] for gp in results if not gp[0].empty])
	gdi_coloc = pd.concat([gd[1] for gd in results if not gd[1].empty])

	# progress status
	print("")

	return (gtpase_coloc, gdi_coloc)

#--- Combine channels to single output
def combine_channels(gtpase_coloc, gdi_coloc):
	gtpase_coloc["CHANNEL"] = "GTPase"
	gdi_coloc["CHANNEL"] = "GDI"
	combine_channels = pd.concat([gtpase_coloc, gdi_coloc])

	return combine_channels

#--- Write output files
def dataOUT(data_frame, outname):
	# write parameters as meta data
	with open(outname, "w") as fh:
		fh.write("# {:=^40}\n".format(" Meta-data lines "))
		for arg in vars(args):
			fh.write("# {}: {}\n".format(arg, getattr(args, arg)))
		fh.write("# {:=^40}\n".format(" Colocalization data lines "))

	# write CSV file for colocalization events
	data_frame.to_csv(outname, index=False, float_format="%.3f", mode="a")

#--- Main function
def main():
	global args

	args = get_args()					# input arguments
	gtpase_data = dataIN(args.gtpase)	# GTPase data
	gdi_data = dataIN(args.gdi)			# GDI data

	# progress status
	print("# {:>20s} : {:^50s}".format("GTPase file", args.gtpase))
	print("# {:>20s} : {:^50s}".format("GDI file", args.gdi))
	print("# Data imported")

	# Add pseudo track IDs
	gtpase_data = add_PsedoTrackID(gtpase_data)
	gdi_data = add_PsedoTrackID(gdi_data)

	# progress status
	print("# Psedo Track IDs assigned")

	# remove spot data for GTPase channel
	if args.gtpase_track_min_length > 1:
		gtpase_data = remove_spots(gtpase_data)
		gtpase_data = remove_short_tracks(gtpase_data)

	# filter by first and last frame
	gtpase_data	= filter_frames(gtpase_data)
	gdi_data = filter_frames(gdi_data)

	# progress status
	print("# Frames filtered")

	# filter if custom field of view is set
	if args.field != 1.0:
		gtpase_data_fov = filter_fov(gtpase_data)	# filter field of view - GTPase
		gdi_data_fov = filter_fov(gdi_data)			# filter field of view - GDI
	else:
		gtpase_data_fov = gtpase_data
		gdi_data_fov = gdi_data

	# progress status
	print("# Field of view filtered")

	# Get colocalization
	gtpase_coloc, gdi_coloc = get_coloc(gtpase_data_fov, gdi_data_fov)

	# progress status
	print("# Colocalizations computed")

	# combine output
	combined_coloc = combine_channels(gtpase_coloc, gdi_coloc)

	# Write colocalization file
	dataOUT(combined_coloc, args.outfile)

	# progress status
	print("# Output written to: {:^50s}".format(args.outfile))


#--- Run main function
if __name__ == '__main__':
	main()

# Ankit Roy
# 11th November, 2020
# 11th April, 2021		
#	--> Added support for colocalization detection in multiple frames.
#	--> Colocalized spots are written to output files.
#	--> Improved code for filtering field of view.
# 15th April, 2021
#	--> Added support for specifying first/last frame.
#	--> Now reports progress.
# 16th April, 2021
#	--> Added support for chosing output files.
#	--> Parallelized colocalization calculations.
# 27th April, 2021
#	--> Removed multiple outputs to streamline script use.
#	--> Script now produces only 1 output file.
#	--> Added ability to specify output file name.
# 29th April, 2021
#	--> Added feature to assign unique colocalization id.
#	--> Colocalization id is written in output file.
# 9th August, 2021
#	--> Now adds input arguments as meta data to output files.
# 17th August, 2021
#	--> Now eliminates tracks that originate before first frame
#	--> Added option to keep tracks that originate before first frame
# 24th November, 2021
#	--> Explicitly stated conditions for first and last frame to be non-zero during frame filtering.
#	--> Fixed bug that previously eliminated frame 0 from results.
#	--> Added functionality to handle control cases separately.
#	--> For control cases only tracks originating in the first few frames are considered.
#	--> Users have option to specify the number of initial frames to consider.
# 25th November, 2021
#	--> Now uses only track data for the GTPase channel in control datasets.
# 26th November, 2021
#	--> Now uses only long tracks for the GTPase channel in control datasets.
#	--> Track length to define long tracks can be specified by the user.
# 28th January, 2022
#	--> Universal GTPase track length parameter added
#	--> Now eliminates entire tracks that move out of the field of view even momentarily
