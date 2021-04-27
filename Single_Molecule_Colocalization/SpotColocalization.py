#!/Users/roy/anaconda3/bin/python3

import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp

__author__ = "Ankit Roy"
__copyright__ = "Copyright 2021, Ankit Roy"
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

	# Output file name
	parser.add_argument("--outfile",
								help = "Output file name",
								default = "Colocalization.csv")

	args = parser.parse_args()

	return args

#--- Get data
def dataIN(filename):
	data = pd.read_csv(filename)
	return data

#--- Keep specified frames
def filter_frames(data):

	# filter frames according to user specifed first and last frame
	# keep frames between first and last
	if args.first_frame and args.last_frame:
		data_filtered = data[(data["FRAME"] >= args.first_frame) & (data["FRAME"] < args.last_frame)]
	# keep all frames after first
	elif args.first_frame:
		data_filtered = data[data["FRAME"] >= args.first_frame]
	# keep all frames up till the last
	else:
		data_filtered = data[data["FRAME"] < args.last_frame]

	return data_filtered

#--- Exclude spots outside field of view
def filter_fov(data):

	# field of view threshold
	threshold = args.image_size * args.pixel_size * args.field

	# filter spots inside field of view
	data_fov = data[(data["POSITION_X"] < threshold) & (data["POSITION_Y"] < threshold)]

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
	gtpase_coloc["COLOCALIZATION"] = False
	
	# colocalized GDI spots
	gdi_coloc = gdi_data_frame
	gdi_coloc["COLOCALIZATION"] = False

	# Calculate colocalizations in a single frame
	for gtpase_index in range(len(gtpase_data_frame)):
		for gdi_index in range(len(gdi_data_frame)):
			# GTPase spot coordinates
			x1 = gtpase_data_frame.iloc[gtpase_index,]["POSITION_X"]
			y1 = gtpase_data_frame.iloc[gtpase_index,]["POSITION_Y"]

			# GDI spot coordinates
			x2 = gdi_data_frame.iloc[gdi_index,]["POSITION_X"]
			y2 = gdi_data_frame.iloc[gdi_index,]["POSITION_Y"]

			# spot distance
			d = calc_dist(x1, y1, x2, y2)

			# store colocalized spots
			if d <= dist:
				# Colocalized GTPase spot
				gtpase_coloc.loc[(gtpase_coloc["POSITION_X"] == x1) & (gtpase_coloc["POSITION_Y"] == y1), "COLOCALIZATION"] = True

				# Colocalized GDI spot
				gdi_coloc.loc[(gdi_coloc["POSITION_X"] == x2) & (gdi_coloc["POSITION_Y"] == y2), "COLOCALIZATION"] = True

	return [gtpase_coloc, gdi_coloc]

#--- Get all colocalizations
def get_coloc(gtpase_data, gdi_data):

	# total number of frames
	total_frames = min(max(gtpase_data["FRAME"]), max(gdi_data["FRAME"]))

	# processes for parallelization
	processes = []

	# Calculate colocalization for every pair of spots per frame
	for frame in range(1, total_frames):

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
	# write CSV file for colocalization events
	data_frame.to_csv(outname, index=False, float_format="%.3f")

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

	# gtpase_count = len(gtpase_data_fov)			# count GTPase spots
	# gdi_count = len(gdi_data_fov)				# count GDI spots
	# coloc_count = len(gtpase_coloc)				# count colocalized spots

	# # fraction of colocalized GTPase spots
	# gtpase_coloc_frac = round((coloc_count/gtpase_count)*100, 2)
	# # fraction of colocalized GDI spots
	# gdi_coloc_frac = round((coloc_count/gdi_count)*100, 2)

	# print(args.gtpase, args.gdi, gtpase_count, gdi_count, coloc_count, gtpase_coloc_frac, gdi_coloc_frac)


	# # Write colocalization files
	# dataOUT(gtpase_coloc, args.gtpase)
	# dataOUT(gdi_coloc, args.gdi)

	# combine output
	combined_coloc = combine_channels(gtpase_coloc, gdi_coloc)

	# Write colocalization file
	dataOUT(combined_coloc, args.outfile)

	# progress status
	print("# Output files written")


#--- Run main function
if __name__ == '__main__':
	main()

# Ankit Roy
# 11th November, 2020
# 11th April, 2021		
#	--> added support for colocalization detection in multiple frames.
#	--> colocalized spots are written to output files.
#	--> improved code for filtering field of view.
# 15th April, 2021
#	--> added support for specifying first/last frame.
#	--> now reports progress.
# 16th April, 2021
#	--> added support for chosing output files.
#	--> Parallelized colocalization calculations.
# 27th April, 2021
#	--> removed multiple outputs to streamline script use.
#	--> script now produces only 1 output file.
#	--> added ability to specify output file name.