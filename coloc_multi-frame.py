#!/Users/roy/anaconda3/bin/python3

import argparse
import pandas as pd
import numpy as np

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

	args = parser.parse_args()

	return args

#--- Get data
def dataIN(filename):
	data = pd.read_csv(filename)
	return data

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

#--- Get colocalization
def get_coloc(gtpase_data, gdi_data):

	# colocalized GTPase spots
	gtpase_coloc = pd.DataFrame(columns = gtpase_data.columns.values)
	# colocalized GDI spots
	gdi_coloc = pd.DataFrame(columns = gdi_data.columns.values)

	# total number of frames
	total_frames = min(max(gtpase_data["FRAME"]), max(gdi_data["FRAME"]))

	# Calculate colocalization for every pair of spots per frame
	for frame in range(1, total_frames):
		gtpase_data_frame = gtpase_data[gtpase_data["FRAME"] == frame]
		gdi_data_frame = gdi_data[gdi_data["FRAME"] == frame]
		
		# Skip frames that do not have any spots
		if (len(gtpase_data_frame) == 0) or (len(gdi_data_frame) == 0):
			continue
		
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
				if d <= args.dist:
					# Store colocalized GTPase spot)
					gtpase_coloc = gtpase_coloc.append(gtpase_data_frame.iloc[gtpase_index,], ignore_index=True)
					# Store colocalized GDI spot
					gdi_coloc = gdi_coloc.append(gdi_data_frame.iloc[gdi_index,], ignore_index=True)

	return (gtpase_coloc, gdi_coloc)

#--- Write output files
def dataOUT(data_frame, filename):
	outname = "{}_coloc.csv".format(filename[:-4])

	# write CSV file for colocalization events
	data_frame.to_csv(outname, index=False)

#--- Main function
def main():
	global args

	args = get_args()					# input arguments
	gtpase_data = dataIN(args.gtpase)	# GTPase data
	gdi_data = dataIN(args.gdi)			# GDI data

	# filter if custom field of view is set
	if args.field != 1.0:
		gtpase_data_fov = filter_fov(gtpase_data)	# filter field of view - GTPase
		gdi_data_fov = filter_fov(gdi_data)			# filter field of view - GDI
	else:
		gtpase_data_fov = gtpase_data
		gdi_data_fov = gdi_data

	# Get colocalization
	gtpase_coloc, gdi_coloc = get_coloc(gtpase_data_fov, gdi_data_fov)

	gtpase_count = len(gtpase_data_fov)			# count GTPase spots
	gdi_count = len(gdi_data_fov)				# count GDI spots
	coloc_count = len(gtpase_coloc)				# count colocalized spots

	# fraction of colocalized GTPase spots
	gtpase_coloc_frac = round((coloc_count/gtpase_count)*100, 2)
	# fraction of colocalized GDI spots
	gdi_coloc_frac = round((coloc_count/gdi_count)*100, 2)

	print(gtpase_count, gdi_count, coloc_count, gtpase_coloc_frac, gdi_coloc_frac)

	# Write colocalization files
	dataOUT(gtpase_coloc, args.gtpase)
	dataOUT(gdi_coloc, args.gdi)


#--- Run main function
main()

# Ankit Roy
# 11th November, 2020
# 11th April, 2021		
#	--> added support for colocalization detection in multiple frames.
#	--> Colocalized spots are written to output files.
#	--> Improved code for filtering field of view.