#!/Users/roy/anaconda3/bin/python3

import argparse
import pandas as pd
import numpy as np

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

	args = parser.parse_args()

	return args

#--- Get data
def dataIN(filename):
	data = pd.read_csv(filename)
	return data

#--- Exclude spots outside field of view
def filter_fov(data):

	# spot data in the field of view
	data_fov = pd.DataFrame(columns = data.columns.values)

	# field of view threshold
	threshold = args.image_size * args.pixel_size * args.field

	# filter spots inside field of view
	for index in range(len(data)):
		x = data.loc[index,]["POSITION_X"]
		y = data.loc[index,]["POSITION_Y"]

		# ignore spots outside field of view
		if (x > threshold) or (y > threshold):
			continue

		# add spots within field of view
		data_fov = data_fov.append(data.loc[index,], ignore_index=True)

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

	# Calculate colocalization for every pair of spots
	for gtpase_index in range(len(gtpase_data)):
		for gdi_index in range(len(gdi_data)):
			# GTPase spot coordinates
			x1 = gtpase_data.loc[gtpase_index,]["POSITION_X"]
			y1 = gtpase_data.loc[gtpase_index,]["POSITION_Y"]

			# GDI spot coordinates
			x2 = gdi_data.loc[gdi_index,]["POSITION_X"]
			y2 = gdi_data.loc[gdi_index,]["POSITION_Y"]

			# spot distance
			d = calc_dist(x1, y1, x2, y2)

			# store colocalized spots
			if d <= args.dist:
				# Store colocalized GTPase spot
				gtpase_coloc = gtpase_coloc.append(gtpase_data.loc[gtpase_index,], ignore_index=True)
				# Store colocalized GDI spot
				gdi_coloc = gdi_coloc.append(gdi_data.loc[gdi_index,], ignore_index=True)

	return (gtpase_coloc, gdi_coloc)

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


#--- Run main function
main()

# Ankit Roy
# 11th November, 2020
