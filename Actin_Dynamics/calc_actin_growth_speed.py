#!/Users/roy/anaconda3/bin/python

#--- This script is used to calculate actin polymerization growth speeds
#--- Growth speeds are calculated using slope angles derived from kymographs
#--- Growth speeds are reported in monomers/s

import argparse
import pandas as pd
import numpy as np

#--- Fetch arguments
def get_args():
	parser = argparse.ArgumentParser()

	# Positional argument for file name
	parser.add_argument("filename",
						help = "Input file name")

	# Optional argument to set pixel size; default = 0.27 µm
	parser.add_argument("-ps", "--pixel_size",
						help = "(default = 0.27 µm) Pixel size ",
						type = float,
						default = 0.27)

	# Optional argument to set monomer size; default = 2.65 nm
	parser.add_argument("-ms", "--monomer_size",
						help = "(default = 2.65 nm) Monomer size ",
						type = float,
						default = 2.65)

	# Optional argument to set time interval; default = 5 s
	parser.add_argument("-t", "--time_interval",
						help = "(default = 5 s) Time interval between images ",
						type = float,
						default = 5)

	args = parser.parse_args()

	return args

#--- Get data
def dataIN(filename):
	data = pd.read_csv(filename)
	return data

#--- Calculate grwoth speeds
def calc_speed(data):

	# calculate absolute value of angle in radians
	data["absAngle_Radians"] = np.abs(data["Angle"] * np.pi/180)
	# calculate slope in pixel/frame
	data["Slope"] = 1/np.tan(data["absAngle_Radians"])
	# calculate rate in µm/s
	data["Rate_micron_per_min"] = data["Slope"] * args.pixel_size / args.time_interval
	# calculate rate in monomer/s
	data["Rate_monomer_per_min"] = round(data["Rate_micron_per_min"] * 1000 / args.monomer_size, 2)

	return data

#--- Write data
def dataOUT(data):
	# output file name
	outname = args.filename[:-4]
	outname = "{}_monomerRate.csv".format(outname)

	data.to_csv(outname, float_format='%.3f')

#--- Main function
def main():
	global args

	args = get_args()				# input arguments
	data = dataIN(args.filename)	# data
	data = calc_speed(data)			# calculate growth speed
	dataOUT(data)					# write output data file


#--- Run main
main()

# Ankit Roy
# 30th October, 2020
