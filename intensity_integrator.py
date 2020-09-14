#!/usr/bin/python3

#--- This script is used to integrate fluorescence intensity across a line
#--- Input data needs to have the first 3 rows with description
#--- The first column should have the distance
#--- Subsequent columns should have fluorescence intensity
#--- Input file is a CSV with values seperated by ;
#--- Script written for integrating Johanna's fluorescence intensity data acquired from a line scan in Fiji

import sys
from numpy import trapz

#--- Get data in
def dataIN(filename):
	with open(filename, 'r') as fh:
		data = fh.readlines()
	return data

#--- Parse data
def parser(data):

	description = []			# store first 3 lines
	parsed_data = {}			# store parsed data in the form of a list of fluorescence intensities for every column

	#--- Parse data
	for row_index in range(len(data)):
		# Store description line as comma seperated line
		if row_index < 3:
			description.append(data[row_index].replace(';', ','))
			continue

		row = data[row_index]			# row instance
		columns = row.split(';')		# all columns

		for col_index in range(len(columns)):
			# Initiate column index
			if col_index not in parsed_data.keys():
				parsed_data[col_index] = []

			# Add 0 for missing values
			try:
				parsed_data[col_index].append(float(columns[col_index]))
			except:
				parsed_data[col_index].append(0.0)

	return parsed_data, description

#--- Integrate over distance
def integrate(data):

	ncols = len(data)					# number of columns
	distance = data[0]					# distance column

	all_aucs = {}						# store area under the curves

	dx = distance[1] - distance[0]		# delta x

	#--- Integrate
	for col in range(1, ncols):
		auc = trapz(data[col], dx = dx)	# integrate using the Trapezoidal rule
		all_aucs[col] = round(auc, 2)

	return all_aucs

#--- Write data out
def dataOUT(data, description, auc, filename):

	# Output file name with the extension _integrated.csv
	outname = filename[:-4]
	outname = "{}_integrated.csv".format(outname)

	with open(outname, "w") as fh:
		ncols = len(data)				# number of columns
		nrows = len(data[0])			# number of rows

		fh.write(''.join(description))	# write description

		# Write each row
		for row in range(nrows):
			for col in range(ncols):
				fh.write("{},".format(data[col][row]))
			fh.write("\n")

		# Write AUCs for every column
		for col in range(ncols):
			if col == 0:
				fh.write("AUC,")
			else:
				fh.write("{},".format(auc[col]))

#--- Main function
def main():
	filename = sys.argv[1]
	data = dataIN(filename)

	parsed_data, description = parser(data)
	auc = integrate(parsed_data)

	dataOUT(parsed_data, description, auc, filename)

#--- Run main function
main()

# Ankit Roy
# 11th September, 2020
# 14th September, 2020 			>	writes description in the output file
