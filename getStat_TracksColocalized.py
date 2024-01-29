#!/Users/roy/anaconda3/bin/python3

import argparse
import pandas as pd
import numpy as np
import sys
from collections import Counter

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

	# Spot colocalization file
	required_args.add_argument("-cf", "--colocalization_file",
								help = "Spot colocalization file name.",
								required = True)

	# Limit GDI free frames
	parser.add_argument("--limit_free_gdi",
								help = "(default = True) Limit the number of frames GDI spots can remain un-colocalized.",
								choices = ['True', 'False'],
								default = 'True')

	# Maximum free frames in GDI channel
	parser.add_argument("--max_gdi_free_frames",
								help = "(default = 3) Set the maximum number of frames GDI spots can remain un-colocalized.",
								type = int,
								default = 3)

	# Set the number of frames to consider for a recruitment event
	parser.add_argument("--recruitment_frames",
								help = "(default = 3) Set the number of frames from the start of a GTPase track in which a GDI colocalization is considered a true recruitment event.",
								type = int,
								default = 3)

	# Set the number of frames to consider for an extraction event
	parser.add_argument("--extraction_frames",
								help = "(default = 3) Set the number of frames before the end of a GTPase track in which a GDI colocalization is considered a true extraction event.",
								type = int,
								default = 3)

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

	# Time resolution
	parser.add_argument("-tr", "--time_resolution",
								help = "(default = 0.022 s) Time resolution",
								default = 0.022,
								type = float)
	
	args = parser.parse_args()

	return args

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

#--- Modify field
def modifyField(data, condition, name, value):
	data.loc[condition, name] = value
	return data

#--- Get colocalized tracks
def get_ColocalizedTracks(data):
	data = data.copy()
	data["COLOCALIZED_TRACK"] = False

	# group data with pseudo track IDs and channels
	grouping = data.groupby(["PSEUDO_TRACK_ID", "CHANNEL"])
	
	# identify colocalized tracks
	# any track with at least one colocalized spot is considered a colocalized track
	for gid, group in grouping:
		if group["COLOCALIZED_SPOT"].sum() > 0:
			# group condition
			condition = (data["PSEUDO_TRACK_ID"] == gid[0]) & (data["CHANNEL"] == gid[1])
			# modify group field "COLOCALIZED_TRACK" to True
			data = modifyField(data, condition, "COLOCALIZED_TRACK", True)

	return data

#--- Subset colocalized tracks
def subsetData(data):
	data = data.copy()
	data = data[data["COLOCALIZED_TRACK"]]
	return data

#--- Count the number of colocalized frames
def count_ColocalizedFrames(data):
	data = data.copy()
	data["TOTAL_FRAME_COUNT"] = 0
	data["COLOCALIZED_FRAME_COUNT"] = 0

	# group data with pseudo track IDs and channels
	grouping = data.groupby(["PSEUDO_TRACK_ID", "CHANNEL"])

	for gid, group in grouping:
		f_count = len(group)						# total frame count
		cf_count = group["COLOCALIZED_SPOT"].sum()  # colocalized frame count

		# group condition
		condition = (data["PSEUDO_TRACK_ID"] == gid[0]) & (data["CHANNEL"] == gid[1])

		# add total frame count
		data = modifyField(data, condition, "TOTAL_FRAME_COUNT", f_count)
		# add colocalized frame count
		data = modifyField(data, condition, "COLOCALIZED_FRAME_COUNT", cf_count)


	# add free frame count
	data["FREE_FRAME_COUNT"] = data.apply(lambda x: x["TOTAL_FRAME_COUNT"] - x["COLOCALIZED_FRAME_COUNT"], axis=1)
	# add colocalized frame fraction
	data["COLOCALIZED_FRAME_FRACTION"] = data.apply(lambda x: x["COLOCALIZED_FRAME_COUNT"]/x["TOTAL_FRAME_COUNT"], axis=1)

	return data

#--- Channel based filter
def filter_channel(data, channel, pseudo_ids):
	data = data.copy()
	channel_filter = data["CHANNEL"] == channel					# filter CHANNEL
	pseudo_id_filter = data["PSEUDO_TRACK_ID"].isin(pseudo_ids) # filter PSEUDO_TRACK_ID

	data = data.loc[channel_filter & pseudo_id_filter, ]		# combine filters

	return data

#--- Filter GDI tracks based on un-colocalized frames
def filter_freeFrames(data, count):
	data = data.copy()

	# COLOCALIZATION_ID of GDI tracks that remain un-colocalized for <= user defined number of frames
	coloc_ids = set(data.loc[(data["FREE_FRAME_COUNT"] <= count) & (data["CHANNEL"] == "GDI"), "COLOCALIZATION_ID"])
	coloc_ids = coloc_ids - {np.nan}							# remove NaN
	
	# GTPase and GDI ids of desired COLOCALIZATION_ID
	gtpase_ids = [i.split('-')[0] for i in coloc_ids]
	gdi_ids = [i.split('-')[1] for i in coloc_ids]

	# Apply filters and combine
	gtpase_filtered = filter_channel(data, "GTPase", gtpase_ids)
	gdi_filtered = filter_channel(data, "GDI", gdi_ids)
	data_filtered = pd.concat([gtpase_filtered, gdi_filtered])

	return data_filtered

#--- Identify true recruitment/extraction events
def get_trueEvents(data, coloc_id, condition, frame_threshold):
	# GTPase PSEUDO_TRACK_ID
	gtpase_id = coloc_id.split('-')[0]
	
	# Identify recruitment events
	if condition == "recruitment":
		# Error handeled for orphan spots
		try:
			# GTPase enters
			enter_gtpase = min(data.loc[(data["PSEUDO_TRACK_ID"] == gtpase_id) & (data["CHANNEL"] == "GTPase"), "FRAME"])
			# First colocalization event of GTPase
			first_coloc = min(data.loc[(data["COLOCALIZATION_ID"] == coloc_id) & (data["CHANNEL"] == "GTPase"), "FRAME"])
			frame_diff = first_coloc - enter_gtpase

			# True if colocalization within first few frames of GTPase
			return (first_coloc - enter_gtpase) < frame_threshold
		except:
			return False

	else:
		# Error handeled for orphan spots
		try:
			# GTPase exits
			exit_gtpase = max(data.loc[(data["PSEUDO_TRACK_ID"] == gtpase_id) & (data["CHANNEL"] == "GTPase"), "FRAME"])
			# Last colocalization event of GTPase
			last_coloc = max(data.loc[(data["COLOCALIZATION_ID"] == coloc_id) & (data["CHANNEL"] == "GTPase"), "FRAME"])

			# True if colocalization within last few frames of GTPase
			return (exit_gtpase - last_coloc) < frame_threshold
		except:
			return False

#--- Annotate spot
def annotateSpots(data, recruitment_events, extraction_events):
	data["ANNOTATION_SPOT"] = np.nan						# initiate spot annotation

	# Annotate extracted spots
	data["ANNOTATION_SPOT"] = data.apply(lambda x: "Extraction" if x["COLOCALIZATION_ID"] in extraction_events else x["ANNOTATION_SPOT"], axis = 1)
	# Annotate recruited spots
	data["ANNOTATION_SPOT"] = data.apply(lambda x: "Recruitment" if x["COLOCALIZATION_ID"] in recruitment_events else x["ANNOTATION_SPOT"], axis = 1)

	return data

#--- Get PSEUDO_TRACK_ID from COLOCALIZATION_ID
def getTrack_from_Colocalization(coloc_ids, molecule):
	# Split COLOCALIZATION_ID to get PSEUDO_TRACK_ID
	# molecule 0 --> GTPase
	# molecule 1 --> GDI
	pseudo_track_ids = [pair.split('-')[molecule] for pair in coloc_ids]

	return pseudo_track_ids

#--- Annotate tracks
def annotateTracks(data, recruitment_events, extraction_events):
	data["ANNOTATION_TRACK"] = "None"

	# Recruited and/or extracted tracks
	recruited_gtpase = set(getTrack_from_Colocalization(recruitment_events, 0))
	extracted_gtpase = set(getTrack_from_Colocalization(extraction_events, 0))
	recruited_gdi = set(getTrack_from_Colocalization(recruitment_events, 1))
	extracted_gdi = set(getTrack_from_Colocalization(extraction_events, 1))

	# with pd.option_context('display.max_rows', None, 'display.max_columns', None):
	#	 print(data.loc[(data["CHANNEL"] == "GTPase") & (data["PSEUDO_TRACK_ID"].isin(recruited_gtpase)), ].groupby(["PSEUDO_TRACK_ID"]).agg('count')["FRAME"])
	
	# Recruited and extracted tracks
	intersection_gtpase = recruited_gtpase.intersection(extracted_gtpase)
	intersection_gdi = recruited_gdi.intersection(extracted_gdi)

	# print(len(recruited_gtpase), len(recruited_gdi))
	# print(len(extracted_gtpase), len(extracted_gdi))
	# testthis = sorted([c.split('-')[1] for g in recruited_gtpase for c in recruitment_events if g == c.split('-')[0]])
	# print(Counter(testthis))

	# Exclusively recruited or extracted tracks
	recruited_gtpase = recruited_gtpase - intersection_gtpase
	extracted_gtpase = extracted_gtpase - intersection_gtpase
	recruited_gdi = recruited_gdi - intersection_gdi
	extracted_gdi = extracted_gdi - intersection_gdi

	# Annotate tracks
	# Recruited and extracted GTPases
	data["ANNOTATION_TRACK"] = data.apply(lambda x: "Recruitment and Extraction" if ((x["CHANNEL"] == "GTPase") and (x["PSEUDO_TRACK_ID"] in intersection_gtpase)) else x["ANNOTATION_TRACK"], axis = 1)
	# Only recruited GTPases
	data["ANNOTATION_TRACK"] = data.apply(lambda x: "Recruitment" if ((x["CHANNEL"] == "GTPase") and (x["PSEUDO_TRACK_ID"] in recruited_gtpase)) else x["ANNOTATION_TRACK"], axis = 1)
	# Only extracted GTPases
	data["ANNOTATION_TRACK"] = data.apply(lambda x: "Extraction" if ((x["CHANNEL"] == "GTPase") and (x["PSEUDO_TRACK_ID"] in extracted_gtpase)) else x["ANNOTATION_TRACK"], axis = 1)
	# Recruited and extracted GDI
	data["ANNOTATION_TRACK"] = data.apply(lambda x: "Recruitment and Extraction" if ((x["CHANNEL"] == "GDI") and (x["PSEUDO_TRACK_ID"] in intersection_gdi)) else x["ANNOTATION_TRACK"], axis = 1)
	# Only recruited GDI
	data["ANNOTATION_TRACK"] = data.apply(lambda x: "Recruitment" if ((x["CHANNEL"] == "GDI") and (x["PSEUDO_TRACK_ID"] in recruited_gdi)) else x["ANNOTATION_TRACK"], axis = 1)
	# Only extracted GDI
	data["ANNOTATION_TRACK"] = data.apply(lambda x: "Extraction" if ((x["CHANNEL"] == "GDI") and (x["PSEUDO_TRACK_ID"] in extracted_gdi)) else x["ANNOTATION_TRACK"], axis = 1)

	# testthis = sorted([c.split('-')[1] for g in recruited_gtpase for c in recruitment_events if g == c.split('-')[0]])
	# print(testthis)
	# print(len(testthis))
	# print(len(set(testthis)))

	# print("GTPase: Recruited - {}, Extracted - {}, Both - {}".format(len(recruited_gtpase), len(extracted_gtpase), len(intersection_gtpase)))
	# print("GDI: Recruited - {}, Extracted - {}, Both - {}".format(len(recruited_gdi), len(extracted_gdi), len(intersection_gdi)))

	return data

#--- Annotate events as recruitment or extraction events
def annotateEvents(data, recruitment_frame_threshold, extraction_frame_threshold):
	data = data.copy()

	# COLOCALIZATION_IDs of all GDI spots
	coloc_ids = set(data.loc[data["CHANNEL"] == "GDI", "COLOCALIZATION_ID"])
	coloc_ids = coloc_ids - {np.nan}

	# All true recruitment events
	recruitment_events = [c for c in coloc_ids if get_trueEvents(data, c, "recruitment", recruitment_frame_threshold)]

	# All true extraction events
	extraction_events = [c for c in coloc_ids if get_trueEvents(data, c, "extraction", extraction_frame_threshold)]

	# Annotate spots
	data = annotateSpots(data, recruitment_events, extraction_events)

	# Annotate tracks
	data = annotateTracks(data, recruitment_events, extraction_events)

	# with pd.option_context('display.max_rows', None, 'display.max_columns', None):
	#	 print(data.loc[(data["CHANNEL"] == "GTPase") & (data["ANNOTATION_TRACK"] == "Recruitment and Extraction"), ].groupby("PSEUDO_TRACK_ID").agg('count')[["FRAME", "ANNOTATION_TRACK"]])

	return data

#--- Calculate landing rate
def calcLandingRate(data):
	fov = args.field				# field of view
	image_size = args.image_size	# image size in pixels
	pixel_size = args.pixel_size	# pixel size in µm
	time_resolution = args.time_resolution	# time resolution in s

	# Total number of frames
	total_frames = max(set(data["FRAME"]))
	
	# Number of landing events
	gtpase_landing = len(set(data.loc[data["CHANNEL"] == "GTPase", "PSEUDO_TRACK_ID"]))
	gdi_landing = len(set(data.loc[data["CHANNEL"] == "GDI", "PSEUDO_TRACK_ID"]))

	# Landing rate
	gtpase_landing_rate = gtpase_landing/(total_frames * time_resolution * (image_size * pixel_size * fov)**2)
	gdi_landing_rate = gdi_landing/(total_frames * time_resolution * (image_size * pixel_size * fov)**2)
	
	return (gtpase_landing_rate, gdi_landing_rate)

#--- All colocalization statistics:
def getStat(all_data, subset_data, input_file):
	# Count all GTPase tracks
	total_gtpase_tracks = len(set(all_data.loc[all_data["CHANNEL"] == "GTPase", "PSEUDO_TRACK_ID"]))
	# Count all GDI tracks
	total_gdi_tracks = len(set(all_data.loc[all_data["CHANNEL"] == "GDI", "PSEUDO_TRACK_ID"]))
	# Count colocalized GTPase tracks
	coloc_gtpase_tracks = len(set(subset_data.loc[subset_data["CHANNEL"] == "GTPase", "PSEUDO_TRACK_ID"]))
	# Count colocalized GDI tracks
	coloc_gdi_tracks = len(set(subset_data.loc[subset_data["CHANNEL"] == "GDI", "PSEUDO_TRACK_ID"]))
	# Percentage of colocalized GTPase tracks
	percent_coloc_gtpase = coloc_gtpase_tracks/total_gtpase_tracks * 100
	# Percentage of colocalized GDI tracks
	percent_coloc_gdi = coloc_gdi_tracks/total_gdi_tracks * 100
	# Count recruited GTPases
	recruited_gtpase = len(set(subset_data.loc[(subset_data["CHANNEL"] == "GTPase") & (subset_data["ANNOTATION_TRACK"] == "Recruitment"), "PSEUDO_TRACK_ID"]))
	# Count extracted GTPases
	extracted_gtpase = len(set(subset_data.loc[(subset_data["CHANNEL"] == "GTPase") & (subset_data["ANNOTATION_TRACK"] == "Extraction"), "PSEUDO_TRACK_ID"]))
	# Count recruited and extracted GTPases
	intersection_gtpase = len(set(subset_data.loc[(subset_data["CHANNEL"] == "GTPase") & (subset_data["ANNOTATION_TRACK"] == "Recruitment and Extraction"), "PSEUDO_TRACK_ID"]))
	# Count recruited GDI
	recruited_gdi = len(set(subset_data.loc[(subset_data["CHANNEL"] == "GDI") & (subset_data["ANNOTATION_TRACK"] == "Recruitment"), "PSEUDO_TRACK_ID"]))
	# Count extracted GTPases
	extracted_gdi = len(set(subset_data.loc[(subset_data["CHANNEL"] == "GDI") & (subset_data["ANNOTATION_TRACK"] == "Extraction"), "PSEUDO_TRACK_ID"]))
	# Count recruited and extracted GTPases
	intersection_gdi = len(set(subset_data.loc[(subset_data["CHANNEL"] == "GDI") & (subset_data["ANNOTATION_TRACK"] == "Recruitment and Extraction"), "PSEUDO_TRACK_ID"]))
	# Percentage of recruited GTPases
	percentage_recruited_gtpase = recruited_gtpase/total_gtpase_tracks * 100
	# Percentage of extracted GTPases
	percentage_extracted_gtpase = extracted_gtpase/total_gtpase_tracks * 100
	# Percentage of recruited and extracted GTPases
	percentage_intersection_gtpase = intersection_gtpase/total_gtpase_tracks * 100
	# Percentage of recruited GDI
	percentage_recruited_gdi = recruited_gdi/total_gdi_tracks * 100
	# Percentage of extracted GDI
	percentage_extracted_gdi = extracted_gdi/total_gdi_tracks * 100
	# Percentage of recruited and extracted GDI
	percentage_intersection_gdi = intersection_gdi/total_gdi_tracks * 100

	# Calculate landing rate
	gtpase_landing_rate, gdi_landing_rate = calcLandingRate(all_data)

	# Display stats
	header = "# {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format("File_Name",
													"GTPase_Count",
													"GDI_Count",
													"Colocalized_GTPase_Count",
													"Colocalized_GDI_Count",
													"Percent_Colocalized_GTPase",
													"Percent_Colocalized_GDI",
													"GTPase_Landing_Rate",
													"GDI_Landing_Rate",
													"Recruited_GTPase_Count",
													"Extracted_GTPase_Count",
													"Recruited_and_Extracted_GTPase_Count",
													"Recruited_GDI_Count",
													"Extracted_GDI_Count",
													"Recruited_and_Extracted_GDI_Count",
													"Recruited_GTPase_Percent",
													"Extracted_GTPase_Percent",
													"Recruited_and_Extracted_GTPase_Percent",
													"Recruited_GDI_Percent",
													"Extracted_GDI_Percent",
													"Recruited_and_Extracted_GDI_Percent")
	stat_line = "# {},{},{},{},{},{:.2f},{:.2f},{:.2E},{:.2E},{},{},{},{},{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}\n".format(input_file,
													total_gtpase_tracks,
													total_gdi_tracks,
													coloc_gtpase_tracks,
													coloc_gdi_tracks,
													percent_coloc_gtpase,
													percent_coloc_gdi,
													gtpase_landing_rate,
													gdi_landing_rate,
													recruited_gtpase,
													extracted_gtpase,
													intersection_gtpase,
													recruited_gdi,
													extracted_gdi,
													intersection_gdi,
													percentage_recruited_gtpase,
													percentage_extracted_gtpase,
													percentage_intersection_gtpase,
													percentage_recruited_gdi,
													percentage_extracted_gdi,
													percentage_intersection_gdi)

	return (header, stat_line)

#--- Write output files
def dataOUT(data_frame, outname, header, stat_line):
	outname = "{}_subset.csv".format(outname.split('.csv')[0])
	
	# write header and statistics
	with open(outname, 'w') as fh:
		fh.write("# {:=^40}\n".format(" Meta-data lines "))
		for arg in vars(args):
			fh.write("# {}: {}\n".format(arg, getattr(args, arg)))
		fh.write("# {:=^40}\n".format(" Summary lines "))
		fh.write(header)
		fh.write(stat_line)
		fh.write("# {:=^40}\n".format(" Colocalization data lines "))

	# write CSV file for colocalization events
	data_frame.to_csv(outname, index=False, float_format="%.3f", mode = 'a')

#--- Main function
def main():
	global args
	pd.set_option('display.max_columns', None)

	args = get_args()								# input arguments
	coloc_data = dataIN(args.colocalization_file)	# Colocalization data

	# Get colocalized tracks
	coloc_data = get_ColocalizedTracks(coloc_data)

	# Subset of colocalized tracks
	sub_coloc_data = subsetData(coloc_data)

	# Exit if no colocalizations are found
	if sub_coloc_data.empty:
		print("No colocalization found!")
		sub_coloc_data["TOTAL_FRAME_COUNT"] = 0
		sub_coloc_data["COLOCALIZED_FRAME_COUNT"] = 0
		sub_coloc_data["FREE_FRAME_COUNT"] = 0
		sub_coloc_data["COLOCALIZED_FRAME_FRACTION"] = 0
		sub_coloc_data["ANNOTATION_SPOT"] = 0
		sub_coloc_data["ANNOTATION_TRACK"] = 0

	# Get number of colocalized frames
	sub_coloc_data = count_ColocalizedFrames(sub_coloc_data)

	# Filter GDI spots which stay un-colocalized for more than a threshold number of frames	
	if args.limit_free_gdi == 'True':
		# Filter GDI tracks based on the number of uncolocalized frames
		sub_coloc_data = filter_freeFrames(sub_coloc_data, args.max_gdi_free_frames)

	# Annotate recruitment events
	sub_coloc_data = annotateEvents(sub_coloc_data, args.recruitment_frames, args.extraction_frames)

	# Show colocalization statistics
	header, stat_line = getStat(coloc_data, sub_coloc_data, args.colocalization_file)

	# Wtite output file
	dataOUT(sub_coloc_data, args.colocalization_file, header, stat_line)

#--- Run main function
if __name__ == '__main__':
	main()

# Ankit Roy
# 29th April, 2021
# 3rd May, 2021
#   --> Added ability to filter based on the number of free frames in GDI tracks
#   --> Added ability to filter based on the fraction of colocalized frames in GDI tracks
#   --> Displays track based colocalization statistics
# 6th April, 2021
#   --> Now calculates landing rate in terms of tracks per µm^2
#   --> Stores colocalization statistics along with data in output CSV file
# 9th August, 2021
#   --> Colocalization file argument set to a required argument
#   --> Ignores comment lines from colocalization file
#   --> Now adds input arguments as meta-data to output files
# 18th August, 2021
#   --> Fixed function to filter the number of frames a GDI spot can remain un-colocalized for
# 19th August, 2021
#   --> Added functionality to annotate colocalization events as recruitment/extraction events
#   --> Number of frames from the start/end of GTPase track to consider for recruitment/extraction can now be supplied via arguments
# 21st August, 2021
#   --> Now writes the number and percentages of recruitment and/or extraction events in the output file
# 14th January, 2022
#	--> Fixed an issue with orphan tracks
# 17th January, 2022
#	--> Execution stops if no colocalization is found
# 20th January, 2022
#	--> Complete script is executed and a subset file is generated even if no colocalizations are found
# 9th February, 2022
#	--> Input data types for some columns are explicitly mentioned to avoid pandas from guessing data types
#	--> Annotates recruitment spots after annotating extraction spots to prevent extraction events from overwriting recruitments
# 29th February, 2024
#	--> BUG: Forgot to square the image dimensions to calculated landing rate in /frame/µm^2. This has now been fixed.
#	--> Updated calcLandingRate function to calculate landing rate in units of /s/µm^2 instead of /frame/µm^2.