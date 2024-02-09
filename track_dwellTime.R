library(dplyr)
library(argparse)

# This script is used to calculate the 1-ECDF distrubution of single molecule track dwell time.

# Create argument parser
parser = ArgumentParser()

parser$add_argument('-i', '--in_file',
					help = "Input file name")
parser$add_argument('-o', '--out_file',
					help = "Output file name")

args = parser$parse_args()

# TrackMate track data
track_data = read.csv(args$in_file)

# Time resolution of a frame
time_resolution = 0.022

# Remove spots and keep only tracks
track_data_subset = track_data %>%
	subset(TRACK_ID != "None")

# Calculate lifetime
lifetime_data = track_data_subset %>%
	group_by(TRACK_ID) %>%
	summarise(TRACK_ID, FRAME_COUNT = n()) %>%
	ungroup() %>%
	unique() %>%
	mutate(DWELL_TIME = FRAME_COUNT * time_resolution)

# Calculate ECDF function from lifetime data
ecdf_fn = ecdf(lifetime_data$DWELL_TIME)

# Create dwell time data frame
dwell_data = data.frame(
	DWELL_TIME = sort(unique(lifetime_data$DWELL_TIME))
	)

# Store 1-ECDF for dwell time data
dwell_data$ECDF_INV = 1 - ecdf_fn(dwell_data$DWELL_TIME)

# Write output file
write.csv(dwell_data,
		  file = args$out_file,
		  row.names = FALSE,
		  quote = FALSE)

# Ankit Roy
# 9th February, 2024