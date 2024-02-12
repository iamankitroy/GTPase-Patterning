library(ggplot2)
library(dplyr)
library(argparse)

# Create argument parser
parser = ArgumentParser()

parser$add_argument('-i', '--in_file',
					help = "Input file name")
parser$add_argument('-o', '--out_file',
					help = "Output file name")

args = parser$parse_args()

# Read input position specific colocalization probability file
posProb_data = read.csv(args$in_file)

# Keep first 100 frames
posProb_data_subset = posProb_data %>%
	subset(FRAME <= 100)

# Calculate mean (mu) and standard deviation (sigma) of data
posProb_mu = mean(posProb_data_subset$PROBS)
posProb_sigma = sd(posProb_data_subset$PROBS)

# Calculate Z score for each data point
posProb_data_subset$ZSCORE = (posProb_data_subset$PROBS - posProb_mu)/posProb_sigma

# Write output file
write.csv(posProb_data_subset,
		  file = args$out_file,
		  row.names = FALSE,
		  quote = FALSE)

# Ankit Roy
# 12th February, 2024