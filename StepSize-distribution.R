suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(argparse))

# Create argument parser
parser = ArgumentParser()

parser$add_argument('-i', '--in_file',
					help = "Input file name")

args = parser$parse_args()

# Output file prefix
args$out_file = substr(args$in_file, 1, nchar(args$in_file) - 4)
# Step size data file
args$msd_file = paste(args$out_file, "_stepsize-data.csv", sep = "")
# Step size distribution plot file
args$plot_file = paste(args$out_file, "_stepsize-dist_plot.png", sep = "")

# Read TrackMate file
track_data = read.csv(args$in_file, comment.char = "#")

# Remove spots that are not part of tracks
track_data = track_data %>%
	subset(TRACK_ID != "None")

# Calculate step size 
stepsize_data = track_data %>%
	group_by(TRACK_ID) %>%
	arrange(FRAME) %>%
	summarise(TRACK_ID,
			  FRAME,
			  POSITION_X,
			  POSITION_Y,
			  DEL_X = POSITION_X - lag(POSITION_X, default = first(POSITION_X)),
			  DEL_Y = POSITION_Y - lag(POSITION_Y, default = first(POSITION_Y))) %>%
	ungroup() %>%
	mutate(SQ_DISP = DEL_X^2 + DEL_Y^2) %>%
	mutate(STEP_SIZE = (SQ_DISP)^0.5)

# Write step size data file
write.csv(stepsize_data,
		  file = args$msd_file,
		  row.names = FALSE)

# Plot step size distribution
stepsize_plot = ggplot(stepsize_data, aes(STEP_SIZE)) +
	geom_density(color="darkgrey", fill="grey", alpha=0.5) +
	xlab(bquote("Step size ["*µm*"]")) +
	ylab("Probability density") +
	ggtitle(args$out_file) +
	scale_x_continuous(limits = c(0, 0.81), expand = c(0, 0)) +
	scale_y_continuous(limits = c(0, 7.5), expand = c(0, 0), breaks = seq(0, 7, 1)) +
	coord_equal(0.05) +
	theme_classic() +
	theme(panel.grid = element_blank())

# Save plot
ggsave(args$plot_file,
	   plot = stepsize_plot,
	   height = 3,
	   width = 6,
	   dpi = 300)

# Ankit Roy
# 20th March, 2024