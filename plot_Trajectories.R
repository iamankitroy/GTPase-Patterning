suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(argparse))

# This script is used to plot Single Molecule trajectories
# Trajectories of all tracks within a certain frames are plotted

# Create argument parser
parser = ArgumentParser()

parser$add_argument('-i', '--in_file',
					help = "Input file name")
parser$add_argument('-o', '--out_file',
					help = "Output file name")
parser$add_argument('-f', '--first_frame',
					default = 0,
					help = "First frame (default = 0)")
parser$add_argument('-l', '--last_frame',
					default = 200,
					help = "Last frame (default = 200)")
parser$add_argument('-c', '--channel',
					default = "GTPase",
					help = "Fluorescence channel (default = GTPase)")

args = parser$parse_args()

dt = 0.022			# time resolution in s
im_size = 512		# image size in px
px_size = 0.178		# pixel size in µm
fov = 0.75			# field of view

# Image field of view threshold
threshold = im_size * px_size * fov

# Crop region in pixels
crop_size = 168		# box length in px
crop_x = 97			# box x-min in px
crop_y = 106		# box y-min in px

# Convert crop values to µm
crop_size = crop_size * px_size
crop_x = crop_x * px_size
crop_y = crop_y * px_size

# Track data from colocalisation files
track_data = read.csv(args$in_file, comment.char = "#")

# Get track data subset for a certin time span in crop region
track_data_subset = track_data %>%
	subset(CHANNEL == args$channel) %>%									# specific channel
	subset(FRAME >= args$first_frame & FRAME <= args$last_frame) %>%	# specific time frame
	subset(POSITION_X >= crop_x & POSITION_X <= crop_x+crop_size) %>%	# crop region in x dimension
	subset(POSITION_Y >= crop_y & POSITION_Y <= crop_y+crop_size) %>%	# crop region in y dimension
	group_by(TRACK_ID) %>%												# group by TRACK_ID
	summarise(TRACK_ID,
			  FRAME,
			  POSITION_X,
			  POSITION_Y,
			  LIFETIME = FRAME - min(FRAME)) %>%						# calculate track lifetime
	ungroup()

# Convert time from frame to s
track_data_subset$FRAME_S = track_data_subset$FRAME * dt
track_data_subset$LIFETIME_S = track_data_subset$LIFETIME * dt

# Find track endpoints
end_points = track_data_subset %>%
	group_by(TRACK_ID) %>%
	summarise(TRACK_ID,
			  FRAME,
			  POSITION_X,
			  POSITION_Y,
			  MAX_FRAME = max(FRAME)) %>%
	ungroup() %>%
	subset(FRAME == MAX_FRAME)

# Trajectory plot
traj_plot = ggplot(track_data_subset,
				   aes(POSITION_X, POSITION_Y, group=factor(TRACK_ID), color=LIFETIME_S)) +
	geom_path(size = 1.2,
			  lineend = "round") +
	geom_point(data = end_points,
			   aes(POSITION_X, POSITION_Y, group=factor(TRACK_ID)), color='black', size=3) +
	scale_x_continuous(expand = c(0,0),
					   limits = c(crop_x, crop_x + crop_size)) +
	scale_y_continuous(expand = c(0,0),
					   limits = c(crop_y + crop_size, crop_y),
					   trans = "reverse") +
	scale_color_gradient(low = "#edf8b1",
						 high = "#253494",
						 limits = c(0, 4),
						 breaks = seq(0, 4, 1),
						 name = "Lifetime (s)") +
	coord_fixed() +
	theme_bw() +
	theme(panel.grid = element_blank(),
		  axis.title = element_blank(),
		  axis.text = element_blank(),
		  axis.ticks = element_blank(),
		  panel.background = element_blank(),
		  plot.background = element_rect(fill = "transparent",
		  							   colour = NA_character_),
		  legend.background = element_rect(fill = "transparent"),
		  legend.box.background = element_rect(fill = "transparent"),
		  legend.key = element_rect(fill = "transparent"),
		  legend.position = "none")

ggsave(filename = args$out_file, traj_plot)

# Ankit Roy
# 21st February, 2024