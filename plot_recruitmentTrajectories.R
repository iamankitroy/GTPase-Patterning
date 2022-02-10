library(ggplot2)
library(dplyr)

dt = 0.022
im_size = 512
px_size = 0.178
fov = 0.75

threshold = im_size * px_size * fov

# path to plot files
path = "/Users/roy/Sandbox/recruitmentTrajectories_2022.02.10/all_files/"

# all plottable data files
files = list.files(path = path,
				   pattern = "*_5_subset.csv")

for (filename in files) {

	# plot file name
	plotname = paste(substr(filename, 0, nchar(filename)-4), "_trajPlot.png", sep = "")

	data = read.csv(paste(path, filename, sep = ""), comment.char = "#")

	data = data %>%
		subset(CHANNEL == "GTPase")

	data$TIME = data$FRAME * dt

	data_recruitment_tracks = data %>%
		subset(ANNOTATION_TRACK == "Recruitment")

	data_recruitment_spots = data %>%
		subset(ANNOTATION_SPOT == "Recruitment")

	trajplot = ggplot() +
		geom_path(data = data_recruitment_tracks, aes(x = POSITION_X, y = POSITION_Y, group = TRACK_ID, color = TIME)) +
		scale_color_gradient(low = "red", high = "blue") +
		geom_point(data = data_recruitment_spots, aes(x = POSITION_X, y = POSITION_Y), size = 2) +
		scale_x_continuous(limits = c(0, threshold),
						   expand = c(0, 0)) +
		scale_y_continuous(limits = c(0, threshold),
						   expand = c(0, 0)) +
		coord_fixed() +
		ggtitle(filename) +
		xlab("x-coordinate") +
		ylab("y-coordinate") +
		theme_bw()
	
	ggsave(plot = trajplot,
		   path = path,
		   filename = plotname)
}

# Ankit Roy
# 10th February, 2022