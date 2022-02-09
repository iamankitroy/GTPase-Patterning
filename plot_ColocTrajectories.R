library(ggplot2)
library(dplyr)

data = read.csv("/Users/roy/Sandbox/plotProbabilities_positionSpecific_2022.02.02/MM_100_pM_complex_GDP_1_colocalisation_5_subset.csv", comment.char = "#")

data = data %>%
	subset(CHANNEL == "GTPase")

data_subset = data %>%
	subset(CHANNEL == "GTPase" &
		   	(ANNOTATION_TRACK == "Recruitment" |
		   	 	ANNOTATION_TRACK == "Recruitment and Extraction")) %>%
	subset(COLOCALIZED_SPOT == "True")

image_size = 512
px_size = 0.178
fov = 0.75

field_limit = image_size * px_size * fov

complete_trajectories = ggplot(data, aes(POSITION_X, POSITION_Y, color=FRAME, group=PSEUDO_TRACK_ID)) +
	geom_point() +
	geom_path() +
	scale_color_gradient2(midpoint = 500) +
	theme_bw() +
	scale_x_continuous(limits = c(0, field_limit), expand = c(0,0)) +
	scale_y_continuous(limits = c(0, field_limit), expand = c(0,0)) +
	xlab("x-coordinate") +
	ylab("y-coordinate")

colocalization_trajectories = ggplot(data_subset, aes(POSITION_X, POSITION_Y, color=FRAME, group=PSEUDO_TRACK_ID)) +
	geom_point() +
	geom_path() +
	scale_color_gradient2(midpoint = 500) +
	theme_bw() +
	scale_x_continuous(limits = c(0, field_limit), expand = c(0,0)) +
	scale_y_continuous(limits = c(0, field_limit), expand = c(0,0)) +
	xlab("x-coordinate") +
	ylab("y-coordinate")

complete_trajectories
colocalization_trajectories

# Ankit Roy
# 3rd February, 2022