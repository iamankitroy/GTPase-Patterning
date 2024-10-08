library(ggplot2)
library(dplyr)

# path to plot files
path = "/Users/roy/Sandbox/plotProbabilities_positionSpecific_2022.02.02/selected_files_28.03.2022/combined_replicates/"

# all plottable data files
files = list.files(path = path,
				   pattern = "*_combined.csv")

for (filename in files) {
	
	# plot file name
	plotname = paste(substr(filename, 0, nchar(filename)-4), ".svg", sep = "")
	
	data = read.csv(paste(path, filename, sep = ""), comment.char = "#")
	
	data = data %>% group_by(FRAME) %>% summarise(FRAME, PROBS, NOBS, FILE, SD=sd(PROBS), MEAN=mean(PROBS))
	
	probPlot = ggplot(data, aes(FRAME, PROBS, color = NOBS)) +
		geom_errorbar(aes(ymin = MEAN-SD, ymax = MEAN+SD),
					  position = position_dodge(0.9),
					  color = "black") +
		geom_point(aes(x = FRAME, y = MEAN),
				   shape=18,
				   color="black") +
		geom_point() +
		scale_x_continuous(limits = c(0, 100),
						   breaks = seq(0, 100, 5)) +
		scale_y_continuous(limits = c(-0.01, 0.05)) +
		xlab("Track Lifetime (frame)") +
		ylab("Colocalisation Probability") +
		ggtitle(filename) +
		scale_color_gradient(low = "blue",
							 high = "red",
							 name = "Number of Spots") +
		theme_bw() +
		theme(legend.position = "bottom",
			  legend.key.width = unit(2.5, "cm"))
	
	ggsave(filename = plotname,
		   path = path,
		   plot = probPlot,
		   height = 8,
		   width = 12,
		   units = "in")
	
}

# Ankit Roy
# Pi day, 2022