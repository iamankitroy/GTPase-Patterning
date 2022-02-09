library(ggplot2)

path = "/Users/roy/Sandbox/plotProbabilities_2022.01.31/"			# file path

# all plottable data files
files = list.files(path = path,
				   pattern = "*_probPlot.csv")

for (filename in files) {

	# plot file name
	plotname = paste(substr(filename, 0, nchar(filename)-3), ".png", sep = "")
	#plotname = paste(substr(filename, 0, nchar(filename)-3), "_full.png", sep = "")
	
	# probability data
	probData = read.csv(paste(path, filename, sep = ""))
	probData$Classes = factor(probData$Classes)

	# plot
	probPlot_complete = ggplot(probData, aes(Classes, Probabilities)) +
		geom_boxplot() +
		stat_summary(fun.y=mean, geom="point", shape=20, size=8, color="red", fill="red") +
		geom_jitter(width = 0.2) +
		xlab("Colocalization Type") +
		ylab("Probability (colocalization/frame)") +
		ggtitle(filename)

	probPlot_truncated = probPlot_complete +
		scale_y_continuous(limits = c(0, 0.2))

	#ggsave(plot = probPlot_complete,
	#	   path = path,
	#	   filename = plotname)
	ggsave(plot = probPlot_truncated,
		   path = path,
		   filename = plotname)

}

# Ankit Roy
# 31st January, 2022