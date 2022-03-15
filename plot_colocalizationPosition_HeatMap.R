library(ggplot2)
library(dplyr)
library(forcats)

path = "/Users/roy/Sandbox/headtMap_2022.02.11/all_files/"			# file path

# all plottable data files
files = list.files(path = path,
				   pattern = "*_heatPlotData.csv")

for (filename in files) {
	
	# plot file name
	plotname_norm_start = paste(substr(filename, 0, nchar(filename)-4), "_normStart.svg", sep = "")
	plotname_norm_end = paste(substr(filename, 0, nchar(filename)-4), "_normEnd.svg", sep = "")
	
	# colocalization heat map data
	data = read.csv(paste(path, filename, sep = ""))
	data = data %>% subset(CHANNEL == "GTPase")
	data = data %>% mutate(name = fct_reorder(PSEUDO_TRACK_ID, desc(TRACK_LENGTH)))
	
	# plot with normalized start positions
	start_norm_plot = ggplot(data, aes(NORM_START, name, fill = COLOCALIZATION_STATUS)) +
		geom_tile() +
		scale_x_continuous(limits = c(-1, 100)) +
		scale_fill_manual(values = c("#4056F4", "#B8B42D", "#FFFCE8", "#3E363F")) +
		labs(fill = "Colocalisation Status") +
		theme_bw() +
		theme(axis.ticks.y = element_blank(),
			  axis.text.y = element_blank(),
			  legend.position = "bottom") +
		facet_wrap(vars(ANNOTATION_TRACK)) +
		ylab("GTPase Track") +
		xlab("Track Lifetime (frames)") +
		ggtitle(filename) +
		coord_fixed()
	
	# plot with normalized end positions
	end_norm_plot = ggplot(data, aes(NORM_END, name, fill = COLOCALIZATION_STATUS)) +
		geom_tile() +
		scale_x_continuous(limits = c(-100, 1)) +
		scale_fill_manual(values = c("#4056F4", "#B8B42D", "#FFFCE8", "#3E363F")) +
		labs(fill = "Colocalisation Status") +
		theme_bw() +
		theme(axis.ticks.y = element_blank(),
			  axis.text.y = element_blank(),
			  legend.position = "bottom") +
		facet_wrap(vars(ANNOTATION_TRACK)) +
		ylab("GTPase Track") +
		xlab("Track Lifetime (frames)") +
		ggtitle(filename) +
		coord_fixed()
	
	ggsave(filename = paste(path, plotname_norm_start, sep = ""), plot = start_norm_plot)
	ggsave(filename = paste(path, plotname_norm_end, sep = ""), plot = end_norm_plot)

}

# Ankit Roy
# 1st March, 2022
# 2nd March, 2022	>>	Fixed coordinates for plots.