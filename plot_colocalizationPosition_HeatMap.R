library(ggplot2)
library(dplyr)
library(forcats)

data = read.csv("/Users/roy/Sandbox/headtMap_2022.02.11/hello.csv")

data = data %>% subset(CHANNEL == "GTPase")

data = data %>% mutate(name = fct_reorder(PSEUDO_TRACK_ID, desc(TRACK_LENGTH)))

start_norm_plot = ggplot(data, aes(NORM_START, name, fill = COLOCALIZATION_STATUS)) + geom_tile() + scale_x_continuous(limits = c(-1, 100)) + scale_fill_manual(values = c("#4056F4", "#B8B42D", "#FFFCE8", "#3E363F")) + labs(fill = "Colocalisation Status") + theme_bw() + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + facet_wrap(vars(ANNOTATION_TRACK)) + ylab("GTPase Track") + xlab("Track Lifetime (frames)")

end_norm_plot = ggplot(data, aes(NORM_END, name, fill = COLOCALIZATION_STATUS)) + geom_tile() + scale_x_continuous(limits = c(-100, 1)) + scale_fill_manual(values = c("#4056F4", "#B8B42D", "#FFFCE8", "#3E363F")) + labs(fill = "Colocalisation Status") + theme_bw() + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") + facet_wrap(vars(ANNOTATION_TRACK)) + ylab("GTPase Track") + xlab("Track Lifetime (frames)")

start_norm_plot
end_norm_plot

# Ankit Roy
# 1st March, 2022