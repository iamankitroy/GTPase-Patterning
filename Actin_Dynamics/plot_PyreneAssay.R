#/usr/local/bin/R

#--- This script is used to plot Pyrene Assay fluorometry data
#--- The script requires a CSV file with time in one column and fluorescence intensities in others

library(reshape2)
library(ggplot2)
library(dplyr)
library(drc)
library(aomisc)

fluordata = read.csv("/Users/roy/Documents/Project/PyreneAssays/May2020/Profilin-titration_07.05.2020.csv")

# Convert time to minutes
fluordata$Time = fluordata$Time/60

# Long table
fluordata_long = melt(fluordata, id.vars = "Time", variable.name = "Well", value.name = "Fluorescence")

# Map wells to concentrations
fluordata_long = fluordata_long %>%
  mutate(Conc = case_when(Well == "A1" | Well == "A7" ~ 0,
                          Well == "A2" | Well == "A8" ~ 2,
                          Well == "A3" | Well == "A9" ~ 3,
                          Well == "A4" | Well == "A10" ~ 4,
                          Well == "A5" | Well == "A11" ~ 5,
                          Well == "A6" | Well == "A12" ~ 6))

# Convert concentration to factor
fluordata_long$Conc = factor(fluordata_long$Conc)

#--- Get maximum fluorescence for normalization
#--- Inputs: <fluorescence-data>, <reference-conc>, <saturation-start-time>, <saturation-end-time>
normalize_fluor = function(data, ref_conc, ref_timeStart, ref_timeEnd){
  
  # Subset of points with saturated fluorescence intensities
  ref_data = data %>%
    filter(Conc == ref_conc &
             Time > ref_timeStart &
             Time <= ref_timeEnd)
  
  # Get average initial fluorescence intensity
  avg_min_fluor = mean(filter(data, Time == 0)$Fluorescence)
  
  # Get average of saturated fluorescence intensities
  ref_data$Fluorescence = ref_data$Fluorescence - avg_min_fluor
  avg_max_fluor = mean(ref_data$Fluorescence)
  
  # Normalize
  norm_data = data
  norm_data$Fluorescence = (data$Fluorescence - avg_min_fluor) / avg_max_fluor
  
  return(norm_data)
}

norm_fluordata = normalize_fluor(fluordata_long, 0, 25, 60)

# Logistic growth curve fit
model = drm(Fluorescence ~ Time,
            curveid = Conc,
            fct = L.4(),
            data = norm_fluordata)

# Predict fit from model
norm_fluordata$pred = predict(model)

fluorplot = ggplot(norm_fluordata, aes(Time, Fluorescence, color = Conc)) +
  geom_point(alpha = 0.25, size = 0.25) +
  geom_line(aes(y = pred), size = 1) +
  ggtitle("5µM Actin +") +
  xlab("Time (min)") +
  ylab("Normalized Fluorescence Intensity") +
  scale_x_continuous(breaks = seq(0, 60, 10)) +
  scale_y_continuous(expand = c(0.025, 0.025)) +
  guides(col = guide_legend(title = "Profilin (µM)",
                            nrow = 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 20),
        legend.title = element_text(face = "bold",
                                    size = 18),
        legend.position = "top",
        legend.box = "horizontal",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16))

fluorplot

ggsave("/Users/roy/Documents/Project/PyreneAssays/May2020/Profilin-titration_07.05.2020.jpeg",
       plot = fluorplot,
       height = 10,
       width = 12,
       units = "in",
       dpi = 300)

# Ankit Roy
# 13th September, 2020
