library(dplyr)
library(ggplot2)

data = read.csv("/Users/roy/Sandbox/Redo_26.01.22/PC_GEF_100_pm_complex_GTP_1_colocalisation.csv", comment.char = "#")

dt = 0.022
protein = "GTPase"

#--- sub-setting data for dwell time calculations
data_subset_dwell = data %>%
	subset(CHANNEL == protein) %>%
	group_by(PSEUDO_TRACK_ID) %>%
	summarise(PSEUDO_TRACK_ID,
			  FRAME_COUNT = max(FRAME, na.rm = TRUE) - min(FRAME, na.rm = TRUE) + 1)

data_subset_dwell$DWELL_TIME = data_subset_dwell$FRAME_COUNT * dt
data_subset_dwell = unique(data_subset_dwell)

trackFreq_vs_dwellTime = ggplot(data_subset_dwell, aes(DWELL_TIME)) +
	geom_histogram(binwidth = dt, color = 'black') +
	scale_x_continuous(limits = c(0, 5)) +
	xlab("Dwell Time (s)") +
	ylab("Frequency")

trackFreq_vs_dwellTime

#--- sub-setting data for empirical CDF calculations
ecdf_fit = ecdf(data_subset_dwell$DWELL_TIME)

data_subset_ecdf = data_subset_dwell %>%
	group_by(DWELL_TIME) %>%
	summarise(DWELL_TIME,
			  CDF = ecdf_fit(DWELL_TIME))

data_subset_ecdf = unique(data_subset_ecdf)

cdf_plot = ggplot(data_subset_ecdf, aes(DWELL_TIME, CDF)) +
	geom_point() +
	geom_step() +
	scale_y_continuous(limits = c(0, 1)) +
	xlab("Time (s)") +
	ylab("CDF")

cdf_plot

# 1 - CDF
data_subset_ecdf$CDF_inv = 1 - data_subset_ecdf$CDF

cdf_inv_plot = ggplot(data_subset_ecdf, aes(DWELL_TIME, CDF_inv)) +
	geom_point() +
	geom_step() +
	scale_y_continuous(limits = c(0, 1)) +
	xlab("Time (s)") +
	ylab("1-CDF")

cdf_inv_plot

# log(1 - CDF)
log_cdf_inv_plot = ggplot(data_subset_ecdf, aes(DWELL_TIME, CDF_inv)) +
	geom_point() +
	scale_y_continuous(trans = "log10", limits = c(0.001,1)) +
	scale_x_continuous(breaks = seq(0, 15, 2)) +
	xlab("Time (s)") +
	ylab(bquote(log[10](1-CDF)))

log_cdf_inv_plot

#--- sub-setting data for step-size calculations
data_subset_step_size = data %>%
	subset(CHANNEL == protein) %>%
	group_by(PSEUDO_TRACK_ID) %>%
	summarise(PSEUDO_TRACK_ID,
			  FRAME,
			  dX = (POSITION_X - lag(POSITION_X))^2,
			  dY = (POSITION_Y - lag(POSITION_Y))^2)

data_subset_step_size$STEP_SIZE = (data_subset_step_size$dX + data_subset_step_size$dY)^0.5
data_subset_step_size = unique(data_subset_step_size)

step_size_freq = ggplot(data_subset_step_size, aes(STEP_SIZE)) +
	geom_histogram(binwidth = 0.05, colour = 'black') +
	scale_x_continuous(breaks = seq(0, 3, 0.2)) +
	xlab("Step Size (µm)") +
	ylab("Frequency")

step_size_freq

#--- sub-setting data for MSD calculations
data_subset_msd = data %>%
	subset(CHANNEL == protein) %>%
	group_by(PSEUDO_TRACK_ID) %>%
	summarise(PSEUDO_TRACK_ID,
			  FRAME,
			  dX = (POSITION_X - lag(POSITION_X))^2,
			  dY = (POSITION_Y - lag(POSITION_Y))^2)

data_subset_msd$d2 = data_subset_msd$dX + data_subset_msd$dY

data_subset_msd = data_subset_msd %>%
	group_by(PSEUDO_TRACK_ID) %>%
	summarise(PSEUDO_TRACK_ID,
			  FRAME_COUNT = max(FRAME, na.rm = TRUE) - min(FRAME, na.rm = TRUE) + 1,
			  d2_sum = sum(d2, na.rm = TRUE))

data_subset_msd = unique(data_subset_msd)

data_subset_msd = data_subset_msd %>%
	group_by(FRAME_COUNT) %>%
	summarise(FRAME_COUNT,
			  MSD = mean(d2_sum),
			  STD_MSD = sd(d2_sum),
			  DWELL = FRAME_COUNT * dt)

data_subset_msd = unique(data_subset_msd)

msd_time_evolution = ggplot(data_subset_msd, aes(DWELL, MSD)) +
	geom_point() +
	geom_errorbar(aes(ymin = MSD - STD_MSD, ymax = MSD + STD_MSD)) +
	scale_x_continuous(limits = c(0, 2.5)) +
	scale_y_continuous(limits = c(0, 40)) +
	xlab("Dwell Time (s)") +
	ylab(bquote(MSD~(µm^2)))

msd_time_evolution

# Ankit Roy
# 9th November, 2021