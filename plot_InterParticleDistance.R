library(ggplot2)
library(gganimate)

infile = "/Users/roy/Sandbox/bug_fix/PLL_Glass_crosslinked_1001_42_PB_All_Spot_statistics_dist.csv"
plot_title = "GTPase"
static_plot_outfile = "interparticle_distance_frequency.jpeg"
animated_plot_outfile = "interparticle_distance_frequency.mp4"

# Interparticle distance data
data = read.csv(infile)

# Interparticle distance static plot
ipd_plot = ggplot(data, aes(DIST)) +
	geom_histogram(binwidth = 0.1, alpha = 0.5) +
	scale_x_continuous(limits = c(0, 15),
					   breaks = seq(0, 15, 2)) +
	geom_vline(xintercept = 1,
			   linetype = 'dashed',
			   size = 0.5) +
	xlab("Minimum inter-particle distance (µm)") +
	ylab("Spot frequency") +
	ggtitle(plot_title)

# Interparticle distance plot animated over frames
ipd_plot_anim = ipd_plot +
	transition_manual(FRAME) +
	labs(subtitle = "FRAME: {current_frame}")

# Animate interparticle distance plot
anim = animate(ipd_plot_anim,
			   nframes = 500,
			   renderer = av_renderer(),
			   height = 5,
			   width = 5,
			   units = "in",
			   res = 150)

# Save animated interparticle distance plot
anim_save(animated_plot_outfile, anim)

# Ankit Roy
# 8th December, 2021