# *------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: .R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  
# | 
# |
# *------------------------------------------------------------------


###########################################################################
###  Load functions
###########################################################################
require(tidyverse)
require(here)

### For custom MLE functions
require(spibayes)

### To save in SVG
require(svglite)
require(viridis)

### Packages for spi
require(fitdistrplus)
require(lubridate)

select <- dplyr::select

theme_set(theme_classic(8))

###########################################################################
## Set the Paths
###########################################################################
### Set here path 
here_path <- here::here()

### Path for Data and Output	
data_path <- file.path(here_path, "./data")
output_path <- file.path(here_path, "./output")

### Set up output folders
write_output_path <- file.path(output_path, "mle_fit")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(output_path, "figures/mle_fit")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)


###########################################################################
## Load Functions
###########################################################################

file.sources = list.files(file.path(here_path, "./functions"), pattern="*.R", recursive=TRUE)
sapply(file.path(file.path(here_path, "./functions"), file.sources),source)

###########################################################################
## Set initial values
###########################################################################
### Set seed so everyone gets the same random numbers (reproducible example)
set.seed(7890)


###########################################################################
## Load synthetic precipitation time series
###########################################################################
load(file.path(output_path, "synth_precip/true_param.RData"))
load(file.path(output_path, "synth_precip/synth_df.RData"))

ls()


###########################################################################
###  Fit nonstationary data using MLE
###########################################################################

mle_synth_non_stat_full <- mle_gamma_fit(jdate = synth_non_stat_df$jdate, values = synth_non_stat_df$precip, n = 10000)

ref_df <- synth_non_stat_df %>%
	filter(date >= as.Date("1961-01-01")  & date <= as.Date("1990-12-31"))

mle_synth_non_stat_ref <- mle_gamma_fit(jdate = ref_df$jdate, values = ref_df$precip, n = 10000)


###########################################################################
###  Plot the results
###########################################################################
true_plot <- true_param_non_stat %>% 
	select(-month_day, -spi) %>%
	gather("param", "value", c(-jdate,-date, -year) ) %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta")))

### Plot the MLE result
mle_plot <- mle_synth_non_stat_full$estimate %>%
	pivot_wider(names_from = ci, values_from = value) %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta")))

p <- ggplot(mle_plot %>% filter(param != "Dispersion" & param != "Rate"), aes(x=jdate)) %>%
	+ geom_ribbon(aes(ymin = lower_ci, ymax=upper_ci), fill="grey50", colour="grey50", alpha=0.5) %>%
	+ geom_line(data = true_plot %>% filter(param != "Dispersion" & param != "Rate"), aes(y=value, colour= year, group = year), alpha = 0.2) %>%
	+ geom_line(colour="red", aes(y=estimate))  %>%
	+ scale_colour_viridis(name = "Year", option="inferno") %>%
	+ facet_grid(param ~ ., scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Parameter Value") %>%
	+ coord_cartesian(xlim=c(1, 365))

p
### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_non_stat_subset.png"), p, width =5, height = 7, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_non_stat_subset.pdf"), p, width =5, height = 7)
ggsave(file.path(write_figures_path, "mle_fit_non_stat_subset.svg"), p, width =5, height = 7)


###########################################################################
###  Compare confidence intervals for different reference periods
###########################################################################
### Need to do this



###########################################################################
###  Calculate K-S
###########################################################################

### Combine Estimate dataframe
estimate_nonstat <- mle_synth_non_stat_full$estimate %>%
	mutate(ref = "full") %>%
	bind_rows(mle_synth_non_stat_ref$estimate %>%
	mutate(ref = "ref"))

### Combine Draws dataframe
draws_nonstat <- mle_synth_non_stat_full$draws %>%
	mutate(ref = "full") %>%
	bind_rows(mle_synth_non_stat_ref$draws %>%
	mutate(ref = "ref"))

### At some point, I should make this a function
model_list <-  unique(estimate_nonstat$ref)

for(j in seq(1, length(model_list))){

model_j <- model_list[j]

estimate_j <- estimate_nonstat %>% 
	filter(ref == model_j)

draws_j <- draws_nonstat %>% 
	filter(ref == model_j)

for (i in seq(1, 365)){
cat(i)
cat("\n")

	### Create True CDF
	true_cdf  <- true_param_non_stat %>%
		filter(jdate == i)

	### Extract just the theta values from the estimate (no draws) and then combine with simulated draws
	theta_temp <- estimate_j %>% 
		filter(ci == "estimate") %>% 
		pivot_wider(names_from = "param", values_from = "value") %>%
		select(jdate, theta)	

	sim_cdf <- draws_j %>% 
		filter(jdate == i) %>%
		left_join(theta_temp, by = "jdate")

	### Create a common sequence of precipitation to use for both
	max_precip <- quantile(qgamma(0.999, scale = sim_cdf$scale, shape = sim_cdf$shape), 0.75, na.rm=TRUE, warnings = FALSE)
	precip_seq <- seq(0,max_precip, length.out = 600)

	### Calculate cumulative probability for each precipitation
	true_cdf <- true_cdf %>%
		left_join(data.frame(jdate = i, precip = precip_seq), by = "jdate") %>%
		filter(precip > 0)  %>% 
		mutate(p_pos = pgamma(precip, shape=shape, scale = scale)) %>%
		mutate(scale_factor = 1/(1-theta)) %>%
		mutate(p_overall = p_pos / scale_factor + theta)

	### Cut out a chunk of draws
	sim_cdf <- sim_cdf %>%
		left_join(data.frame(jdate = i, precip = precip_seq), by = "jdate") %>%
		filter(precip > 0)  %>% 
		filter(draw <= 500) %>%
		mutate(p_pos = pgamma(precip, shape=shape, scale = scale)) %>%
		mutate(scale_factor = 1/(1-theta)) %>%
		mutate(p_overall = p_pos / scale_factor + theta)

	### Combine simulated and true distributions by the common precipitation
	comb_df <- sim_cdf %>%
		select(jdate, draw, precip, p_overall) %>%
		rename(p_sim = p_overall) %>%
		full_join(true_cdf %>% select(jdate, date, precip, p_overall), by = c("jdate", "precip")) %>%
		rename(p_true = p_overall) %>%
		mutate(diff = p_sim - p_true) %>%
		mutate(abs_diff = abs(diff)) %>%
		drop_na(diff, abs_diff)

	### Calculate K-S metric and a few other statistics related to it
	ks_temp <- comb_df %>%
		group_by(jdate, date, draw) %>%
		summarise(ks_metric = max(abs_diff), ks_precip = precip[which.max(abs_diff)], ks_p_true = p_true[which.max(abs_diff)], ks_p_sim = p_sim[which.max(abs_diff)])

	if (i == 1) {
		ks_df <- ks_temp
	} else {
		ks_df <- ks_df %>%
			bind_rows(ks_temp)
	}

}

if(j == 1){
	ks_df <- ks_df %>%
		mutate(ref = model_j)
	ks_model <- ks_df
} else {
	ks_df <- ks_df %>%
		mutate(ref = model_j)
	ks_model <- ks_model %>%
		bind_rows(ks_df)
}

}

ks_model <- ks_model %>% 
	mutate(year = year(date))

ggplot(ks_model, aes(x=jdate, y=ks_metric, group = jdate)) + geom_boxplot() + facet_grid(ref ~ .)

ggplot(ks_model, aes(x=jdate, y=ks_metric)) + geom_point() + facet_grid(ref ~ jdate)

ggplot(ks_model, aes(x=date, y=ks_metric)) + geom_point() + facet_grid(ref ~ jdate)

ggplot(ks_model %>% filter(jdate == 50), aes(x=date, y=ks_metric)) + geom_point() + facet_grid(ref ~ .)

ggplot(ks_model %>% filter(jdate == 50), aes(x=year, y=ks_metric, group = year)) + geom_boxplot() + facet_grid(ref ~ .)


ggplot(ks_model %>% filter(jdate == 150), aes(x=date, y=ks_metric)) + geom_point() + facet_grid(ref ~ .)


ks_summary_byyear <- ks_model %>% 
	mutate(year = year(date)) %>%
	group_by(ref, year) %>%
	summarise(ks_median = median(ks_metric), ks_95_low = quantile(ks_metric, 0.025), ks_95_high = quantile(ks_metric, 0.975))

p <- ggplot(ks_summary_byyear, aes(x=year)) + geom_ribbon(data = ks_summary_byyear, aes(ymin =ks_95_low, ymax = ks_95_high, fill=ref), alpha = 0.1) + geom_line(data = ks_summary_byyear, aes(y=ks_median, colour = ref)) + geom_vline(xintercept = c(1961, 1991), linetype = "longdash")

p

### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_nonstat_ks_byyear.png"), p, width =8, height = 6, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_nonstat_ks_byyear.pdf"), p, width =8, height = 6)
ggsave(file.path(write_figures_path, "mle_fit_nonstat_ks_byyear.svg"), p, width =8, height = 6)


ks_summary <- ks_model %>% 
	mutate(year = year(date)) %>%
	group_by(ref) %>%
	summarise(ks_median = median(ks_metric), ks_95_low = quantile(ks_metric, 0.025), ks_95_high = quantile(ks_metric, 0.975))

###########################################################################
###  Save results from as an RDS file for analysis
###########################################################################

### Save results for next step
save(ks_summary_byyear, ks_summary, file = file.path(write_output_path, "ks_nonstat.RData"))

### For CSV
write.csv(ks_summary, file = file.path(write_output_path, "mle_nonstat_ks_summary.csv"), row.names = FALSE)
write.csv(ks_summary_byyear, file = file.path(write_output_path, "mle_nonstat_ks_summary_byyear.csv"), row.names = FALSE)




