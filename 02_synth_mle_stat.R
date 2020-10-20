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
###  Fit synthetic data using MLE
###########################################################################

## Use fitdist (MLE)
### Fit each day with MLE
mle_synth_stat <- mle_gamma_fit(jdate = synth_stat_df$jdate, values = synth_stat_df$precip, n = 3000, method = "param", ci = 0.95)

###########################################################################
###  Plot the results
###########################################################################
true_plot <- true_param_stat %>% 
	gather("param", "value", -jdate) %>%
	mutate(origin = "true") %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta")))

### Plot the MLE result
mle_plot <- mle_synth_stat$estimate %>%
	pivot_wider(names_from = ci, values_from = value) %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta")))

p <- ggplot(mle_plot, aes(x=jdate)) %>%
	+ geom_ribbon(aes(ymin = lower_ci, ymax=upper_ci), fill="grey50", colour="grey50", alpha=0.5) %>%
	+ geom_line(colour="red", aes(y=estimate))  %>%
	+ geom_line(data = true_plot, aes(y=value), colour = "black") %>%
	+ facet_grid(param ~ ., scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Parameter Value") #%>%
#	+ coord_cartesian(xlim=c(1, 365))

p

### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_full.png"), p, width =5, height = 7, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_full.pdf"), p, width =5, height = 7)
ggsave(file.path(write_figures_path, "mle_fit_stat_full.svg"), p, width =5, height = 7)

p <- ggplot(mle_plot %>% filter(param != "Dispersion" & param != "Rate"), aes(x=jdate)) %>%
	+ geom_ribbon(aes(ymin = lower_ci, ymax=upper_ci), fill="grey50", colour="grey50", alpha=0.5) %>%
	+ geom_line(colour="red", aes(y=estimate))  %>%
	+ geom_line(data = true_plot %>% filter(param != "Dispersion" & param != "Rate"), aes(y=value), colour = "black") %>%
	+ facet_grid(param ~ ., scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Parameter Value") %>%
	+ coord_cartesian(xlim=c(1, 365))

p
### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_full_subset.png"), p, width =5, height = 7, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_full_subset.pdf"), p, width =5, height = 7)
ggsave(file.path(write_figures_path, "mle_fit_stat_full_subset.svg"), p, width =5, height = 7)


###########################################################################
###  Fit synthetic data using MLE with different data lengths
###########################################################################

start_list <- c(1945, 1970, 1990)
year_list <- c(75, 50, 30)

### Save the model for later with reference column
mle_model_df <- mle_synth_stat

mle_model_df$estimate <- mle_model_df$estimate %>%
	mutate(ref = "stationary_100_years")

mle_model_df$draws <- mle_model_df$draws %>%
	mutate(ref = "stationary_100_years")


for (j in seq(1,3)){
	short_df <- synth_stat_df %>% 
		filter(date > as.Date(paste0(start_list[j], "-01-01")))
	
	### Assume the first 92 days rolling mean creates 91 days of NA
	short_df$precip[1:91] <- NA

	### Fit model
	model_temp <- mle_gamma_fit(jdate = short_df$jdate, values = short_df$precip, n = 3000)

	### Add reference column
	mle_temp <- model_temp
	mle_temp$estimate <- mle_temp$estimate %>%
		mutate(ref = paste0("stationary_", year_list[j], "_years"))

	mle_temp$draws <- mle_temp$draws %>%
		mutate(ref = paste0("stationary_", year_list[j], "_years"))

	### Bind rows
	mle_model_df$estimate <- mle_model_df$estimate %>%
		bind_rows(mle_temp$estimate)	

	mle_model_df$draws <- mle_model_df$draws %>%
		bind_rows(mle_temp$draws)	
}


### Plot the MLE result
mle_plot <- mle_model_df$estimate %>%
	pivot_wider(names_from = ci, values_from = value) %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta"))) %>%
	mutate(ref = factor(ref, levels = c("stationary_30_years", "stationary_50_years", "stationary_75_years", "stationary_100_years"), labels = c("30 Years", "50 Years", "75 Years", "100 Years")))

p <- ggplot(mle_plot %>% filter(param == "Mean" | param == "Scale" | param == "Shape"), aes(x=jdate)) %>%
	+ geom_ribbon(aes(ymin = lower_ci, ymax=upper_ci, fill = ref), colour ="grey40") %>%
#	+ geom_line(colour="red", aes(y=estimate))  %>%
	+ geom_line(data = true_plot %>% filter(param == "Mean" | param == "Scale" | param == "Shape"), aes(y=value), colour = "black", size =1) %>%
	+ facet_grid(param ~ ., scales="free_y") %>%
	+ theme_bw(8) %>%
	#+ scale_colour_brewer(type="seq", palette = "OrRd", direction = -1) %>%
	+ scale_fill_brewer(name = "Date", type="seq", palette = "OrRd", direction = -1) %>%
	+ theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Parameter Value") #%>%
#	+ coord_cartesian(xlim=c(1, 365))

p

### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_ci_ts.png"), p, width =6, height = 6, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_ci_ts.pdf"), p, width =6, height = 6)
ggsave(file.path(write_figures_path, "mle_fit_stat_ci_ts.svg"), p, width =6, height = 6)

###########################################################################
###  Calculate confidence interval width
###########################################################################
### Reorganize to put the confidence intervals in columns and then calculate the confidence width
ci_comparison <- mle_model_df$estimate %>%
	pivot_wider(names_from = ci, values_from = value) %>%
	mutate(ci_width = upper_ci - lower_ci) %>% 
	mutate(ref = factor(ref, levels = c("stationary_30_years", "stationary_50_years", "stationary_75_years", "stationary_100_years"), labels = c("30 Years", "50 Years", "75 Years", "100 Years")))

head(ci_comparison)

### A couple plots for reference
ggplot(ci_comparison, aes(x=ci_width, colour=ref)) + geom_density(aes(fill=ref), alpha = 0.5) + facet_wrap(.~param, scales = "free") + theme_bw()

plot_df <- ci_comparison %>% 
	filter(param == "mean" | param == "scale" | param == "shape") %>%
	mutate(param = factor(param, levels = c("mean", "scale", "shape"), labels = c("Mean", "Scale", "Shape")))

p <- ggplot(plot_df, aes(x=ref, y=ci_width)) %>%
	+ geom_boxplot(aes(fill=ref)) %>%
	+ facet_wrap(.~param, scales = "free_y") %>%
	+ theme_bw() %>%
	+ theme(legend.position = "none") %>%
	+ scale_y_continuous(name = "Confidence Interval Width")
p
### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_ci_density.png"), p, width =6, height = 4, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_ci_density.pdf"), p, width =6, height = 4)
ggsave(file.path(write_figures_path, "mle_fit_stat_ci_density.svg"), p, width =6, height = 4)

ci_comparison <- ci_comparison %>%
	mutate(ci_relative = ci_width/estimate)

p <- ggplot(ci_comparison %>% filter(param == "mean" | param == "scale" | param == "shape"), aes(x=ref, y=ci_relative)) + geom_boxplot(aes(fill=ref)) + facet_wrap(.~param) + theme_bw()

### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_ci_box.png"), p, width =6, height = 4, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_ci_box.pdf"), p, width =6, height = 4)
ggsave(file.path(write_figures_path, "mle_fit_stat_ci_box.svg"), p, width =6, height = 4)


ci_comparison_summary <- ci_comparison %>%
	group_by(ref, param) %>%
	summarize(ci_width = median(ci_width)) %>%
	pivot_wider(names_from = param, values_from = ci_width) %>%
	select(ref, mean, scale, shape) %>%
	as.data.frame()

ci_comparison_summary


###########################################################################
###  Calculate RMSE
###########################################################################
stat_year <- synth_stat_df %>%
	filter(year == 1975) %>% 
	select(jdate, mean, shape, scale, rate, theta, disp) %>%
	pivot_longer(cols=c(-jdate), names_to = "param", values_to = "values")

rmse_comparison <- mle_model_df$draws %>%
	select(-type) %>%
	pivot_longer(cols=c(-jdate, -ref, -draw), names_to = "param", values_to = "values") %>%
	left_join(stat_year, by = c("jdate", "param"))


rmse_comparison <- rmse_comparison %>%
	mutate(error = values.y - values.x) %>%
	mutate(error_sq = error^2)

sna <- rmse_comparison %>%
	group_by(ref, param) %>%
	summarize(rmse = sqrt(mean(error_sq, na.rm=TRUE)))

sna <- sna %>%
	pivot_wider(names_from = ref, values_from = rmse)

sna

#crps_sample()
#CRPS generalizes the mean absolute error




###########################################################################
###  Calculate K-S
###########################################################################
### At some point, I should make this a function
model_list <-  unique(mle_model_df$draws$ref)

for(j in seq(1, length(model_list))){

model_j <- model_list[j]

estimate_j <- mle_model_df$estimate %>% 
	filter(ref == model_j)

draws_j <- mle_model_df$draws %>% 
	filter(ref == model_j)

for (i in seq(1, 365)){
cat(i)
cat("\n")

	### Create True CDF
	true_cdf  <- true_param_stat %>%
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

	sim_cdf <- sim_cdf %>%
		left_join(data.frame(jdate = i, precip = precip_seq), by = "jdate") %>%
		filter(precip > 0)  %>% 
		mutate(p_pos = pgamma(precip, shape=shape, scale = scale)) %>%
		mutate(scale_factor = 1/(1-theta)) %>%
		mutate(p_overall = p_pos / scale_factor + theta)

	### Combine simulated and true distributions by the common precipitation
	comb_df <- sim_cdf %>%
		select(jdate, draw, precip, p_overall) %>%
		rename(p_sim = p_overall) %>%
		left_join(true_cdf %>% select(jdate, precip, p_overall), by = c("jdate", "precip")) %>%
		rename(p_true = p_overall) %>%
		mutate(diff = p_sim - p_true) %>%
		mutate(abs_diff = abs(diff)) %>%
		drop_na(diff, abs_diff)

	### Calculate K-S metric and a few other statistics related to it
	ks_temp <- comb_df %>%
		group_by(jdate, draw) %>%
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

rm(ks_df)

ks_model$ref <- factor(ks_model$ref, levels = c("stationary_100_years", "stationary_75_years", "stationary_50_years", "stationary_30_years"))

ggplot(ks_model, aes(x=jdate, group=jdate)) + geom_boxplot(aes(y=ks_metric)) + facet_grid(ref ~ .)

ks_summary_jdate <- ks_model %>% 
	group_by(ref, jdate) %>%
	summarise(ks_median = median(ks_metric), ks_95_low = quantile(ks_metric, 0.025), ks_95_high = quantile(ks_metric, 0.975))

p <- ggplot(ks_model, aes(x=jdate)) + geom_jitter(aes(y=ks_metric), alpha =0.05) + geom_ribbon(data = ks_summary_jdate, aes(ymin =ks_95_low, ymax = ks_95_high), fill = "grey50") + geom_line(data = ks_summary_jdate, aes(y=ks_median), colour="red") + facet_grid(ref ~ .)

### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_facets.png"), p, width =6, height = 8, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_facets.pdf"), p, width =6, height = 8)
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_facets.svg"), p, width =6, height = 8)


### Couple more plots
plot_df <- ks_summary_jdate
plot_df$ref <- factor(plot_df$ref, levels = c("stationary_30_years", "stationary_50_years", "stationary_75_years", "stationary_100_years"))
plot_df <- plot_df %>% arrange(ref)

p <- ggplot(plot_df, aes(x=jdate)) + geom_ribbon(data = plot_df, aes(ymin =ks_95_low, ymax = ks_95_high, fill = ref), alpha=0.7) + geom_line(data = plot_df, aes(y=ks_median, colour = ref))
p

### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_overplot.png"), p, width =6, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_overplot.pdf"), p, width =6, height = 4.5)
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_overplot.svg"), p, width =6, height = 4.5)


p <- ggplot(ks_model, aes(x=ks_metric, fill = ref)) + geom_density(alpha = 0.8) + scale_fill_brewer(type="qual") + coord_cartesian(xlim=c(0,0.22))
p

### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_density.png"), p, width =5, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_density.pdf"), p, width =5, height = 4.5)
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_density.svg"), p, width =5, height = 4.5)


p <- ggplot(ks_model, aes(y=ks_metric, fill = ref, x = ref)) + geom_boxplot() + scale_fill_brewer(type="qual")
### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_box.png"), p, width =5, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_box.pdf"), p, width =5, height = 4.5)
ggsave(file.path(write_figures_path, "mle_fit_stat_ks_box.svg"), p, width =5, height = 4.5)


ks_summary <- ks_model %>% 
	group_by(ref) %>%
	summarise(ks_median = median(ks_metric), ks_95_low = quantile(ks_metric, 0.025), ks_95_high = quantile(ks_metric, 0.975))

ks_summary


###########################################################################
###  Save results from as an RDS file for analysis
###########################################################################

### Save results for next step
save(mle_synth_stat, file = file.path(write_output_path, "mle_synth_stat.rda"))

### Save results for next step
save(ci_comparison_summary, ks_summary, ks_summary_jdate, file = file.path(write_output_path, "mle_stat_summary.RData"))
save(ks_model, ks_model, file = file.path(write_output_path, "mle_stat_all.RData"))

### For CSV
write.csv(ci_comparison_summary, file = file.path(write_output_path, "mle_stat_ci_summary.csv"), row.names = FALSE)
write.csv(ks_summary, file = file.path(write_output_path, "mle_stat_ks_summary.csv"), row.names = FALSE)



