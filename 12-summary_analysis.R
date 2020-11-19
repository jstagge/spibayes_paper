





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
###  Compare confidence intervals for different reference periods
###########################################################################
### Need to do this



###########################################################################
###  Calculate K-S
###########################################################################

### Combine Estimate dataframe
estimate_nonstat <- mle_synth_non_stat_long$estimate %>%
	mutate(ref = "full") %>%
	bind_rows(mle_synth_non_stat_ref$estimate %>%
	mutate(ref = "ref"))

### Combine Draws dataframe
draws_nonstat <- mle_synth_non_stat_long$draws %>%
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



### Save results for next step
save(ci_comparison_summary, ks_summary, ks_summary_jdate, file = file.path(write_output_path, "mle_stat_summary.RData"))
save(ks_model, ks_model, file = file.path(write_output_path, "mle_stat_all.RData"))

### For CSV
write.csv(ci_comparison_summary, file = file.path(write_output_path, "mle_stat_ci_summary.csv"), row.names = FALSE)
write.csv(ks_summary, file = file.path(write_output_path, "mle_stat_ks_summary.csv"), row.names = FALSE)




