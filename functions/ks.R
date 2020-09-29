
spi_cdf <- function(shape, scale, theta, p_seq = seq(0.0001,0.9999,0.0001)){

	cdf_df <- data.frame(p_overall = p_seq, shape = shape, scale = scale, theta = theta) %>%
		mutate(precip_notheta = qgamma(p_overall, shape = shape, scale = scale))

	cdf_df <- cdf_df %>% 
		mutate(spi_thresh = qnorm(theta)) %>%
		mutate(precip_zero = case_when(
			p_overall < theta ~ 0,
			TRUE ~ NA_real_)) %>%
		mutate(p_precip =  case_when(
			is.na(precip_zero) ~ p_overall - theta,
			TRUE ~ NA_real_)) %>%
		mutate(scale_factor = 1/(1-theta)) %>%
		mutate(p_precip = p_precip*scale_factor) %>%
		mutate(precip = case_when(
			precip_zero == 0 ~ 0,
			is.na(precip_zero) ~ qgamma(p_precip, shape = shape, scale = scale),
			TRUE ~ NA_real_)
		) 
	return(cdf_df)
}



ks_metric <- function(sim, true){

	### Calculate the CDF based on the true parameters
	true_cdf <- spi_cdf(shape = true$shape, scale = true$scale, theta=true$theta)

	true_zero <- true_cdf[1:2,] %>%
		mutate(precip = c(0,0), p_overall = c(0,true$theta))

	true_cdf <- true_zero %>%
		bind_rows(true_cdf %>% filter(precip > 0))

	### Calculate the CDF based on the simulated parameters
	### This needs to be done in a slightly different way so that the precip values are equivalent to the true
	### It makes calculating the K-S distance much easier
	sim_cdf <- true_cdf %>%
		select(precip ) %>%
		mutate(shape = sim$shape, scale = sim$scale, theta=sim$theta)

	sim_pos <- sim_cdf  %>%
		filter(precip > 0)  %>% 
		mutate(p_pos = pgamma(precip, shape=shape, scale = scale)) %>%
		mutate(scale_factor = 1/(1-theta)) %>%
		mutate(p_overall = p_pos / scale_factor + theta)

	sim_zero <- sim_pos[1:2,] %>%
		mutate(precip = c(0,0), p_overall = c(0,sim$theta))

 	sim_cdf <- sim_zero %>%
		bind_rows(sim_pos)

	### Calculate the K-S distance
	sim_cdf <- sim_cdf %>%
		mutate(true_p_overall = true_cdf$p_overall) %>%
		mutate(diff = p_overall - true_p_overall) %>%
		mutate(abs_diff = abs(diff))

	ks_metric <- max(sim_cdf$abs_diff)

	ks_loc <- which.max(sim_cdf$abs_diff)
	ks_p_true <- sim_cdf$true_p_overall[ks_loc]
	
	return(list(ks_metric = ks_metric, ks_p_sim = ks_p_true, sim_cdf = sim_cdf))

}

