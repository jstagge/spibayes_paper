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

### To access GHCND
require(rnoaa)

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
write_output_base <- file.path(output_path, "synth_precip")
dir.create(write_output_base, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_base <- file.path(output_path, "figures/synth_precip")
dir.create(write_figures_base, recursive=TRUE, showWarnings = FALSE)


###########################################################################
## Set initial values
###########################################################################
### Set seed so everyone gets the same random numbers (reproducible example)
set.seed(7890)

### SPI-3 3 months equals 92 days
n_roll <- 92



###########################################################################
###  Create a true SPI time series
###########################################################################
### Assume 100 years from 1920 to 2019, will then take subsets of this
spi_true <- data.frame(date = seq(as.Date("1920-01-01"), as.Date("2019-12-31"), by = "1 day")) %>%
	mutate(jdate = yday(date)) %>%
	mutate(month_day = paste0(month(date),"-",day(date))) %>%
	mutate(year = year(date))

n_days <- dim(spi_true)[1]

### Use a 92 day moving average MA(91) with innovations of sqrt(92)/n_roll, which produces N(0,1)
### Technically, the moving average would be on precip, not SPI, but this is a very close approximation
innov_c <- rnorm(n_days, 0, sqrt(n_roll))

### Generate SPI series using an MA(91) model with coef of 1. Need to remove the first 91 values
sim_spi <- arima.sim(list(order = c(0,0,(n_roll - 1)), ma = rep(1,(n_roll -1))), n = n_days, innov=innov_c/n_roll)
sim_spi[seq(1,(n_roll-1))] <- NA

spi_true <- spi_true %>%
	mutate(innov = innov_c/n_roll) %>%
	mutate(spi = c(sim_spi))

ggplot(spi_true, aes(x=date, y=spi)) + geom_line() + theme_classic()

mean(spi_true$spi, na.rm=TRUE)
sd(spi_true$spi, na.rm=TRUE)

###########################################################################
###  Make up synthetic gamma distribution (stationary)
###########################################################################
### Make some daily synthetic distributions that vary smoothly through time
true_param_stat <- tibble(jdate = seq(1,365), shape = rep(1.7, 365)) %>%
	mutate(theta = case_when( jdate <=182 ~ 0.05 + (0.05/182)*jdate, 
		jdate > 182 ~ 0.1 - (0.05/183)*(jdate - 182))
	) %>%
	mutate(mean = sin(((2*pi)/365)*jdate)+5+0.4*sin(((10*pi)/365)*jdate)) %>%
	mutate(scale = mean / shape) %>%
	mutate(rate = 1/scale) %>%
	mutate(disp = 1/shape) %>%
	select(jdate, shape, scale, scale, rate, theta, mean, disp)


###########################################################################
###  Plot the underlying distribution
###########################################################################
write_output_path <- file.path(write_output_base, "stationary")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(write_figures_base, "stationary")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

### Plot how this looks
plot_df <- true_param_stat %>% 
	gather("param", "value", -jdate) %>%
	mutate(origin = "true") %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta")))

### Create the figure
p <- ggplot(plot_df, aes(x=jdate, y=value)) %>%
	+ geom_line() %>%
	+ facet_grid(param ~ ., scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "True Parameter Value") %>%
	+ coord_cartesian(xlim=c(1, 365))

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_true_params.png"), p, width =5, height = 7, dpi = 300)
ggsave(file.path(write_figures_path, "synth_true_params.pdf"), p, width =5, height = 7)
ggsave(file.path(write_figures_path, "synth_true_params.svg"), p, width =5, height = 7)



### Create the figure
p <- ggplot(plot_df %>% filter(param != "Dispersion" & param != "Rate"), aes(x=jdate, y=value)) %>%
	+ geom_line() %>%
	+ facet_grid(param ~ ., scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme( panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "True Parameter Value") %>%
	+ coord_cartesian(xlim=c(1, 365))

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_true_params_subset.png"), p, width =5, height = 7, dpi = 300)
ggsave(file.path(write_figures_path, "synth_true_params_subset.pdf"), p, width =5, height = 7)
ggsave(file.path(write_figures_path, "synth_true_params_subset.svg"), p, width =5, height = 7)


###########################################################################
###  Sample from this distribution
###########################################################################
### Apply day 365 to day 366, join with true parameters, and then return to the original format
synth_stat_df <- spi_true %>%
	mutate(jdate = case_when(
			jdate == 366 ~ 365,
			TRUE ~ jdate)
	) %>%
	left_join(true_param_stat) %>%
	mutate(jdate = yday(date)) 

### First transform SPI into percentile. Then set zeros when percentile is below theta threshold
### Stretch the remaining probability space so that it fills the space (0 to 1)
### Then sample from gamma distribution and combine with zeros
synth_stat_df <- synth_stat_df %>% 
	mutate(p_overall = pnorm(spi)) %>%
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
	) %>%
	select(-p_overall, precip_zero, p_precip, scale_factor) %>%
	mutate(zero = precip == 0)



###########################################################################
###  Plot the stationary synthetic distribution
###########################################################################
p <- ggplot(synth_stat_df, aes(x=jdate, y=precip, colour=zero, group = year)) %>%
#	+ geom_jitter(alpha=0.2) %>%
	+ geom_line(alpha = 0.2) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(1,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "3-month Precip") %>%
	+ theme(legend.position = "none")%>%
	+ coord_cartesian(xlim=c(1, 365))

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_precip_jdate.png"), p, width =5, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "synth_precip_jdate.pdf"), p, width =5, height = 4.5)
ggsave(file.path(write_figures_path, "synth_precip_jdate.svg"), p, width =5, height = 4.5)


p <- ggplot(synth_stat_df %>% filter(zero == FALSE), aes(x=date, y=precip, colour=zero)) %>%
	+ geom_line() %>%
	+ geom_point(data = synth_stat_df %>% filter(zero == TRUE)) %>%
#	+ scale_x_continuous(name = "Julian Date", breaks=seq(1,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "3-month Precip") %>%
	+ scale_colour_manual(values = c("grey30", "#e41a1c")) %>%
	+ theme(legend.position = "none")

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_precip_fullseries.png"), p, width =7, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "synth_precip_fullseries.pdf"), p, width =7, height = 4.5)
ggsave(file.path(write_figures_path, "synth_precip_fullseries.svg"), p, width =7, height = 4.5)


p <- ggplot(synth_stat_df, aes(x=date, y=spi)) %>%
	+ geom_hline(yintercept = 0, size=0.15, colour="grey80") %>%
	+ geom_ribbon(aes(ymin = -99, ymax = spi_thresh), fill = "red", alpha=0.2) %>%
	+ geom_line(aes(colour=zero, group = year), size = 0.15) %>%
#	+ geom_point(data = y_df %>% filter(zero == TRUE)) %>%
	+ scale_x_date(name = "Date", breaks=seq(as.Date("1920-01-01"), as.Date("2020-01-01"), by = "10 years"), expand=c(0,0), date_labels = "%Y") %>%
	+ scale_y_continuous(name = "SPI") %>%
	+ scale_colour_manual(values = c("grey30", "#e41a1c")) %>%
	+ theme(legend.position = "none") %>%
	+ coord_cartesian(ylim=c(-3.5, 3.5))

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_spi_fullseries.png"), p, width =6, height = 3, dpi = 300)
ggsave(file.path(write_figures_path, "synth_spi_fullseries.pdf"), p, width =6, height = 3)
ggsave(file.path(write_figures_path, "synth_spi_fullseries.svg"), p, width =6, height = 3)











###########################################################################
###  Make up non-stationary time series
###########################################################################
### Make some daily synthetic distributions that vary smoothly through time
### Have to do a few things to make sure 1975 is identical for stationary and non-stationary
true_param_non_stat <- spi_true %>% 
	mutate( year_frac = decimal_date(date)) %>%
	mutate(shape = 1.75) %>%
	mutate(theta = case_when( jdate <=182 ~ 0.05 + (0.05/182)*jdate+0.015, 
		jdate > 182 ~ 0.1 - (0.05/183)*(jdate - 182)+.015)
	) %>%
	mutate(theta = case_when( year >= 1960 ~ theta - 0.001 *(year - 1960), 
		TRUE ~ theta
	)) %>%
	mutate(mean = sin(((2*pi)/365)*jdate)+0.4*sin(((10*pi)/365)*jdate)) %>%
	mutate(mean = 20 + mean + 0.05* mean*(year_frac-1900)) %>%  ### This line gives the mean a trend over long periods
	mutate(mean = mean*0.2094296 + 0.8132083) %>%
	mutate(scale = mean / shape) %>%
	mutate(rate = 1/scale) %>%
	mutate(disp = 1/shape) %>%
	select(-year_frac)


### Make sure that the stationary and nonstationary line up at year 1975
### Quick check plots
check_plot <- true_param_non_stat %>% left_join(true_param_stat, by ="jdate")

ggplot(check_plot, aes(x=date, y=mean.x)) + geom_line() + geom_line(aes(y=mean.y), colour="red")
ggplot(check_plot, aes(x=date, y=theta.x)) + geom_line() + geom_line(aes(y=theta.y), colour="red")

ggplot(check_plot %>% filter(year == 1975), aes(x=date, y=mean.x)) + geom_line() + geom_line(aes(y=mean.y), colour="red")

ggplot(check_plot %>% filter(year == 1975), aes(x=date, y=theta.x)) + geom_line() + geom_line(aes(y=theta.y), colour="red")

ggplot(check_plot, aes(x=jdate, y=mean.x, group = year)) + geom_line() + geom_line(aes(y=mean.y), colour="red")
ggplot(check_plot, aes(x=jdate, y=theta.x, group = year)) + geom_line() + geom_line(aes(y=theta.y), colour="red")


###########################################################################
###  Plot the underlying distribution
###########################################################################
write_output_path <- file.path(write_output_base, "nonstationary")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(write_figures_base, "nonstationary")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

### Plot how this looks	
p <- ggplot(true_param_non_stat, aes(x=jdate, y=mean, group=year)) %>%
	+  geom_line(aes(colour=year)) %>%
	+ scale_colour_viridis(name = "Year", option = "cividis") %>%
	+ scale_y_continuous(name = "True Mean") %>%
	+ scale_x_continuous(name = "Julian Date", breaks = round(seq(1, 365, length.out=13)), expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1, 365)) 
	
p	

### Plot surface for Mean
plot_df <- true_param_non_stat %>% filter(jdate <= 365)

p <- ggplot(plot_df, aes(x=jdate, y=year)) %>%
	+ geom_tile(aes(fill=mean)) %>%
	+ scale_fill_viridis(name = "True\nMean") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Year", expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1,365), ylim=c(1920,2020))

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_nonstat_mean_true_surface.png"), p, width =6.5, height = 6, dpi = 300)
ggsave(file.path(write_figures_path, "synth_nonstat_true_params_subset.pdf"), p, width =6.5, height = 6)
ggsave(file.path(write_figures_path, "synth_nonstat_true_params_subset.svg"), p, width =6.5, height = 6)

### Plot surface for Theta
plot_df <- true_param_non_stat %>% filter(jdate <= 365)

p <- ggplot(plot_df, aes(x=jdate, y=year)) %>%
	+ geom_tile(aes(fill=theta)) %>%
	+ scale_fill_viridis(name = "True\nTheta") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Year", expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1,365), ylim=c(1920,2020))

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_nonstat_theta_true_surface.png"), p, width =6.5, height = 6, dpi = 300)
ggsave(file.path(write_figures_path, "synth_nonstat_theta_params_subset.pdf"), p, width =6.5, height = 6)
ggsave(file.path(write_figures_path, "synth_nonstat_theta_params_subset.svg"), p, width =6.5, height = 6)


### Create a jdate plot
plot_df <- true_param_non_stat %>% 
	select(jdate, year, mean, shape, scale, theta) %>%
	gather("param", "value", -jdate, -year)

single_line_shape <- plot_df %>% filter(param == "shape" & year == 1920)

plot_df <- plot_df %>%
	filter(param != "shape") %>%
	bind_rows(single_line_shape) %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta")))
	
p <- ggplot(plot_df, aes(x=jdate, y=value, group=year)) %>% 
    + geom_line(aes(colour=year)) %>%
	+ scale_colour_viridis(name = "Year", option="inferno") %>%
	+ facet_grid(param ~ ., scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "True Parameter Value") %>%
	+ coord_cartesian(xlim=c(1, 365))

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_nonstat_true_params.png"), p, width =6, height = 7, dpi = 300)
ggsave(file.path(write_figures_path, "synth_nonstat_true_params.pdf"), p, width =6, height = 7)
ggsave(file.path(write_figures_path, "synth_nonstat_true_params.svg"), p, width =6, height = 7)



###########################################################################
###  Sample from this distribution
###########################################################################
### Apply day 365 to day 366, join with true parameters, and then return to the original format
synth_non_stat_df <- spi_true %>%
	mutate(jdate = case_when(
			jdate == 366 ~ 365,
			TRUE ~ jdate)
	) %>%
	left_join(true_param_non_stat) %>%
	mutate(jdate = yday(date)) 

### First transform SPI into percentile. Then set zeros when percentile is below theta threshold
### Stretch the remaining probability space so that it fills the space (0 to 1)
### Then sample from gamma distribution and combine with zeros
synth_non_stat_df <- synth_non_stat_df %>% 
	mutate(p_overall = pnorm(spi)) %>%
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
	) %>%
	select(-p_overall, precip_zero, p_precip, scale_factor) %>%
	mutate(zero = precip == 0)


###########################################################################
###  Plot the non-stationary syntehtic distribution
###########################################################################
p <- ggplot(synth_non_stat_df, aes(x=jdate, y=precip, colour=zero, group = year)) %>%
#	+ geom_jitter(alpha=0.2) %>%
	+ geom_line(alpha = 0.2) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(1,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "3-month Precip") %>%
	+ theme(legend.position = "none")%>%
	+ coord_cartesian(xlim=c(1, 365))

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_precip_jdate.png"), p, width =5, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "synth_precip_jdate.pdf"), p, width =5, height = 4.5)
ggsave(file.path(write_figures_path, "synth_precip_jdate.svg"), p, width =5, height = 4.5)


p <- ggplot(synth_non_stat_df %>% filter(zero == FALSE), aes(x=date, y=precip, colour=zero)) %>%
	+ geom_line() %>%
	+ geom_point(data = synth_non_stat_df %>% filter(zero == TRUE)) %>%
#	+ scale_x_continuous(name = "Julian Date", breaks=seq(1,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "3-month Precip") %>%
	+ scale_colour_manual(values = c("grey30", "#e41a1c")) %>%
	+ theme(legend.position = "none")

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_precip_fullseries.png"), p, width =7, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "synth_precip_fullseries.pdf"), p, width =7, height = 4.5)
ggsave(file.path(write_figures_path, "synth_precip_fullseries.svg"), p, width =7, height = 4.5)


p <- ggplot(synth_non_stat_df, aes(x=date, y=spi)) %>%
	+ geom_hline(yintercept = 0, size=0.15, colour="grey80") %>%
	+ geom_ribbon(aes(ymin = -99, ymax = spi_thresh), fill = "red", alpha=0.2) %>%
	+ geom_line(aes(colour=zero, group = year), size = 0.15) %>%
#	+ geom_point(data = y_df %>% filter(zero == TRUE)) %>%
	+ scale_x_date(name = "Date", breaks=seq(as.Date("1920-01-01"), as.Date("2020-01-01"), by = "10 years"), expand=c(0,0), date_labels = "%Y") %>%
	+ scale_y_continuous(name = "SPI") %>%
	+ scale_colour_manual(values = c("grey30", "#e41a1c")) %>%
	+ theme(legend.position = "none") %>%
	+ coord_cartesian(ylim=c(-3.5, 3.5))

p

### Save Figure
ggsave(file.path(write_figures_path, "synth_spi_fullseries.png"), p, width =6, height = 3, dpi = 300)
ggsave(file.path(write_figures_path, "synth_spi_fullseries.pdf"), p, width =6, height = 3)
ggsave(file.path(write_figures_path, "synth_spi_fullseries.svg"), p, width =6, height = 3)



###########################################################################
###  Save results from as an RDS file for analysis
###########################################################################

### Save results for next step
save(true_param_stat, true_param_non_stat, file = file.path(write_output_base, "true_param.RData"))
save(synth_stat_df, synth_non_stat_df, file = file.path(write_output_base, "synth_df.RData"))

write.csv(true_param_stat, file = file.path(write_output_base, "stationary/true_param_stat.csv"), row.names = FALSE)
write.csv(true_param_non_stat, file = file.path(write_output_base, "nonstationary/true_param_non_stat.csv"), row.names = FALSE)

write.csv(true_param_stat, file = file.path(write_output_base, "stationary/synth_stat_df.csv"), row.names = FALSE)
write.csv(true_param_non_stat, file = file.path(write_output_base, "nonstationary/synth_non_stat_df.csv"), row.names = FALSE)






