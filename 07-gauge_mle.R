# *------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: .R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  Test the fundamental models, for a single distribution. Use a hurdle model to incorporate zeros
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
require(zoo)

### To run in parallel
library(furrr)

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
write_output_path <- file.path(output_path, "gauge_mle")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(output_path, "figures/gauge_mle")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)


###########################################################################
## Load Functions
###########################################################################

file.sources = list.files(file.path(here_path, "./functions"), pattern="*.R", recursive=TRUE)
sapply(file.path(file.path(here_path, "./functions"), file.sources),source)

###################################################################
###  Set initial values
###################################################################
### Calculate SPI 3, use 92 day rolling mean
n_roll <- 92

### How many NAs are acceptible in this rolling mean. Use 9 days or 10%
n_roll_min <- n_roll - 9

###########################################################################
## Load gauge precipitation time series
###########################################################################
load(file.path(output_path, "gauge_ghcnd/ghcnd_data.rda"))

ls()


###################################################################
###  Calculate accumualated precip
###################################################################
### Fill in all missing dates in the series and then create some columns for jdate, month_day, year
precip_df <- precip_df %>%
	group_by(short_name) %>%
	complete(date = seq.Date(min(date), max(date), by="1 day")) %>%
	arrange(date) %>%
	mutate(jdate = yday(date)) %>%
	mutate(year = year(date))

### Calculate the rolling mean
### First mutate calculates the rolling mean, second mutate calculates the number of non-na values, third mutate gives an NA where there are less observations than the threshold
accum_df <- precip_df %>%
	select(short_name, date, prcp_mm) %>%
	mutate(roll_mean_3 = rollmeanr(x=prcp_mm, k=n_roll, fill=NA, na.rm=TRUE)) %>%
	mutate(roll_mean_3_notna = rollsumr(x=!is.na(prcp_mm), k=n_roll, fill=NA, na.rm=TRUE)) %>%
	mutate(roll_mean_3 = case_when(roll_mean_3_notna  > n_roll_min ~ roll_mean_3,
		 TRUE ~ NA_real_)
	) %>%
	select(-roll_mean_3_notna, -prcp_mm)

head(accum_df)

### Reorganize
accum_df <- accum_df %>%
	arrange(short_name, date) %>%
	ungroup() %>%
	mutate(jdate = yday(date)) %>%
	mutate(year = year(date)) %>%
	rename(precip = roll_mean_3)

###################################################################
###  Test plot
###################################################################
p <- ggplot(accum_df, aes(x=date, y=precip, colour = short_name)) %>%
	+ geom_line() %>%
	+ scale_x_date(name = "Date") %>%
	+ scale_y_continuous(name = "3-Month Mean Precip (mm)") %>%
	+ facet_wrap(.~short_name) %>%
	+ scale_colour_brewer(name = "Gauge", type = "qual", palette = "Paired") %>%
	+ theme_bw(9)
p

p <- ggplot(accum_df, aes(x=jdate, y=precip, colour = short_name, group = year)) %>%
	+ geom_hline(yintercept = 0, colour = "grey20", size = 0.15) %>%
	+ geom_line(size=0.15) %>%
	+ scale_x_continuous(name = "Julian Day", breaks = seq(0, 365, 60)) %>%
	+ scale_y_continuous(name = "3-Month Mean Precip (mm)") %>%
	+ facet_wrap(.~short_name, scales = "free_y") %>%
	+ scale_colour_brewer(name = "Gauge", type = "qual", palette = "Paired") %>%
	+ theme_bw(9) %>%
	+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
### Save Figure
ggsave(file.path(write_figures_path, "precip_facet.png"), p, width =7, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "precip_facet.pdf"), p, width =7, height = 4.5)
ggsave(file.path(write_figures_path, "precip_facet.svg"), p, width =7, height = 4.5)

p <- ggplot(accum_df, aes(x=jdate, y=precip, colour = short_name, group = year)) %>%
	+ geom_hline(yintercept = 0, colour = "grey20", size = 0.15) %>%
	+ geom_line(size=0.15) %>%
	+ scale_x_continuous(name = "Julian Day", breaks = seq(0, 365, 60)) %>%
	+ scale_y_continuous(name = "3-Month Mean Precip (mm)") %>%
	+ facet_wrap(.~short_name) %>%
	+ scale_colour_brewer(name = "Gauge", type = "qual", palette = "Paired") %>%
	+ theme_bw(9)%>%
	+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

### Save Figure
ggsave(file.path(write_figures_path, "precip_facet_commonscale.png"), p, width =7, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "precip_facet_commonscale.pdf"), p, width =7, height = 4.5)
ggsave(file.path(write_figures_path, "precip_facet_commonscale.svg"), p, width =7, height = 4.5)


###########################################################################
###  Fit synthetic data using MLE
###########################################################################
### For supercomputer
#cl <- parallel::makeCluster(c("n1", "n2", "n3"))
#plan(cluster, workers = cl)
#parallel::stopCluster(cl)

### Set up the cores to run in parallel
availableCores()
future::plan(multiprocess)

### Function to fit each site
fit_each <- function(data){
	### Fit the gamma distribution using MLE	
	mle_output <- mle_gamma_fit(jdate = data$jdate, values = data$precip, n = 3000, method = "param", ci = 0.95)

	### Return the results
	return(mle_output)
}

### Split and run for each site, then combine
mle_gauge_long <- accum_df %>%
	split(.$short_name) %>% 
	future_map(~ fit_each(.x), .progress = TRUE) %>%
	transpose() %>%
	map(bind_rows, .id = "short_name") 

### Short uses the former WMO 30 year climatology
mle_gauge_short <- accum_df %>%
	filter(year >= 1961 & year <= 1990) %>%
	split(.$short_name) %>% 
	future_map(~ fit_each(.x), .progress = TRUE) %>%
	transpose() %>%
	map(bind_rows, .id = "short_name") 

###########################################################################
###  Plot the results
###########################################################################

### Plot the MLE result
mle_plot <- mle_gauge_long$estimate %>%
	pivot_wider(names_from = ci, values_from = value) %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta")))

p <- ggplot(mle_plot, aes(x=jdate)) %>%
	+ geom_ribbon(aes(ymin = lower_ci, ymax=upper_ci), fill="grey50", colour="grey50", alpha=0.5) %>%
	+ geom_line(colour="red", aes(y=estimate))  %>%
	+ facet_grid(param ~ short_name, scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Parameter Value") #%>%
#	+ coord_cartesian(xlim=c(1, 365))

p

### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_full.png"), p, width =10, height = 7, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_full.pdf"), p, width =10, height = 7)
ggsave(file.path(write_figures_path, "mle_fit_stat_full.svg"), p, width =10, height = 7)

p <- ggplot(mle_plot %>% filter(param != "Dispersion" & param != "Rate"), aes(x=jdate)) %>%
	+ geom_ribbon(aes(ymin = lower_ci, ymax=upper_ci), fill="grey50", colour="grey50", alpha=0.5) %>%
	+ geom_line(colour="red", aes(y=estimate))  %>%
	+ facet_grid(param ~ short_name, scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Parameter Value") %>%
	+ coord_cartesian(xlim=c(1, 365))

p
### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_subset.png"), p, width =10, height = 7, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_subset.pdf"), p, width =10, height = 7)
ggsave(file.path(write_figures_path, "mle_fit_stat_subset.svg"), p, width =10, height = 7)


### I could also remove the scale

###########################################################################
###  Save results from as an RDS file for analysis
###########################################################################

### Save results for next step
save(accum_df, file = file.path(output_path, "gauge_ghcnd/accum_df.rda"))
save(mle_gauge_long, mle_gauge_short, file = file.path(write_output_path, "mle_fit_gauge.rda"))






###########################################################################
###  Quick comparison of short and long
###########################################################################

### Plot the MLE result
mle_plot <- mle_gauge_long$estimate %>%
	mutate(length = "long") %>%
	bind_rows(mle_gauge_short$estimate %>% mutate(length = "short")) %>%
	pivot_wider(names_from = ci, values_from = value) %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta")))

p <- ggplot(mle_plot %>% filter(param != "Dispersion" & param != "Rate" & param != "Scale"), aes(x=jdate)) %>%
	+ geom_ribbon(aes(ymin = lower_ci, ymax=upper_ci, colour = length), fill = NA, alpha=0.5, size = 0.2) %>%
	+ geom_line(aes(y=estimate, colour = length))  %>%
	+ facet_grid(param ~ short_name, scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Parameter Value") %>%
	+ coord_cartesian(xlim=c(1, 365))

p


