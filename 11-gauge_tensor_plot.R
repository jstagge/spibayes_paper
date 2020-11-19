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

require(ncdf4)

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
write_output_path <- file.path(output_path, "gauge_tensor")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(output_path, "figures/gauge_tensor")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)


###########################################################################
## Load Functions
###########################################################################
file.sources = list.files(file.path(here_path, "./functions"), pattern="*.R", recursive=TRUE)
sapply(file.path(file.path(here_path, "./functions"), file.sources),source)




disp_to_shape <- function(disp){
	1/ exp(-7 + log(1 + exp(disp)))
}


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

###########################################################################
###  Loop through and save parameter estimates
###########################################################################
### Extract the site names
site_list <- levels(precip_df$short_name)
site_list <- site_list[site_list != "NYC"]

### Loop through each site
for (j in seq(1, length(site_list))){

short_name_j <- site_list[j]

cat(j)
cat(short_name_j)

### Quick plot of Initial Values
month_breaks <- c(yday(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month")), 365)
month_labels <- c(as.character(month(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month"), label=TRUE)), "Jan")

### Find all files within folder
#write_output_path <- "/run/media/jhstagge/Seagate Backup Plus Drive/spibayes_paper/output/gauge_tensor"
site_folder <- file.path(write_output_path, short_name_j)

### Read in Marginal
marginal_nc <- nc_open(file.path(site_folder, paste0(short_name_j, "_marginal.nc")))

date_vec <- ncvar_get(marginal_nc, "time") + as.Date("1900-01-01")
year_vec <- ncvar_get(marginal_nc, "year") 
jdate_vec <- ncvar_get(marginal_nc, "jdate") 

mean_int <- ncvar_get(marginal_nc, "mean", start = c(1,1,1), count = c(-1,-1,1)) 
mean_jdate <- ncvar_get(marginal_nc, "mean", start = c(1,1,2), count = c(-1,-1,1)) 
mean_year <- ncvar_get(marginal_nc, "mean", start = c(1,1,3), count = c(-1,-1,1)) 
mean_tensor <- ncvar_get(marginal_nc, "mean", start = c(1,1,4), count = c(-1,-1,1)) 
mean_estimate <- ncvar_get(marginal_nc, "mean", start = c(1,1,5), count = c(-1,-1,1)) 
mean_sigma<- ncvar_get(marginal_nc, "mean", start = c(1,1,6), count = c(-1,-1,1)) 

nc_close(marginal_nc)

### Do jdate
mean_marginal <- exp(mean_int + mean_jdate)
mean_lower <- exp(mean_int + mean_jdate + qnorm(0.025) * mean_sigma)
mean_upper <- exp(mean_int + mean_jdate + qnorm(0.975) * mean_sigma)

mean_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(mean_marginal, 1, mean))

mean_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(mean_marginal, 1, quantile, 0.025), 
	upper = apply(mean_marginal, 1, quantile, 0.975), 
	measure = "estimate")

mean_marginal_ribbon <- mean_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(mean_lower, 1, mean), 
	upper = apply(mean_upper, 1, mean), 
	measure = "sigma")
	)

plot_ribbon <- mean_marginal_ribbon %>%
	filter(year == 1990) %>%
	mutate(measure = factor(measure, levels = c("sigma", "estimate")))
plot_line <- mean_marginal_line %>%
	filter(year == 1990)

p <- ggplot(plot_line, aes(x=jdate)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = lower, ymax = upper, fill = measure), alpha=0.4) %>%
	+ geom_line(aes(y=estimate), colour = "red") %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>%
	#+ scale_colour_manual(name = "Shape", values = c("red", "black")) %>%
#	+ scale_fill_manual(name = NULL, values = c("grey70")) %>%
	+ scale_fill_manual(name = "95% CI", values = c("grey70", "grey30")) %>%
	+ scale_y_continuous(name="Mean Parameter") %>%
	+ coord_cartesian(xlim=c(1,365)) %>%
	+ theme(legend.position = "none") 
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_mean_jdate"))
ggsave(paste0(save_file,".png"), p, width =4.5, height = 3.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =4.5, height = 3.5)
ggsave(paste0(save_file,".svg"), p, width =4.5, height = 3.5)



p <- ggplot(plot_line, aes(x=jdate)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = lower, ymax = upper, fill = measure), alpha=0.4) %>%
	+ geom_line(aes(y=estimate), colour = "red") %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>%
	#+ scale_colour_manual(name = "Shape", values = c("red", "black")) %>%
#	+ scale_fill_manual(name = NULL, values = c("grey70")) %>%
	+ scale_fill_manual(name = "95% CI", values = c("grey70", "grey30")) %>%
	+ scale_y_continuous(name="Mean Parameter") %>%
	+ coord_cartesian(xlim=c(1,365), ylim=c(0,13)) %>%
	+ theme(legend.position = "none") 
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_mean_jdate_fixed"))
ggsave(paste0(save_file,".png"), p, width =4.5, height = 3.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =4.5, height = 3.5)
ggsave(paste0(save_file,".svg"), p, width =4.5, height = 3.5)


### Do year
mean_marginal <- exp(mean_int + mean_year)
mean_lower <- exp(mean_int + mean_year + qnorm(0.025) * mean_sigma)
mean_upper <- exp(mean_int + mean_year + qnorm(0.975) * mean_sigma)

mean_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(mean_marginal, 1, mean))

mean_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(mean_marginal, 1, quantile, 0.025), 
	upper = apply(mean_marginal, 1, quantile, 0.975), 
	measure = "estimate")

mean_marginal_ribbon <- mean_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(mean_lower, 1, mean), 
	upper = apply(mean_upper, 1, mean), 
	measure = "sigma")
	)

plot_ribbon <- mean_marginal_ribbon %>%
	filter(jdate == 1) %>%
	mutate(measure = factor(measure, levels = c("sigma", "estimate")))
plot_line <- mean_marginal_line %>%
	filter(jdate == 1)

p <- ggplot(plot_line, aes(x=year)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = lower, ymax = upper, fill = measure), alpha=0.4) %>%
	+ geom_line(aes(y=estimate), colour = "red") %>%
	+ scale_x_continuous(name = "Year", breaks=seq(1800,2020,10), expand = c(0,0)) %>%
	+ scale_fill_manual(name = "95% CI", values = c("grey70", "grey30")) %>%
	+ scale_y_continuous(name="Mean Parameter") %>%
	+ coord_cartesian(xlim=c(1880,2020)) %>%
	+ theme(legend.position = "none") 
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_mean_year"))
ggsave(paste0(save_file,".png"), p, width =4.5, height = 3.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =4.5, height = 3.5)
ggsave(paste0(save_file,".svg"), p, width =4.5, height = 3.5)



p <- ggplot(plot_line, aes(x=year)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = lower, ymax = upper, fill = measure), alpha=0.4) %>%
	+ geom_line(aes(y=estimate), colour = "red") %>%
	+ scale_x_continuous(name = "Year", breaks=seq(1800,2020,10), expand = c(0,0)) %>%
	+ scale_fill_manual(name = "95% CI", values = c("grey70", "grey30")) %>%
	+ scale_y_continuous(name="Mean Parameter") %>%
	+ coord_cartesian(xlim=c(1880,2020), ylim=c(0,13)) %>%
	+ theme(legend.position = "none") 
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_mean_year_fixed"))
ggsave(paste0(save_file,".png"), p, width =4.5, height = 3.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =4.5, height = 3.5)
ggsave(paste0(save_file,".svg"), p, width =4.5, height = 3.5)


### Do tensor
mean_marginal <- exp(mean_int + mean_tensor)
mean_lower <- exp(mean_int + mean_tensor + qnorm(0.025) * mean_sigma)
mean_upper <- exp(mean_int + mean_tensor + qnorm(0.975) * mean_sigma)

mean_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(mean_marginal, 1, mean))

mean_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(mean_marginal, 1, quantile, 0.025), 
	upper = apply(mean_marginal, 1, quantile, 0.975), 
	measure = "estimate")

mean_marginal_ribbon <- mean_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(mean_lower, 1, mean), 
	upper = apply(mean_upper, 1, mean), 
	measure = "sigma")
	)

plot_line <- mean_marginal_line %>%
	filter(jdate <= 365)

p <- ggplot(plot_line, aes(x=jdate, y=year)) %>%
	+ geom_tile(aes(fill = estimate)) %>%
	+ scale_fill_viridis(name = "Mean") %>%
	+ theme_bw(8) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_y_continuous(name="Year", breaks = seq(1800,2020,20),  expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1,365), ylim = c(1880, 2020))
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_mean_tensor"))
ggsave(paste0(save_file,".png"), p, width =5.5, height = 5.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =5.5, height = 5.5)
ggsave(paste0(save_file,".svg"), p, width =5.5, height = 5.5)


### Do full
mean_marginal <- exp(mean_int + mean_jdate + mean_year + mean_tensor)
mean_lower <- exp(mean_int + mean_jdate + mean_year + mean_tensor + qnorm(0.025) * mean_sigma)
mean_upper <- exp(mean_int + mean_jdate + mean_year + mean_tensor + qnorm(0.975) * mean_sigma)

mean_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(mean_marginal, 1, mean))

mean_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(mean_marginal, 1, quantile, 0.025), 
	upper = apply(mean_marginal, 1, quantile, 0.975), 
	measure = "estimate")

mean_marginal_ribbon <- mean_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(mean_lower, 1, mean), 
	upper = apply(mean_upper, 1, mean), 
	measure = "sigma")
	)

plot_line <- mean_marginal_line %>%
	filter(jdate <= 365)

p <- ggplot(plot_line, aes(x=jdate, y=year)) %>%
	+ geom_tile(aes(fill = estimate)) %>%
	+ scale_fill_viridis(name = "Mean") %>%
	+ theme_bw(8) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_y_continuous(name="Year", breaks = seq(1800,2020,20),  expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1,365), ylim = c(1880, 2020))
### Plot
p


### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_mean_full"))
ggsave(paste0(save_file,".png"), p, width =5.5, height = 5.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =5.5, height = 5.5)
ggsave(paste0(save_file,".svg"), p, width =5.5, height = 5.5)

rm(mean_marginal_line)
rm(mean_marginal_ribbon)
rm(plot_line)

rm(mean_int)
rm(mean_jdate)
rm(mean_year)
rm(mean_tensor) 
rm(mean_estimate)
rm(mean_sigma)
}



###########################################################################
###  Loop through and save parameter estimates
###########################################################################
### Extract the site names
site_list <- levels(precip_df$short_name)
site_list <- site_list[site_list != "NYC"]

### Loop through each site
for (j in seq(1, length(site_list))){

short_name_j <- site_list[j]

cat(j)
cat(short_name_j)

### Quick plot of Initial Values
month_breaks <- c(yday(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month")), 365)
month_labels <- c(as.character(month(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month"), label=TRUE)), "Jan")

### Find all files within folder
#write_output_path <- "/run/media/jhstagge/Seagate Backup Plus Drive/spibayes_paper/output/gauge_tensor"
site_folder <- file.path(write_output_path, short_name_j)

### Read in Marginal
marginal_nc <- nc_open(file.path(site_folder, paste0(short_name_j, "_marginal.nc")))

date_vec <- ncvar_get(marginal_nc, "time") + as.Date("1900-01-01")
year_vec <- ncvar_get(marginal_nc, "year") 
jdate_vec <- ncvar_get(marginal_nc, "jdate") 

disp_int <- ncvar_get(marginal_nc, "disp", start = c(1,1,1), count = c(-1,-1,1)) 
disp_jdate <- ncvar_get(marginal_nc, "disp", start = c(1,1,2), count = c(-1,-1,1)) 
disp_year <- ncvar_get(marginal_nc, "disp", start = c(1,1,3), count = c(-1,-1,1)) 
disp_tensor <- ncvar_get(marginal_nc, "disp", start = c(1,1,4), count = c(-1,-1,1)) 
disp_estimate <- ncvar_get(marginal_nc, "disp", start = c(1,1,5), count = c(-1,-1,1)) 
disp_sigma<- ncvar_get(marginal_nc, "disp", start = c(1,1,6), count = c(-1,-1,1)) 

nc_close(marginal_nc)

### Do jdate
shape_marginal <- disp_to_shape(disp_int + disp_jdate)
shape_lower <- disp_to_shape(disp_int + disp_jdate + qnorm(0.025) * disp_sigma)
shape_upper <- disp_to_shape(disp_int + disp_jdate + qnorm(0.975) * disp_sigma)

shape_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(shape_marginal, 1, mean))

shape_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(shape_marginal, 1, quantile, 0.025), 
	upper = apply(shape_marginal, 1, quantile, 0.975), 
	measure = "estimate")

shape_marginal_ribbon <- shape_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(shape_lower, 1, mean), 
	upper = apply(shape_upper, 1, mean), 
	measure = "sigma")
	)

plot_ribbon <- shape_marginal_ribbon %>%
	filter(year == 1990) %>%
	mutate(measure = factor(measure, levels = c("sigma", "estimate")))
plot_line <- shape_marginal_line %>%
	filter(year == 1990)

p <- ggplot(plot_line, aes(x=jdate)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = lower, ymax = upper, fill = measure), alpha=0.4) %>%
	+ geom_line(aes(y=estimate), colour = "red") %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>%
	+ scale_fill_manual(name = "95% CI", values = c("grey70", "grey30")) %>%
	+ scale_y_continuous(name="Shape Parameter") %>%
	+ coord_cartesian(xlim=c(1,365)) %>%
	+ theme(legend.position = "none") 
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_shape_jdate"))
ggsave(paste0(save_file,".png"), p, width =4.5, height = 3.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =4.5, height = 3.5)
ggsave(paste0(save_file,".svg"), p, width =4.5, height = 3.5)



### Do year
shape_marginal <- disp_to_shape(disp_int + disp_year)
shape_lower <- disp_to_shape(disp_int + disp_year + qnorm(0.025) * disp_sigma)
shape_upper <- disp_to_shape(disp_int + disp_year + qnorm(0.975) * disp_sigma)


shape_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(shape_marginal, 1, mean))

shape_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(shape_marginal, 1, quantile, 0.025), 
	upper = apply(shape_marginal, 1, quantile, 0.975), 
	measure = "estimate")

shape_marginal_ribbon <- shape_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(shape_lower, 1, mean), 
	upper = apply(shape_upper, 1, mean), 
	measure = "sigma")
	)

plot_ribbon <- shape_marginal_ribbon %>%
	filter(jdate == 1) %>%
	mutate(measure = factor(measure, levels = c("sigma", "estimate")))
plot_line <- shape_marginal_line %>%
	filter(jdate == 1)

p <- ggplot(plot_line, aes(x=year)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = lower, ymax = upper, fill = measure), alpha=0.4) %>%
	+ geom_line(aes(y=estimate), colour = "red") %>%
	+ scale_x_continuous(name = "Year", breaks=seq(1800,2020,10), expand = c(0,0)) %>%
	+ scale_fill_manual(name = "95% CI", values = c("grey70", "grey30")) %>%
	+ scale_y_continuous(name="Shape Parameter") %>%
	+ coord_cartesian(xlim=c(1880,2020)) %>%
	+ theme(legend.position = "none") 
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_shape_year"))
ggsave(paste0(save_file,".png"), p, width =4.5, height = 3.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =4.5, height = 3.5)
ggsave(paste0(save_file,".svg"), p, width =4.5, height = 3.5)


### Do tensor
shape_marginal <- disp_to_shape(disp_int + disp_tensor)
shape_lower <- disp_to_shape(disp_int + disp_tensor + qnorm(0.025) * disp_sigma)
shape_upper <- disp_to_shape(disp_int + disp_tensor + qnorm(0.975) * disp_sigma)

shape_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(shape_marginal, 1, mean))

shape_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(shape_marginal, 1, quantile, 0.025), 
	upper = apply(shape_marginal, 1, quantile, 0.975), 
	measure = "estimate")

shape_marginal_ribbon <- shape_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(shape_lower, 1, mean), 
	upper = apply(shape_upper, 1, mean), 
	measure = "sigma")
	)

plot_line <- shape_marginal_line %>%
	filter(jdate <= 365)

p <- ggplot(plot_line, aes(x=jdate, y=year)) %>%
	+ geom_tile(aes(fill = estimate)) %>%
	+ scale_fill_viridis(name = "Shape") %>%
	+ theme_bw(8) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_y_continuous(name="Year", breaks = seq(1800,2020,20),  expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1,365), ylim = c(1880, 2020))
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_shape_tensor"))
ggsave(paste0(save_file,".png"), p, width =5.5, height = 5.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =5.5, height = 5.5)
ggsave(paste0(save_file,".svg"), p, width =5.5, height = 5.5)


### Do full
shape_marginal <- disp_to_shape(disp_int + disp_jdate + disp_year + disp_tensor)
shape_lower <- disp_to_shape(disp_int + disp_jdate + disp_year + disp_tensor + qnorm(0.025) * disp_sigma)
shape_upper <- disp_to_shape(disp_int + disp_jdate + disp_year + disp_tensor + qnorm(0.975) * disp_sigma)

shape_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(shape_marginal, 1, mean))

shape_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(shape_marginal, 1, quantile, 0.025), 
	upper = apply(shape_marginal, 1, quantile, 0.975), 
	measure = "estimate")

shape_marginal_ribbon <- shape_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(shape_lower, 1, mean), 
	upper = apply(shape_upper, 1, mean), 
	measure = "sigma")
	)

plot_line <- shape_marginal_line %>%
	filter(jdate <= 365)

p <- ggplot(plot_line, aes(x=jdate, y=year)) %>%
	+ geom_tile(aes(fill = estimate)) %>%
	+ scale_fill_viridis(name = "Shape") %>%
	+ theme_bw(8) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_y_continuous(name="Year", breaks = seq(1800,2020,20),  expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1,365), ylim = c(1880, 2020))
### Plot
p


### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_shape_full"))
ggsave(paste0(save_file,".png"), p, width =5.5, height = 5.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =5.5, height = 5.5)
ggsave(paste0(save_file,".svg"), p, width =5.5, height = 5.5)


rm(shape_marginal_line)
rm(shape_marginal_ribbon)
rm(plot_line)

rm(disp_int)
rm(disp_jdate)
rm(disp_year)
rm(disp_tensor) 
rm(disp_estimate)
rm(disp_sigma)
}





###########################################################################
###  Loop through and save parameter estimates
###########################################################################
### Extract the site names
site_list <- levels(precip_df$short_name)
site_list <- site_list[site_list != "NYC"]

### Loop through each site
for (j in seq(1, length(site_list))){

short_name_j <- site_list[j]

cat(j)
cat(short_name_j)

### Quick plot of Initial Values
month_breaks <- c(yday(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month")), 365)
month_labels <- c(as.character(month(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month"), label=TRUE)), "Jan")

### Find all files within folder
#write_output_path <- "/run/media/jhstagge/Seagate Backup Plus Drive/spibayes_paper/output/gauge_tensor"
site_folder <- file.path(write_output_path, short_name_j)

### Read in Marginal
marginal_nc <- nc_open(file.path(site_folder, paste0(short_name_j, "_marginal.nc")))

date_vec <- ncvar_get(marginal_nc, "time") + as.Date("1900-01-01")
year_vec <- ncvar_get(marginal_nc, "year") 
jdate_vec <- ncvar_get(marginal_nc, "jdate") 

theta_int <- ncvar_get(marginal_nc, "theta", start = c(1,1,1), count = c(-1,-1,1)) 
theta_jdate <- ncvar_get(marginal_nc, "theta", start = c(1,1,2), count = c(-1,-1,1)) 
theta_year <- ncvar_get(marginal_nc, "theta", start = c(1,1,3), count = c(-1,-1,1)) 
theta_tensor <- ncvar_get(marginal_nc, "theta", start = c(1,1,4), count = c(-1,-1,1)) 
theta_estimate <- ncvar_get(marginal_nc, "theta", start = c(1,1,5), count = c(-1,-1,1)) 
theta_sigma<- ncvar_get(marginal_nc, "theta", start = c(1,1,6), count = c(-1,-1,1)) 

nc_close(marginal_nc)

### Do jdate
theta_marginal <- exp(theta_int + theta_jdate)/(1 + exp(theta_int + theta_jdate))
theta_lower <- exp(theta_int + theta_jdate + qnorm(0.025) * theta_sigma)/(1 + exp(theta_int + theta_jdate+ qnorm(0.025) * theta_sigma))
theta_upper <-  exp(theta_int + theta_jdate + qnorm(0.975) * theta_sigma)/(1 + exp(theta_int + theta_jdate+ qnorm(0.975) * theta_sigma))

theta_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(theta_marginal, 1, mean))

theta_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(theta_marginal, 1, quantile, 0.025), 
	upper = apply(theta_marginal, 1, quantile, 0.975), 
	measure = "estimate")

theta_marginal_ribbon <- theta_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(theta_lower, 1, mean), 
	upper = apply(theta_upper, 1, mean), 
	measure = "sigma")
	)

plot_ribbon <- theta_marginal_ribbon %>%
	filter(year == 1990) %>%
	mutate(measure = factor(measure, levels = c("sigma", "estimate")))
plot_line <- theta_marginal_line %>%
	filter(year == 1990)

p <- ggplot(plot_line, aes(x=jdate)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = lower, ymax = upper, fill = measure), alpha=0.4) %>%
	+ geom_line(aes(y=estimate), colour = "red") %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>%
	+ scale_fill_manual(name = "95% CI", values = c("grey70", "grey30")) %>%
	+ scale_y_continuous(name="Theta Parameter", labels = scales::percent) %>%
	+ coord_cartesian(xlim=c(1,365)) %>%
	+ theme(legend.position = "none") 
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_theta_jdate"))
ggsave(paste0(save_file,".png"), p, width =4.5, height = 3.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =4.5, height = 3.5)
ggsave(paste0(save_file,".svg"), p, width =4.5, height = 3.5)

p <- p 	+ coord_cartesian(xlim=c(1,365), ylim=c(0,0.6))  
save_file <- file.path(write_figures_path, paste0(short_name_j, "_theta_jdate_fixed"))
ggsave(paste0(save_file,".png"), p, width =4.5, height = 3.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =4.5, height = 3.5)
ggsave(paste0(save_file,".svg"), p, width =4.5, height = 3.5)


### Do year
theta_marginal <- exp(theta_int + theta_year)/(1 + exp(theta_int + theta_year))
theta_lower <- exp(theta_int + theta_year + qnorm(0.025) * theta_sigma)/(1 + exp(theta_int + theta_year+ qnorm(0.025) * theta_sigma))
theta_upper <-  exp(theta_int + theta_year + qnorm(0.975) * theta_sigma)/(1 + exp(theta_int + theta_year+ qnorm(0.975) * theta_sigma))

theta_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(theta_marginal, 1, mean))

theta_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(theta_marginal, 1, quantile, 0.025), 
	upper = apply(theta_marginal, 1, quantile, 0.975), 
	measure = "estimate")

theta_marginal_ribbon <- theta_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(theta_lower, 1, mean), 
	upper = apply(theta_upper, 1, mean), 
	measure = "sigma")
	)

plot_ribbon <- theta_marginal_ribbon %>%
	filter(jdate == 1) %>%
	mutate(measure = factor(measure, levels = c("sigma", "estimate")))
plot_line <- theta_marginal_line %>%
	filter(jdate == 1)

p <- ggplot(plot_line, aes(x=year)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = lower, ymax = upper, fill = measure), alpha=0.4) %>%
	+ geom_line(aes(y=estimate), colour = "red") %>%
	+ scale_x_continuous(name = "Year", breaks=seq(1800,2020,10), expand = c(0,0)) %>%
	+ scale_fill_manual(name = "95% CI", values = c("grey70", "grey30")) %>%
	+ scale_y_continuous(name="Theta Parameter", labels = scales::percent) %>%
	+ coord_cartesian(xlim=c(1880,2020)) %>%
	+ theme(legend.position = "none") 
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_theta_year"))
ggsave(paste0(save_file,".png"), p, width =4.5, height = 3.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =4.5, height = 3.5)
ggsave(paste0(save_file,".svg"), p, width =4.5, height = 3.5)

p <- p 	+ coord_cartesian(xlim=c(1,365), ylim=c(0,0.6))  
save_file <- file.path(write_figures_path, paste0(short_name_j, "_theta_year_fixed"))
ggsave(paste0(save_file,".png"), p, width =4.5, height = 3.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =4.5, height = 3.5)
ggsave(paste0(save_file,".svg"), p, width =4.5, height = 3.5)

### Do tensor
theta_marginal <- exp(theta_int + theta_tensor)/(1 + exp(theta_int + theta_tensor))
theta_lower <- exp(theta_int + theta_tensor + qnorm(0.025) * theta_sigma)/(1 + exp(theta_int + theta_tensor+ qnorm(0.025) * theta_sigma))
theta_upper <-  exp(theta_int + theta_tensor + qnorm(0.975) * theta_sigma)/(1 + exp(theta_int + theta_tensor+ qnorm(0.975) * theta_sigma))

theta_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(theta_marginal, 1, mean))

theta_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(theta_marginal, 1, quantile, 0.025), 
	upper = apply(theta_marginal, 1, quantile, 0.975), 
	measure = "estimate")

theta_marginal_ribbon <- theta_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(theta_lower, 1, mean), 
	upper = apply(theta_upper, 1, mean), 
	measure = "sigma")
	)


plot_line <- theta_marginal_line %>%
	filter(jdate <= 365)

p <- ggplot(plot_line, aes(x=jdate, y=year)) %>%
	+ geom_tile(aes(fill = estimate)) %>%
	+ scale_fill_viridis(name="Theta Parameter", labels = scales::percent) %>%
	+ theme_bw(8) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_y_continuous(name="Year", breaks = seq(1800,2020,20),  expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1,365), ylim = c(1880, 2020))
### Plot
p

### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_theta_tensor"))
ggsave(paste0(save_file,".png"), p, width =5.5, height = 5.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =5.5, height = 5.5)
ggsave(paste0(save_file,".svg"), p, width =5.5, height = 5.5)


### Do full
theta_marginal <- exp(theta_int + theta_jdate + theta_year + theta_tensor)/(1 + exp(theta_int + theta_jdate + theta_year + theta_tensor))
theta_lower <- exp(theta_int + theta_jdate + theta_year + theta_tensor + qnorm(0.025) * theta_sigma)/(1 + exp(theta_int + theta_jdate + theta_year + theta_tensor+ qnorm(0.025) * theta_sigma))
theta_upper <-  exp(theta_int + theta_jdate + theta_year + theta_tensor + qnorm(0.975) * theta_sigma)/(1 + exp(theta_int + theta_jdate + theta_year + theta_tensor + qnorm(0.975) * theta_sigma))

theta_marginal_line <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, estimate = apply(theta_marginal, 1, mean))

theta_marginal_ribbon <- data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(theta_marginal, 1, quantile, 0.025), 
	upper = apply(theta_marginal, 1, quantile, 0.975), 
	measure = "estimate")

theta_marginal_ribbon <- theta_marginal_ribbon %>%
	rbind(data.frame(date = date_vec, 
	year = year_vec, 
	jdate = jdate_vec, 
	lower = apply(theta_lower, 1, mean), 
	upper = apply(theta_upper, 1, mean), 
	measure = "sigma")
	)


plot_line <- theta_marginal_line %>%
	filter(jdate <= 365)

p <- ggplot(plot_line, aes(x=jdate, y=year)) %>%
	+ geom_tile(aes(fill = estimate)) %>%
	+ scale_fill_viridis(name="Theta Parameter", labels = scales::percent) %>%
	+ theme_bw(8) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_y_continuous(name="Year", breaks = seq(1800,2020,20),  expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1,365), ylim = c(1880, 2020))
### Plot
p


### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_theta_full"))
ggsave(paste0(save_file,".png"), p, width =5.5, height = 5.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =5.5, height = 5.5)
ggsave(paste0(save_file,".svg"), p, width =5.5, height = 5.5)



rm(theta_marginal_line)
rm(theta_marginal_ribbon)
rm(plot_line)

rm(theta_int)
rm(theta_jdate)
rm(theta_year)
rm(theta_tensor) 
rm(theta_estimate)
rm(theta_sigma)
}









j <- 1

short_name_j <- site_list[j]

cat(j)
cat(short_name_j)

### Quick plot of Initial Values
month_breaks <- c(yday(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month")), 365)
month_labels <- c(as.character(month(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month"), label=TRUE)), "Jan")

### Find all files within folder
#write_output_path <- "/run/media/jhstagge/Seagate Backup Plus Drive/spibayes_paper/output/gauge_tensor"
site_folder <- file.path(write_output_path, short_name_j)

estimate_nc <- nc_open(file.path(site_folder, paste0(short_name_j, "_estimate.nc")))

date_vec <- ncvar_get(estimate_nc, "time") + as.Date("1900-01-01")
year_vec <- ncvar_get(estimate_nc, "year") 
jdate_vec <- ncvar_get(estimate_nc, "jdate") 

mean_mat <- ncvar_get(estimate_nc, "mean") 
shape_mat <- ncvar_get(estimate_nc, "shape") 
scale_mat <- ncvar_get(estimate_nc, "scale") 

dist_df <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, shape = apply(shape_mat, 1, mean), scale = apply(scale_mat, 1, mean))

x_seq <- seq(0.1,20, 0.1)
x_date <- seq(as.Date("1950-01-01"), as.Date("1950-12-1"), by = "1 month")
first_ofmonth <- yday(x_date)
year_seq <-  seq(1810, 2020, 15)

plot_dist <- tibble(year = NA, jdate = NA, date = NA, precip = NA, density = NA, cumul = NA)

for (k in seq(1, length(year_seq))) {
	for(i in seq(1,length(first_ofmonth))){

	dist_temp <- dist_df %>%
		filter(year == year_seq[k] & jdate %in%  first_ofmonth[i])
	
	if(dim(dist_temp)[1] == 0) {next}
	
	dist_temp_full <- tibble(year = year_seq[k], 
		jdate = first_ofmonth[i], 
		date = x_date[i],
		precip= x_seq, 
		density = dgamma(x_seq, shape = dist_temp$shape[1], scale = dist_temp$scale[1]),
		cumul = pgamma(x_seq, shape = dist_temp$shape[1], scale = dist_temp$scale[1])
	)
	
	plot_dist <- plot_dist %>%
		bind_rows(dist_temp_full)

	rm(dist_temp_full)
	rm(dist_temp)
	}
}

plot_dist <- plot_dist %>%
	drop_na()


ggplot(plot_dist %>% filter(year != 1900 & year != 2020), aes(x=precip, y=density)) + geom_density(stat = "identity", fill = "grey50") + facet_grid(year ~ date) + coord_cartesian(xlim=c(0,6))


ggplot(plot_dist, aes(x=precip, y=density)) + geom_density(stat = "identity", fill = "grey50") + facet_grid(year ~ date)


ggplot(plot_dist, aes(x=precip, y=density, colour = year, group = year)) + geom_line() + facet_grid(. ~ jdate) + scale_colour_viridis(option = "inferno")





### Loop through each site
for (j in seq(1, length(site_list))){

short_name_j <- site_list[j]

cat(j)
cat(short_name_j)

### Quick plot of Initial Values
month_breaks <- c(yday(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month")), 365)
month_labels <- c(as.character(month(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month"), label=TRUE)), "Jan")

### Find all files within folder
#write_output_path <- "/run/media/jhstagge/Seagate Backup Plus Drive/spibayes_paper/output/gauge_tensor"
site_folder <- file.path(write_output_path, short_name_j)


### Quick plot of Initial Values
month_breaks <- c(yday(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month")), 365)
month_labels <- c(as.character(month(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month"), label=TRUE)), "Jan")

### Find all files within folder
#write_output_path <- "/run/media/jhstagge/Seagate Backup Plus Drive/spibayes_paper/output/gauge_tensor"
site_folder <- file.path(write_output_path, short_name_j)

estimate_nc <- nc_open(file.path(site_folder, paste0(short_name_j, "_estimate.nc")))

date_vec <- ncvar_get(estimate_nc, "time") + as.Date("1900-01-01")
year_vec <- ncvar_get(estimate_nc, "year") 
jdate_vec <- ncvar_get(estimate_nc, "jdate") 

mean_mat <- ncvar_get(estimate_nc, "mean") 
shape_mat <- ncvar_get(estimate_nc, "shape") 
scale_mat <- ncvar_get(estimate_nc, "scale") 

dist_df <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, shape = apply(shape_mat, 1, mean), scale = apply(scale_mat, 1, mean))

x_seq <- seq(0.1,20, 0.1)
x_date <- seq(as.Date("1950-01-01"), as.Date("1950-12-1"), by = "1 month")
first_ofmonth <- yday(x_date)
year_seq <-  seq(1810, 2020, 15)


plot_dist <- tibble(year = NA, jdate = NA, date = NA, ymin = NA, lower = NA, middle = NA, upper = NA, ymax = NA)

for (k in seq(1, length(year_seq))) {
	for(i in seq(1,length(first_ofmonth))){

	dist_temp <- dist_df %>%
		filter(year == year_seq[k] & jdate %in%  first_ofmonth[i])
	
	if(dim(dist_temp)[1] == 0) {next}
	
	#dist_temp_full <- tibble(year = year_seq[k], 
	#	jdate = first_ofmonth[i], 
	#	precip = rgamma(50000, shape = dist_temp$shape[1], scale = dist_temp$scale[1])
	#)

	dist_temp_full <- tibble(year = year_seq[k], 
		jdate = first_ofmonth[i], 
		date = x_date[i],
		ymin = qgamma(0.025, shape = dist_temp$shape[1], scale = dist_temp$scale[1]),
		lower =qgamma(0.25, shape = dist_temp$shape[1], scale = dist_temp$scale[1]),
		middle = qgamma(0.5, shape = dist_temp$shape[1], scale = dist_temp$scale[1]),
		upper = qgamma(0.75, shape = dist_temp$shape[1], scale = dist_temp$scale[1]),
		ymax = qgamma(0.975, shape = dist_temp$shape[1], scale = dist_temp$scale[1])
	)
	
	plot_dist <- plot_dist %>%
		bind_rows(dist_temp_full)

	rm(dist_temp_full)
	rm(dist_temp)
	}
}

plot_dist <- plot_dist %>%
	drop_na() %>%
	mutate(plot_month = paste0(month(date, label=TRUE), "-", day(date))) %>%
	mutate(plot_month = factor(plot_month, levels = unique(paste0(month(plot_dist$date, label=TRUE), "-", day(plot_dist$date)))))


p <- ggplot(plot_dist , aes(x=plot_month, lower=lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax, fill = factor(year))) %>%
	+ geom_boxplot(stat = "identity", colour = "grey30", alpha = 0.9, size = 0.35) %>%
	+ scale_fill_viridis(name = "Year", option = "inferno", discrete = TRUE) %>%
	+ scale_y_continuous(name = "3-month precip (in)") %>%
	+ scale_x_discrete(name = "Month") %>%
	+ theme_classic(11)


### Save Figure
save_file <- file.path(write_figures_path, paste0(short_name_j, "_boxplot"))
ggsave(paste0(save_file,".png"), p, width =12.5, height = 4.5, dpi = 300)
ggsave(paste0(save_file,".pdf"), p, width =12.5, height = 4.5)
ggsave(paste0(save_file,".svg"), p, width =12.5, height = 4.5)


}






	%>%
	+ scale_x_date(NULL, date_labels = "%b %y", breaks = "month")



	%>%
	+ 
















plot(date_vec, mean_vec, type="l")

nc_close(estimate_nc)

plot_df <- data.frame(date = date_vec, year = year_vec, jdate = jdate_vec, mean = mean_vec)











### Extract the site names
#site_list <- levels(precip_df$short_name)


#for (j in seq(1, length(site_list))){

#short_name_j <- site_list[j]

#cat(j)
#cat(short_name_j)

### Find all files within folder
#site_folder <- file.path(write_output_path, short_name_j)

### Find the biggest. Could do a newest
#site_files <- Sys.glob(file.path(site_folder, "*.csv"))
#site_files_info <- file.info(site_files)
#file_j <- rownames(site_files_info)[which.max(site_files_info$size)]
#order(site_files_info$mtime)

#load(file.path(site_folder, paste0(short_name_j, "_knot_loc.rda")))
#load(file.path(site_folder, paste0(short_name_j, "_tensor.rda")))

#precip_j <- precip_df %>% filter(short_name == short_name_j) %>% mutate(year = year(date))

### Estimate parameter values from model using initial beta estimates
#newdata_df <- expand.grid(jdate = seq(1,365,1), year = seq(min(precip_j$year),max(precip_j$year),1))
#save_file <- file_j
#file.path(site_folder, file_j)

#param_est <- predict_vals(tensor_sample, newdata = newdata_df, saved_model = save_file)

### Save everything
#save(param_est, file = file.path(site_folder, paste0(short_name_j, "_paramest.rda")))

#rm(param_est)
#rm(newdata_df)

#}


