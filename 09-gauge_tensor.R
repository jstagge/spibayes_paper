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


###########################################################################
###  Create knots
###########################################################################
### Check the length
accum_df %>% group_by(short_name) %>% summarize(start = min(date), end = max(date))


### Create knots
n_knots_jdate <- 15

knots_jdate <- seq(1,365,length=n_knots_jdate)
knots_year <- seq(1855,2020,by = 15)

n_knots_year <- length(knots_year)

knot_loc <- list(jdate = knots_jdate, year = knots_year)
knot_loc


###########################################################################
###  Extract data to be fitted
###########################################################################
### Cut to only 365 days
fitting_df <- accum_df %>%
	filter(jdate <= 365) %>%
	drop_na(precip) %>%
	mutate(zero = precip == 0) %>%
	select(short_name, date, jdate, year, precip, zero)

### Check
head(fitting_df)
knot_loc

###########################################################################
###  Fit synthetic data using MLE
###########################################################################
### For supercomputer
#cl <- parallel::makeCluster(c("n1", "n2", "n3"))
#plan(cluster, workers = cl)
#parallel::stopCluster(cl)

### Set up the cores to run in parallel
availableCores()
#future::plan(multiprocess)
#future::plan(multiprocess(workers = 4))
future::plan(sequential)

### Function to fit each site
fit_each <- function(data, knot_loc, write_output_path, short_name){

	require(mgcv)
	require(spibayes)
	
	# Create an output folder
	write_output_path <- file.path(write_output_path, short_name)
	dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

	### Tweak the knots for simplicity
	year_extremes <- data %>%
		summarize(min = year(min(date)), max = year(max(date)))

	### Find the closest knot that is less than the start, go one knot before
	min_year_diff <- knot_loc$year-year_extremes$min
	min_year_diff[min_year_diff > 0] <- 999999
	min_knot <- max(c(which.min(abs(min_year_diff)) - 1, 1))

	### Find the closest knot that is more than the start, go one knot after
	max_year_diff <- knot_loc$year-year_extremes$max
	max_year_diff[max_year_diff < 0] <- 999999
	max_knot <- min(c(which.min(abs(max_year_diff)) + 1, length(knot_loc$year)))

	knot_loc$year <- knot_loc$year[seq(min_knot, max_knot)]
	save(knot_loc, file = file.path(write_output_path, paste0(short_name, "_knot_loc.rda")))
	
	### Preprocess
	tensor_init <- pre_proc(data = data, type = "tensor", knot_loc = knot_loc)

	### Run the nonBayes
	tensor_optimize <- model_fit(spi_input = tensor_init, iter = 2000, engine = "optimize", output_dir = write_output_path)
	#save(tensor_optimize, file = file.path(write_output_path, paste0(short_name, "_tensor_optimize.rda")))

	newdata_df <- expand.grid(jdate = seq(1,365,1), year = seq(year_extremes$min,year_extremes$max,1))
	param_est_optimize <- predict_vals(tensor_optimize, newdata = newdata_df)	
	output_file_loc <- tensor_optimize$model_fit$output_files()

	save(param_est_optimize, output_file_loc, file = file.path(write_output_path, paste0(short_name, "_tensor_optimize_paramest.rda")))

	### Run the full model
	tensor_sample <- model_fit(spi_input = tensor_init, n_chains = 1, cores = 1, iter = 800, engine = "sample", output_dir = write_output_path)
	#save(tensor_sample, file = file.path(write_output_path, paste0(short_name, "_tensor.rda")))

	param_est_sample <- predict_vals(tensor_sample, newdata = newdata_df)	
	output_file_loc <- tensor_sample$model_fit$output_files()
	save(param_est_sample, output_file_loc, file = file.path(write_output_path, paste0(short_name, "_tensor_sample_paramest.rda")))

	### Return the results
	return(list(file = output_file_loc, param_est_optimize = param_est_optimize, param_est_sample = param_est_sample))
}

### Split and run for each site, then combine
tensor_gauge_long <- fitting_df %>%
	#filter(short_name == "ABD" | short_name == "ALB") %>%
	split(.$short_name) %>% 
	future_map(~ fit_each(.x, knot_loc = knot_loc, write_output_path = write_output_path, short_name = .x$short_name[1]), .progress = FALSE) 

### Save everything
save(tensor_gauge_long, file = file.path(write_output_path, "tensor_gauge_long.rda"))



