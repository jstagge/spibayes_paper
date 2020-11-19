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
write_output_path <- file.path(output_path, "gauge_cyclic")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(output_path, "figures/gauge_gauge_cyclic")
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
### Create knots
n_knots_jdate <- 15

knot_loc <- list(jdate = seq(1,365,length=n_knots_jdate))

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
future::plan(multiprocess)

### Function to fit each site
fit_each <- function(data, knot_loc, write_output_path, short_name){

	require(mgcv)
	require(spibayes)
	
	#print("Hostname: ", Sys.info(), "\n")
	#print(short_name)
	write_output_path <- file.path(write_output_path, short_name)
	dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

	### Preprocess
	cyclic_init <- pre_proc(data = data, type = "cyclic", knot_loc = knot_loc)


	write.csv(x=data.frame(progress = c("starting_optimize")), file=file.path(write_output_path, paste0(short_name, "_starting_optimize.csv")))

	### Run the nonBayes
	cyclic_optimize <- model_fit(spi_input = cyclic_init, iter = 2000, engine = "optimize", output_dir = write_output_path)
	save(cyclic_optimize, file = file.path(write_output_path, paste0(short_name, "_cyclic_optimize.rda")))


	write.csv(x=data.frame(progress = c("optimize_complete")), file=file.path(write_output_path, paste0(short_name, "_optimize_complete.csv")))


	### Run the full model
	cyclic_sample <- model_fit(spi_input = cyclic_init, n_chains = 1, cores = 1, iter = 2000, engine = "sample", output_dir = write_output_path)
	save(cyclic_sample, file = file.path(write_output_path, paste0(short_name, "_cyclic.rda")))


	write.csv(x=data.frame(progress = c("sample_complete")), file=file.path(write_output_path, paste0(short_name, "_sample_complete.csv")))


	### Return the results
	return(cyclic_sample)
}

### Split and run for each site, then combine
cyclic_gauge_long <- fitting_df %>%
	split(.$short_name) %>% 
	future_map(~ fit_each(.x, knot_loc = knot_loc, write_output_path = write_output_path, short_name = .x$short_name[1]), .progress = FALSE) 

### Save everything
save(cyclic_gauge_long, file = file.path(write_output_path, "cyclic_gauge_long.rda"))


