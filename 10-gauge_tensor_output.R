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

### Loop through each site
for (j in seq(1, length(site_list))){

short_name_j <- site_list[j]

cat(j)
cat(short_name_j)

### Find all files within folder
#write_output_path <- "/run/media/jhstagge/Seagate Backup Plus Drive/spibayes_paper/output/gauge_tensor"
site_folder <- file.path(write_output_path, short_name_j)

### Find the biggest. Could do a newest
site_files <- Sys.glob(file.path(site_folder, "*.csv"))
site_files_info <- file.info(site_files)
file_j <- rownames(site_files_info)[which.max(site_files_info$size)]
#order(site_files_info$mtime)

load(file.path(site_folder, paste0(short_name_j, "_knot_loc.rda")))
load(file.path(site_folder, paste0(short_name_j, "_tensor.rda")))

### Read in precip to get dates
precip_j <- precip_df %>% filter(short_name == short_name_j) %>% mutate(year = year(date))

year_list <- seq(min(precip_j$year),max(precip_j$year),1)
date_seq <- seq(as.Date(paste0(min(year_list), "-01-01")), as.Date(paste0(max(year_list), "-12-31")), by = "1 day")

dates_df <- data.frame(date = date_seq) %>%
	mutate(year = year(date), jdate = yday(date))

year_seq <- dates_df$year
jdate_seq <- dates_df$jdate
date_seq <- dates_df$date
time_seq <- as.numeric(date_seq - as.Date("1900-01-01"))

rm(precip_j)

cat("Line 123")

### Check the number of succesful draws
### I don't have a good way to do this at the moment, so just filter to mean !
param_est <- predict_vals(tensor_sample, newdata = expand.grid(jdate = 1, year = 2000), saved_model = file_j)

param_est <- param_est$estimate$gamma %>%
  filter(mean != 1) 
draws_seq <- unique(param_est$draw)
rm(param_est)

cat("Line 134")
###########################################################################
## Create dimensions and variables
###########################################################################
mis_val <-  -9.99e+08

# Define temporal dimensions
#date_seq <- seq(as.Date(paste0(year_j, "-01-01")), as.Date(paste0(year_j, "-12-31")), by = "1 day")
time_seq <- as.numeric(date_seq - as.Date("1900-01-01"))

# Define other dimensions
dim_date <- ncdim_def( "time", units="days since 1900-01-01", vals=time_seq, longname="time", unlim=TRUE)
dim_draw <- ncdim_def( "draw", units="", vals=draws_seq, longname="Bayesian draw", create_dimvar=TRUE)
dim_marginal <- ncdim_def( "marginal", units="", vals=seq(1,6), longname="Marginal", create_dimvar=TRUE)


cat("Line 150")

### Define variables
mean_var = ncvar_def(name="mean", units="", dim=list(dim_date, dim_draw), missval=mis_val, longname="Mean parameter from gamma distribution", prec="float", compression=4)
shape_var = ncvar_def(name="shape", units="", dim=list(dim_date, dim_draw), missval=mis_val, longname="Shape parameter from gamma distribution", prec="float", compression=4)
scale_var = ncvar_def(name="scale", units="", dim=list(dim_date, dim_draw), missval=mis_val, longname="Scale parameter from gamma distribution", prec="float", compression=4)
theta_var = ncvar_def(name="theta", units="", dim=list(dim_date, dim_draw), missval=mis_val, longname="Theta parameter from binomial distribution ", prec="float", compression=4)

cat("Line 158")


jdate_var = ncvar_def(name="jdate", units="", dim=list(dim_date), missval=mis_val, longname="Julian date", prec="integer", compression=4)
year_var = ncvar_def(name="year", units="", dim=list(dim_date), missval=mis_val, longname="Year", prec="integer", compression=4)

### Define marginal variables
mean_marginal_var = ncvar_def(name="mean", units="", dim=list(dim_date, dim_draw, dim_marginal), missval=mis_val, longname="Mean parameter marginals", prec="float", compression=4)
disp_marginal_var = ncvar_def(name="disp", units="", dim=list(dim_date, dim_draw, dim_marginal), missval=mis_val, longname="Dispersion parameter marginals", prec="float", compression=4)
theta_marginal_var = ncvar_def(name="theta", units="", dim=list(dim_date, dim_draw, dim_marginal), missval=mis_val, longname="Theta parameter marginals", prec="float", compression=4)

cat("Line 169")

###########################################################################
## Create estimate file
###########################################################################

### Begin to create a new domain file
nc_file <- file.path(site_folder, paste0(short_name_j, "_estimate.nc"))

# Define file
nc = nc_create(nc_file, list(mean_var, shape_var, scale_var, theta_var, jdate_var, year_var))  #

cat("Line 181")

### Put the variables
ncvar_put( nc, "jdate", vals = jdate_seq, start=c(1), count=length(jdate_seq), verbose=FALSE )
ncvar_put( nc, "year", vals = year_seq, start=c(1), count=length(year_seq), verbose=FALSE )

cat("Line 187")

# Add global attributes
ncatt_put(nc, 0, "Title", "Bayesian SPI Estimates from Tensor Product")
ncatt_put(nc, 0, "Site", short_name_j)
#ncatt_put(nc, 0, "Conventions", "CF-1.6")
#ncatt_put(nc, 0, "Source1", "Generated from VICGlobal v 1.6")
#ncatt_put(nc, 0, "Source2", "https://doi.org/10.5281/zenodo.3972341")
ncatt_put(nc, 0, "Note1", "Generated by Jim Stagge")	
cur_time = Sys.time()
ncatt_put(nc, 0, "Note2", paste("Modified", format(Sys.time(), "%Y%m%d_%H%M%S"), "on host", Sys.info()[4]))

nc_close(nc)

cat("Line 201")

###########################################################################
## Create marginal file
###########################################################################
### Begin to create a new domain file
nc_marginal_file <- file.path(site_folder, paste0(short_name_j, "_marginal.nc"))

# Define file
nc_marginal = nc_create(nc_marginal_file, list(mean_marginal_var, disp_marginal_var, theta_marginal_var, jdate_var, year_var))  #

cat("Line 212")

### Put the variables
ncvar_put( nc_marginal, "jdate", vals = jdate_seq, start=c(1), count=length(jdate_seq), verbose=FALSE )
ncvar_put( nc_marginal, "year", vals = year_seq, start=c(1), count=length(year_seq), verbose=FALSE )

cat("Line 218")

# Add global attributes
ncatt_put(nc_marginal, "marginal", "comment", "Marginals presented as follows: intercept, jdate, year, tensor product, estimate, sigma", prec="char" )

ncatt_put(nc_marginal, 0, "Title", "Bayesian SPI Marginals from Tensor Product")
ncatt_put(nc_marginal, 0, "Site", short_name_j)
#ncatt_put(nc, 0, "Conventions", "CF-1.6")
#ncatt_put(nc, 0, "Source1", "Generated from VICGlobal v 1.6")
#ncatt_put(nc, 0, "Source2", "https://doi.org/10.5281/zenodo.3972341")
ncatt_put(nc_marginal, 0, "Note1", "Generated by Jim Stagge")	
cur_time = Sys.time()
ncatt_put(nc_marginal, 0, "Note2", paste("Modified", format(Sys.time(), "%Y%m%d_%H%M%S"), "on host", Sys.info()[4]))

nc_close(nc_marginal)

cat("Line 234")

### Loop year by year to prevent memory overflows
	#for (year_k in unique(year_seq)){

year_cut_start <- seq(min(year_list), max(year_list), by = 5)
year_cut_end <- year_cut_start + 4
year_cut_end[length(year_cut_end)] <- max(year_list)
year_cut_df <- data.frame(start = year_cut_start, end = year_cut_end)

	for (k in seq(1, dim(year_cut_df)[1])){ 

		start_year_k <- year_cut_df$start[k]
		end_year_k <- year_cut_df$end[k]
		cat(start_year_k)
		cat("\n")

		dates_k <- seq(as.Date(paste0(start_year_k, "-01-01")), as.Date(paste0(end_year_k, "-12-31")), by = "1 day")

cat("Line 243")

		newdata_df <- data.frame(date = dates_k, year = year(dates_k), jdate = yday(dates_k))
		param_est <- predict_vals(tensor_sample, newdata = newdata_df, saved_model = file_j)

cat("Line 243")

		### Process Gamma estimaate
		gamma_df <- param_est$estimate$gamma %>%
			filter(mean != 1) %>%	
			mutate(date = as.Date(jdate-1, origin=as.Date(paste0(year, "-01-01")))) #%>%
			#complete(date, newdata_df$date)

cat("Line 256")

		mean_mat <- gamma_df %>%
			select(date, draw, mean) %>%
			pivot_wider(names_from=draw, values_from = mean) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		shape_mat <- gamma_df %>%
			select(date, draw, shape) %>%
			pivot_wider(names_from=draw, values_from = shape) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		scale_mat <- gamma_df %>%
			select(date, draw, scale) %>%
			pivot_wider(names_from=draw, values_from = scale) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 


cat("Line 280")

		### Process Theta estimate
		theta_mat <- param_est$estimate$theta %>%
			filter(theta != 0.5) %>%	
			mutate(date = as.Date(jdate-1, origin=as.Date(paste0(year, "-01-01"))))  %>%
			select(date, draw, theta) %>%
			pivot_wider(names_from=draw, values_from = theta) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		### Figure out where to place year in netcdf
		time_index <- which(date_seq %in% dates_k)
		start_time <- min(time_index)
		count_time <- length(time_index)
		start_draw <- 1
		count_draw <- length(unique(gamma_df$draw))

cat("Line 299")

		### Add to netcdf file
		nc = nc_open(nc_file, write = TRUE)

		ncvar_put( nc, "mean", vals = mean_mat, start=c(start_time,start_draw), count=c(count_time,count_draw), verbose=FALSE )
		nc_close(nc)

cat("Line 307")

		nc = nc_open(nc_file, write = TRUE)
		ncvar_put( nc, "shape", vals = shape_mat, start=c(start_time,start_draw), count=c(count_time,count_draw), verbose=FALSE )
		nc_close(nc)

cat("Line 313")

		nc = nc_open(nc_file, write = TRUE)
		ncvar_put( nc, "scale", vals = scale_mat, start=c(start_time,start_draw), count=c(count_time,count_draw), verbose=FALSE )
		nc_close(nc)

cat("Line 319")

		nc = nc_open(nc_file, write = TRUE)
		ncvar_put( nc, "theta", vals = theta_mat, start=c(start_time,start_draw), count=c(count_time,count_draw), verbose=FALSE )

		nc_close(nc)


cat("Line 327")

		### Create empty object for marginals
		mean_mar_array <- array(NA, c(count_time, count_draw, 6))
		disp_mar_array <- array(NA, c(count_time, count_draw, 6))
		theta_mar_array <- array(NA, c(count_time, count_draw, 6))

		### Process Mean marginal
		mean_marginals <- param_est$marginal$mean %>%
			filter(draw <= count_draw)

cat("Line 338")

		mean_mar_array[,,1] <- mean_marginals %>%
			select(date, draw, mean_0) %>%
			pivot_wider(names_from=draw, values_from = mean_0) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		mean_mar_array[,,2] <- mean_marginals %>%
			select(date, draw, mean_jdate) %>%
			pivot_wider(names_from=draw, values_from = mean_jdate) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		mean_mar_array[,,3] <- mean_marginals %>%
			select(date, draw, mean_year) %>%
			pivot_wider(names_from=draw, values_from = mean_year) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		mean_mar_array[,,4] <- mean_marginals %>%
			select(date, draw, mean_tensor) %>%
			pivot_wider(names_from=draw, values_from = mean_tensor) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		mean_mar_array[,,5] <- mean_marginals %>%
			select(date, draw, mean) %>%
			pivot_wider(names_from=draw, values_from = mean) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		mean_mar_array[,,6] <- mean_marginals %>%
			select(date, draw, sigma_mean) %>%
			pivot_wider(names_from=draw, values_from = sigma_mean) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

	### Process Dispersion marginal
		disp_marginals <- param_est$marginal$disp %>%
			filter(draw <= count_draw)

		disp_mar_array[,,1] <- disp_marginals %>%
			select(date, draw, disp_0) %>%
			pivot_wider(names_from=draw, values_from = disp_0) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		disp_mar_array[,,2] <- disp_marginals %>%
			select(date, draw, disp_jdate) %>%
			pivot_wider(names_from=draw, values_from = disp_jdate) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		disp_mar_array[,,3] <- disp_marginals %>%
			select(date, draw, disp_year) %>%
			pivot_wider(names_from=draw, values_from = disp_year) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		disp_mar_array[,,4] <- disp_marginals %>%
			select(date, draw, disp_tensor) %>%
			pivot_wider(names_from=draw, values_from = disp_tensor) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		disp_mar_array[,,5] <- disp_marginals %>%
			select(date, draw, disp) %>%
			pivot_wider(names_from=draw, values_from = disp) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		disp_mar_array[,,6] <- disp_marginals %>%
			select(date, draw, sigma_disp) %>%
			pivot_wider(names_from=draw, values_from = sigma_disp) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

	### Process Theta marginal
		theta_marginals <- param_est$marginal$theta %>%
			filter(draw <= count_draw)

		theta_mar_array[,,1] <- theta_marginals %>%
			select(date, draw, theta_0) %>%
			pivot_wider(names_from=draw, values_from = theta_0) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		theta_mar_array[,,2] <- theta_marginals %>%
			select(date, draw, theta_jdate) %>%
			pivot_wider(names_from=draw, values_from = theta_jdate) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		theta_mar_array[,,3] <- theta_marginals %>%
			select(date, draw, theta_year) %>%
			pivot_wider(names_from=draw, values_from = theta_year) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		theta_mar_array[,,4] <- theta_marginals %>%
			select(date, draw, theta_tensor) %>%
			pivot_wider(names_from=draw, values_from = theta_tensor) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		theta_mar_array[,,5] <- theta_marginals %>%
			select(date, draw, theta) %>%
			pivot_wider(names_from=draw, values_from = theta) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

		theta_mar_array[,,6] <- theta_marginals %>%
			select(date, draw, sigma_theta) %>%
			pivot_wider(names_from=draw, values_from = sigma_theta) %>%
			arrange(date) %>% 
			select(-date) %>%
			as.matrix() 

cat("Line 474")


		### Add to netcdf file
		nc = nc_open(nc_marginal_file, write = TRUE)

		ncvar_put( nc, "mean", vals = mean_mar_array, start=c(start_time,start_draw, 1), count=c(count_time,count_draw, 6), verbose=FALSE )

		nc_close(nc)
cat("Line 483")

		nc = nc_open(nc_marginal_file, write = TRUE)
		ncvar_put( nc, "disp", vals = disp_mar_array, start=c(start_time,start_draw, 1), count=c(count_time,count_draw, 6), verbose=FALSE )

		nc_close(nc)
cat("Line 489")

		nc = nc_open(nc_marginal_file, write = TRUE)
		ncvar_put( nc, "theta", vals = theta_mar_array, start=c(start_time,start_draw, 1), count=c(count_time,count_draw, 6), verbose=FALSE )

		nc_close(nc)
cat("Line 495")
	}

}















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


