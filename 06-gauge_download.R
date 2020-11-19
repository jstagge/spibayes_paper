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

### To save in SVG
require(svglite)
require(viridis)
require(ggthemes)

### For Dates
require(lubridate)
### Automatically repel text labels
require(ggrepel)
### For mapping
require(rnaturalearth)
require(sf)
library(maps)
## To access GHCND
require(rnoaa)

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
write_output_path <- file.path(output_path, "gauge_ghcnd")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(output_path, "figures/gauge_ghcnd")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

###################################################################
###  Set list of precipitation gauges
###################################################################
### You have to get a token from this website http://www.ncdc.noaa.gov/cdo-web/token and put it here
options(noaakey = "YgQSbiVrRTvNcwJiIehIeXkCtSYCozlZ")

### Set the list of gauges to download for paper
gauge_list <- c("RQW00011641", "USC00450008", "USC00517312", "USW00094728", "USW00093134", "USC00334979", "USC00419522", "USC00215638", "USC00090140")


###################################################################
###  Download station metadata
###################################################################
### Download all available stations
stations <- ghcnd_stations()

### Subset to only the stations we are interested in
stations_df <- stations %>%
	filter(element == "PRCP") %>%
	filter(id %in% gauge_list) %>%
	arrange(state) %>%
	as.data.frame()

### Check the results
stations_df

stations_df <- stations_df %>%
	mutate(short_name = c("USC", "ALB", "PAO", "MOR", "NYC", "MAR", "SAN", "WAX", "ABD"))

### Convert stations to spatial object for future work
stations_sf <- st_as_sf(stations_df, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")

###################################################################
###  Plot a map of the stations
###################################################################
### Download countries and states
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
world <- ne_countries(scale = "medium", returnclass = "sf")
head(states)

p <- ggplot(data = world) %>%
	+ geom_sf(fill = "grey80") %>%
    + geom_sf(data = states, fill = NA, colour="grey60") %>%
	+ geom_sf(fill = NA) %>%
	+ geom_sf(data = stations_sf, size = 4, shape = 23, fill = "#bd0026") %>%
    + coord_sf( crs = 4326, xlim = c(-165, -60), ylim = c(15, 55), expand=FALSE) %>%
	+ geom_text_repel(data = stations_df, aes(x = longitude, y = latitude, label = short_name), fontface = "bold", nudge_x = c(1, -1.5, 2, 2, -1), nudge_y = c(0.25, -0.25, 0.5, 0.5, -0.5)) %>%
#	+ geom_text(data = stations_df, aes(x = longitude, y = latitude, label = short_name), fontface = bold", nudge_x = -5, nudge_y = 1) %>%
	+ theme_bw(8) %>%
	+ scale_x_continuous(name = "Longitude", breaks = seq(-180,100, 10)) %>%
	+ scale_y_continuous(name = "Latitude", breaks = seq(-180,100, 10))
p

### Save Figure
ggsave(file.path(write_figures_path, "map_wgs84.png"), p, width =5, height = 3, dpi = 300)
ggsave(file.path(write_figures_path, "map_wgs84.pdf"), p, width =5, height = 3)
ggsave(file.path(write_figures_path, "map_wgs84.svg"), p, width =5, height = 3)

### Transform to albert equal area projected coordinates
target_crs <- "+proj=aea +lat_1=25 +lat_2=50 +lon_0=-100"
#target_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs"
crop_window <-  st_sfc(st_point(c(-151, -3)), st_point(c(-41, 45)), crs = 4326) %>%
	  sf::st_transform(target_crs)
crop_window_trans <- st_coordinates(crop_window)

### Plot in transformed coordinates
p <- ggplot(data = world %>% st_transform(target_crs)) %>%
	+ geom_sf(fill = "grey80") %>%
    + geom_sf(data = states, fill = NA, colour="grey60") %>%
	+ geom_sf(fill = NA) %>%
    + geom_sf(data = stations_sf, size = 4, shape = 23, fill = "#bd0026") %>%
	+ geom_sf_text(data = stations_sf, aes(label = short_name), fontface = "bold", nudge_x = -400000, nudge_y = 400000) %>%
    + coord_sf(crs = target_crs,  xlim = unlist(crop_window_trans[, 1]), ylim = unlist(crop_window_trans[, 2]),  expand=FALSE) %>%
	+ theme_bw(8) %>%
	+ scale_x_continuous(name = "Longitude", breaks = seq(-180,100, 10)) %>%
	+ scale_y_continuous(name = "Latitude", breaks = seq(-150,100, 10))
p

### Save Figure
ggsave(file.path(write_figures_path, "map_albers.png"), p, width =5, height = 2.5, dpi = 300)
ggsave(file.path(write_figures_path, "map_albers.pdf"), p, width =5, height = 2.5)
ggsave(file.path(write_figures_path, "map_albers.svg"), p, width =5, height = 2.5)

###################################################################
###  Download precipitation data
###################################################################
precip_data_path <- file.path(data_path, "ghcnd/individual")
dir.create(precip_data_path, recursive=TRUE, showWarnings = FALSE)

### Loop through and download all precipitation data
for(j in seq(1, dim(stations_df)[1])){

	### Extract the gauge id
	id_j <- stations_df$id[[j]]

	### Download the data
	precip_j <- meteo_tidy_ghcnd(stationid = id_j, var = "prcp", keep_flags = TRUE)

	### Rename to clarify the units
	precip_j <- precip_j %>%
		mutate(prcp_mm = as.numeric(prcp) * 0.1) %>%
		select(id, date, prcp_mm, mflag_prcp, qflag_prcp, sflag_prcp)	
	
	### Save each station
	write.csv(precip_j, file = file.path(precip_data_path, paste0(stations_df$short_name[[j]], "_", id_j, ".csv")), row.names = FALSE)

	if (j == 1){
		precip_df <- precip_j
	} else {
		precip_df <- bind_rows(precip_df, precip_j)
	}
}

###################################################################
###  Rework the flags
###################################################################

### Read in flags
mflag_df <- read_csv(file.path(data_path, "ghcnd/ghcnd_mflags.csv"), col_names = c("mflag_prcp", "meas_code"))
sflag_df <- read_csv(file.path(data_path, "ghcnd/ghcnd_sflags.csv"), col_names = c("sflag_prcp", "source_code"))
qflag_df <- read_csv(file.path(data_path, "ghcnd/ghcnd_qflags.csv"), col_names = c("qflag_prcp", "qual_code"))

### Add the flags to the dataframe
precip_df <- precip_df %>%
	left_join(mflag_df)%>%
	left_join(sflag_df) %>%
	left_join(qflag_df) %>%
	select(-mflag_prcp, -qflag_prcp, -sflag_prcp)

### Check the results
unique(precip_df$source_code)

### A few questionable values in the early 1900s for USC, Hawaii, and Aberdeen
precip_df %>% drop_na(qual_code) %>% as.data.frame()

### A lot of trace measures
unique(precip_df$meas_code)

###########################################################################
###  Add more processing for the id column
###########################################################################
precip_df <- precip_df %>%
	left_join(select(stations_df, id, short_name) ) %>%
	select(id, short_name, date, prcp_mm, meas_code, source_code, qual_code) %>%
	mutate(id = factor(id)) %>%
	mutate(short_name = factor(short_name))

###################################################################
###  Perform some summary work
###################################################################
### Loop through and download all precipitation data
for(j in seq(1, dim(stations_df)[1])){

	### Extract the gauge id and Download the data
	id_j <- stations_df$id[[j]]
	clim_j <- meteo_tidy_ghcnd(stationid = id_j)

	if (j == 1){
		clim_df <- clim_j
	} else {
		clim_df <- bind_rows(clim_df, clim_j)
	}
}
### Remove any column that has no values
clim_df <- clim_df[colSums(!is.na(clim_df)) > 0]

### Calculate annual sum of precip and snow
annual_sum <- clim_df %>%
	pivot_longer(c(-id, -date), names_to = "measure", values_to = "value") %>%
	mutate(year = lubridate::year(date)) %>%
	filter(measure == "prcp" | measure == "snow") %>%
	group_by(id, year, measure) %>%
	summarise(sum_annual = sum(value, na.rm=TRUE), n = sum(!is.na(value))) %>%
	filter(n > 360)

### Calculate average (across all years) annual precip and snowfall in mm, need to convert prcp from tenths of mm into mm
annual_sum_summary <- annual_sum %>%
	select(-n) %>%
	pivot_wider(names_from = measure, values_from = sum_annual) %>%
	ungroup() %>%
	mutate(prcp_mm = prcp*0.1) %>%
	rename(snow_mm = snow) %>%
	group_by(id) %>%
	summarize(prcp_mm = mean(prcp_mm, na.rm=TRUE), snow_mm = mean(snow_mm, na.rm=TRUE))

### Calculate mean climate values for gauges. Convert temp from tenths degrees
mean_clim <- clim_df %>%
	mutate(tavg = as.numeric(tavg)/10) %>%
	mutate(tmin = as.numeric(tmin)/10) %>%
	mutate(tmax = as.numeric(tmax)/10) %>%
	mutate(yup = case_when(
		!is.na(tavg) ~ tavg,
		is.na(tavg) & !is.na(tmax) ~ 0.5*(tmax + tmin),
		TRUE ~ NA_real_)) %>%
	group_by(id) %>%
	summarise(tmean = mean(tavg, na.rm=TRUE), tsun=mean(tsun, na.rm=TRUE), awnd = mean(awnd, na.rm=TRUE), tmax=mean(tmax, na.rm=TRUE), tmin=mean(tmin, na.rm=TRUE), tavg = mean(yup, na.rm=TRUE), tsun = mean(tsun, na.rm=TRUE), acsh = mean(acsh, na.rm=TRUE), acmh = mean(acmh, na.rm=TRUE))	

mean_clim

### Remove temperatures other than the daily mean
mean_clim <- mean_clim %>% 
	select(id, tavg, tsun, acsh) %>%
	rename(tavg_c = tavg) %>%
	rename(daily_sun_min = tsun) %>%
	rename(cloud_percent = acsh)

mean_clim

### Combine with station data
stations_df <- stations_df %>% 
	select(-element, -wmo_id, -gsn_flag) %>%
	left_join(annual_sum_summary) %>%
	left_join(mean_clim) %>%
	arrange(short_name)


###########################################################################
###  Save results from as an RDS file for analysis
###########################################################################

### Save results for next step
save(precip_df, stations_df, file = file.path(write_output_path, "ghcnd_data.rda"))

write.csv(stations_df, file = file.path(write_output_path, "stations_df.csv"), row.names = FALSE)





