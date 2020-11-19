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
mle_synth_stat_long <- mle_gamma_fit(jdate = synth_stat_df$jdate, values = synth_stat_df$precip, n = 3000, method = "param", ci = 0.95)

### Create WMO 30 year climate normal
short_data <- synth_stat_df %>%
	filter(year >= 1961 & year <= 1990)

mle_synth_stat_short <- mle_gamma_fit(jdate = short_data$jdate, values = short_data$precip, n = 3000, method = "param", ci = 0.95)

###########################################################################
###  Plot the results
###########################################################################
true_plot <- true_param_stat %>% 
	gather("param", "value", -jdate) %>%
	mutate(origin = "true") %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta")))

### Plot the MLE result
mle_plot <- mle_synth_stat_long$estimate %>%
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
###  Compare Long and Short
###########################################################################

### Plot the MLE result
mle_plot <- mle_synth_stat_long$estimate %>%
	mutate(length = "long") %>%
	bind_rows(mle_synth_stat_short$estimate %>% mutate(length = "short")) %>%
	pivot_wider(names_from = ci, values_from = value) %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta")))


plot_ribbon <- mle_plot %>%
	select(-estimate) %>%
	rename(facet = length) %>%
	mutate(line = factor(facet, levels = c("True", "long", "short"), labels = c("True", "Full\n(1920-2019)", "WMO\n(1961-1990)"))) %>%
	mutate(facet = factor(facet, levels = c("long", "short"), labels = c("Full (1920-2019)", "WMO (1961-1990)"))) 

true_line <- true_plot %>%
	mutate(facet = "long") %>%
	bind_rows(true_plot %>% mutate(facet = "short")) %>%
	mutate(line = "True") %>% 
	select(-origin) 

plot_line <- mle_plot %>%
	rename(value = estimate) %>%
	select(-lower_ci, -upper_ci) %>%
	rename(line = length) %>%
	mutate(facet = line) %>%
	bind_rows(true_line) %>%
	mutate(facet = factor(facet, levels = c("long", "short"), labels = c("Full (1920-2019)", "WMO (1961-1990)"))) %>%
	mutate(line = factor(line, levels = c("True", "long", "short"), labels = c("True", "Full\n(1920-2019)", "WMO\n(1961-1990)")))

#rm(true_line)

p <- ggplot(plot_line %>% filter(param != "Dispersion" & param != "Rate"), aes(x=jdate)) %>%
	+ geom_ribbon(data = plot_ribbon %>% filter(param != "Dispersion" & param != "Rate"), aes(ymin = lower_ci, ymax=upper_ci, fill = line, colour = line), alpha=0.3, size = 0.1) %>%
	+ geom_line(aes(y=value, colour = line))  %>%
	+ facet_grid(param ~ facet, scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.minor = element_blank()) %>%
	+ scale_colour_manual(name = "", values = c("black","#66c2a5", "#fc8d62"), drop = FALSE) %>%
	+ scale_fill_manual(name = "", values = c("black","#66c2a5", "#fc8d62"), drop = FALSE) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(0,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Parameter Value") #%>%
#	+ coord_cartesian(xlim=c(1, 365))

p

 
### Save Figure
ggsave(file.path(write_figures_path, "mle_fit_stat_ref_comparison.png"), p, width =5, height = 7, dpi = 300)
ggsave(file.path(write_figures_path, "mle_fit_stat_ref_comparison.pdf"), p, width =5, height = 7)
ggsave(file.path(write_figures_path, "mle_fit_stat_ref_comparison.svg"), p, width =5, height = 7)


###########################################################################
###  Save results from as an RDS file for analysis
###########################################################################

### Save results for next step
save(mle_synth_stat_long, mle_synth_stat_short, file = file.path(write_output_path, "mle_synth_stat.rda"))

