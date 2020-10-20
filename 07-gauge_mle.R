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
require(mgcv)
require(rstan)

### To save in SVG
require(svglite)
require(viridis)
require(ggthemes)

### Packages for spi
require(fitdistrplus)
require(lubridate)

select <- dplyr::select

theme_set(theme_classic(8))

### Set up  number of available cores
avail_cores <- parallel::detectCores()
options(mc.cores = avail_cores)


disp_to_shape <- function(disp){
	1/ exp(-7 + log(1 + exp(disp)))
}

logodds_to_p <- function(logodds){
	exp(logodds)/(1+exp(logodds))
}

###########################################################################
## Set the Paths
###########################################################################
### Set here path 
here_path <- here::here()

### Path for Data and Output	
data_path <- file.path(here_path, "./data")
output_path <- file.path(here_path, "./output")

### Set up output folders
write_output_path <- file.path(output_path, "cyclic_stat")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(output_path, "figures/cyclic_stat")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

###########################################################################
## Load synthetic precipitation time series
###########################################################################
load(file.path(output_path, "synth_precip/true_param.RData"))
load(file.path(output_path, "synth_precip/synth_df.RData"))
load(file.path(output_path, "mle_fit/mle_synth_stat.rda"))
ls()

###########################################################################
###  Create knots
###########################################################################
### Create knots
n_knots_jdate <- 15

knot_loc <- list(jdate = seq(1,365,length=n_knots_jdate))

###########################################################################
###  Create a smooth basis spline for demonstration
###########################################################################

### Create a dataframe with days from 1 to 366 for later and for demonstration
demo_df <- expand.grid(jdate = seq(1,365,0.25))

### Create cyclic spline basis
spline_orig <- mgcv::smoothCon(s(jdate, bs="cc", k = n_knots_jdate), data=demo_df, knots = knot_loc, null.space.penalty = TRUE)

### Create cyclic spline basis with absorbed constraint
spline_reparam <- mgcv::smoothCon(s(jdate, bs="cc", k = n_knots_jdate), data=demo_df, knots = knot_loc, absorb.cons=TRUE, null.space.penalty = TRUE)

### Extract the matrices for basis and penalty term
X_orig <- spline_orig[[1]]$X

### Reparameterize the penalty matrix
s_reparam <- spline_reparam[[1]]$S[[1]]

### Reparameterize both using the QR decomposition following Wood 
### Where Z is the Q matrix without the first column, used to reparameterize
C <- rep(1, nrow(X_orig)) %*% X_orig
qrc <- qr(t(C))
Z <- qr.Q(qrc,complete=TRUE)[,(nrow(C)+1):ncol(C)]

### Calculate reparameterized matrices for basis and penalty
X_reparam <- X_orig%*%Z

### Create an object to plot
demo_plot <- data.frame(demo_df, X_orig)
demo_plot <- demo_plot %>%
	pivot_longer(-jdate, names_to = "id", values_to = "basis") %>%
	mutate(id = as.numeric(sub('.', '', id)))

### Basis plot for demonstration purposes
p <- ggplot(demo_plot, aes(x=jdate, y=basis, colour = id, group = id)) %>%
	+ geom_line() %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(1,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Basis Function") %>%
	+ scale_colour_gradientn(colours = tableau_color_pal(palette = "Classic Cyclic")(13)) %>%
	+ theme(legend.position = "none")%>%
	+ coord_cartesian(xlim=c(1, 365))

p

### Save Figure
ggsave(file.path(write_figures_path, "basis_fig.png"), p, width =5, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "basis_fig.pdf"), p, width =5, height = 4.5)
ggsave(file.path(write_figures_path, "basis_fig.svg"), p, width =5, height = 4.5)


### Calculate the reparameterized basis
reparam_plot <- data.frame(demo_df, X_reparam)
reparam_plot <- reparam_plot %>%
	pivot_longer(-jdate, names_to = "id", values_to = "basis") %>%
	mutate(id = as.numeric(sub('.', '', id)))


### Quick plot for demonstration purposes
p <- ggplot(reparam_plot, aes(x=jdate, y=basis, colour = id, group = id)) %>%
	+ geom_line() %>%
	+ scale_x_continuous(name = "Julian Date", breaks=seq(1,365, 30), expand=c(0,0)) %>%
	+ scale_y_continuous(name = "Basis Function") %>%
	+ scale_colour_gradientn(colours = tableau_color_pal(palette = "Classic Cyclic")(13)) %>%
	+ theme(legend.position = "none")%>%
	+ coord_cartesian(xlim=c(1, 365))

p

### Save Figure
ggsave(file.path(write_figures_path, "basis_fig_reparam.png"), p, width =5, height = 4.5, dpi = 300)
ggsave(file.path(write_figures_path, "basis_fig_reparam.pdf"), p, width =5, height = 4.5)
ggsave(file.path(write_figures_path, "basis_fig_reparam.svg"), p, width =5, height = 4.5)


###########################################################################
###  Calculate basis for data and preprocess to find initial value estimates from MLE
###########################################################################
### Cut to only 365 days
fitting_df <- synth_stat_df %>%
	filter(jdate <= 365) %>%
	drop_na(precip) %>%
	select(date, jdate, year, precip, zero)

### Run initial basis function from spibayes
cyclic_init <- pre_proc(data = fitting_df, type = "cyclic", knot_loc = knot_loc)

###########################################################################
###  Check initial estimates
###########################################################################
### Check initial values
cyclic_init$input$b_0_init
cyclic_init$input$b_init
cyclic_init$input$lambda_init

### Estimate parameter values from model using initial beta estimates
newdata_df <- expand.grid(jdate = seq(1,365,1))
init_est <- predict_vals(cyclic_init, newdata = newdata_df)

### Quick plot of Initial Values
month_breaks <- c(yday(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month")), 365)
month_labels <- c(as.character(month(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month"), label=TRUE)), "Jan")

plot_line <- init_est$estimate$gamma %>%
	#select(jdate, mean) %>%
	mutate(line = "Estimate") %>%
	bind_rows( true_param_stat %>% mutate(line = "True")) %>%
	mutate(line = factor(line, levels = c("True", "Estimate")))

plot_ribbon <- init_est$marginal$mean %>%
	mutate(ymin = exp(mean - qnorm(0.975) * sigma_mean), ymax = exp(mean + qnorm(0.975)*sigma_mean)) %>%
	mutate(fill = "95% CI")

p <- ggplot(plot_line, aes(x=jdate)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = ymin, ymax = ymax, fill = fill), alpha=0.2) %>%
	+ geom_line(aes(y=mean, colour = line)) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_colour_manual(name = "Mean", values = c("red", "black")) %>%
	+ scale_fill_manual(name = NULL, values = c("grey70")) %>%
	+ scale_y_continuous(name="Mean") %>%
	+ coord_cartesian(xlim=c(1,365))
### Plot
p

### Shape parameter
plot_ribbon <- init_est$marginal$disp %>%
	mutate(ymin = disp - qnorm(0.975) * sigma_disp, ymax = disp + qnorm(0.975)*sigma_disp) %>%
	mutate(ymin = disp_to_shape(ymin), ymax = disp_to_shape(ymax)) %>%
	mutate(fill = "95% CI")

p <- ggplot(plot_line, aes(x=jdate)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = ymin, ymax = ymax, fill = fill), alpha=0.2) %>%
	+ geom_line(aes(y=shape, colour = line)) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>%
	+ scale_colour_manual(name = "Shape", values = c("red", "black")) %>%
	+ scale_fill_manual(name = NULL, values = c("grey70")) %>%
	+ scale_y_continuous(name="Shape") %>%
	+ coord_cartesian(xlim=c(1,365))
### Plot
p


### Theta parameter
plot_line <- init_est$estimate$theta %>%
	#select(jdate, mean) %>%
	mutate(line = "Estimate") %>%
	bind_rows( true_param_stat %>% mutate(line = "True")) %>%
	mutate(line = factor(line, levels = c("True", "Estimate")))

plot_ribbon <- init_est$marginal$theta %>%
	mutate(ymin = theta - qnorm(0.975) * sigma_theta, ymax = theta + qnorm(0.975)*sigma_theta) %>%
	mutate(ymin = exp(ymin)/(1+exp(ymin)), ymax = exp(ymax)/(1+exp(ymax))) %>%
	mutate(fill = "95% CI")

p <- ggplot(plot_line, aes(x=jdate)) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin = ymin, ymax = ymax, fill = fill), alpha=0.2) %>%
	+ geom_line(aes(y=theta, colour = line)) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>%
	+ scale_colour_manual(name = "Theta", values = c("red", "black")) %>%
	+ scale_fill_manual(name = NULL, values = c("grey70")) %>%
	+ scale_y_continuous(name="Theta (Zero Precip Proportion)") %>%
	+ coord_cartesian(xlim=c(1,365))
### Plot
p

###########################################################################
###  Run the model to obtain a posterior mode (penalized maximum likelihood) estimate.
###  Without the full Bayesian
###########################################################################
### Run the model
cyclic_result <- model_fit(spi_input = cyclic_init, iter = 2000, engine = "optimize", output_dir = "./output/cyclic_stat")

### Estimate parameter values from model using initial beta estimates
newdata_df <- expand.grid(jdate = seq(1,365,1))
param_est <- predict_vals(cyclic_result, newdata = newdata_df)

head(param_est$gamma)

### Quick check, confirm that optimize is equivalent to MLE (init estimate)
p <- ggplot(param_est$gamma, aes(x=jdate, y = mean)) %>%
	+ geom_line(data = true_param_stat, aes(y=mean), colour = "red") %>%
	+ geom_line(size=1) %>%
	+ geom_line(data = init_est$gamma, aes(y=mean), colour = "white", linetype = "dashed") %>%
	+ theme_classic()%>%
	+ scale_x_continuous(name = "Julian Date", breaks=round(seq(1,365,length.out=20))) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_y_continuous(name="Mean") %>%
	+ coord_cartesian(xlim=c(1,365))
### Plot
p

### Check of exact fitted model
param_est <- predict_vals(cyclic_result)

head(param_est$gamma)

### Raster plot
ggplot(param_est$gamma, aes(x=jdate, y=year, fill=mean)) + geom_raster() + scale_fill_viridis()

### Check everywhere
newdata_df <- expand.grid(jdate = seq(1,365,1), year = seq(min(fitting_df$year, na.rm=TRUE), max(fitting_df$year, na.rm=TRUE), 1))

param_est <- predict_vals(cyclic_result, newdata = newdata_df)

ggplot(param_est$gamma, aes(x=jdate, y=year, fill=mean)) + geom_raster() + scale_fill_viridis()


###########################################################################
###  Fit the full Bayesian model
###########################################################################

### Run the model
model_cyclic <- model_fit(spi_input = cyclic_init, n_chains = 1, cores = 1, iter = 2500, engine = "sample", output_dir = "./output/cyclic_stat")


###########################################################################
###  Save results from as an RDS file for analysis
###########################################################################
### Read in the model run information
model_read <- rstan::read_stan_csv(model_cyclic$model_fit$output_files())

#model_read <- rstan::read_stan_csv("/media/data/Documents/work_folder/projects_research/code/spibayes_paper/output/cyclic_stat/cyclic_ti-202010061934-1-78c360.csv")

#model_read <- rstan::read_stan_csv("/media/data/Documents/work_folder/projects_research/code/spibayes_paper/output/cyclic_stat/cyclic_ti-202010061701-1-2e4b83.csv")

### Save results for next step
save(model_cyclic, file = file.path(write_output_path, "model_cyclic.rda"))
save(model_read, file = file.path(write_output_path, "model_read.rda"))

#load(file.path(write_output_path, "model_cyclic.rda"))
#load(file.path(write_output_path, "model_read.rda"))

###########################################################################
###  Check run
###########################################################################
### plot each trace
plot(model_read, plotfun = "trace", pars = "lambda_mean", inc_warmup = TRUE)
plot(model_read, plotfun = "trace", pars = "lambda_disp", inc_warmup = TRUE)
plot(model_read, plotfun = "trace", pars = "lambda_theta", inc_warmup = TRUE)

### Check the trace plots to confirm the chains converge
plot(model_read, plotfun = "trace", pars = "b_mean_jdate", inc_warmup = TRUE)
plot(model_read, plotfun = "trace", pars = "b_disp_jdate", inc_warmup = TRUE)
plot(model_read, plotfun = "trace", pars = "b_theta_jdate", inc_warmup = TRUE)

### Check chains without warmup
plot(model_read, plotfun = "trace", pars = "b_mean_jdate")
plot(model_read, plotfun = "trace", pars = "b_disp_jdate")
plot(model_read, plotfun = "trace", pars = "b_theta_jdate")

plot(model_read, plotfun = "trace", pars = "sigma_mean", inc_warmup = FALSE)
plot(model_read, plotfun = "trace", pars = "sigma_disp", inc_warmup = FALSE)
plot(model_read, plotfun = "trace", pars = "sigma_theta", inc_warmup = FALSE)


### Check the distributions of beta values
plot(model_read, show_density = TRUE, ci_level = 0.5, pars = "b_mean_jdate", fill_color = "lightblue") + theme_classic()
plot(model_read, show_density = TRUE, ci_level = 0.5, pars = "b_disp_jdate", fill_color = "lightblue") + theme_classic()
plot(model_read, show_density = TRUE, ci_level = 0.5, pars = "b_theta_jdate", fill_color = "lightblue") + theme_classic()

###########################################################################
###  Convert parameters to long form and summarize
###########################################################################
### Estimate parameter values from model using initial beta estimates
newdata_df <- expand.grid(jdate = seq(1,365,1))
param_est <- predict_vals(model_cyclic, newdata = newdata_df)

### Convert to long format
gamma_long <-  param_est$estimate$gamma %>%
	select(-chain, -iteration, -draw) %>%
	pivot_longer(
		cols = c(-jdate),
		names_to = "param",
		values_to = "value")

theta_long <-  param_est$estimate$theta %>%
	select(-chain, -iteration, -draw) %>%
	pivot_longer(
		cols = c(-jdate),
		names_to = "param",
		values_to = "value")

param_long <- bind_rows(gamma_long, theta_long)

rm(gamma_long, theta_long)

### Summarize by jdate
param_summary <- param_long %>%
	group_by(jdate, param) %>%
	summarise(median = median(value, na.rm=TRUE), q_25 = quantile(value, 0.25), q_75 = quantile(value, 0.75), q_2_5 = quantile(value, 0.025), q_97_5 = quantile(value, 0.975)) %>%
	ungroup()

### Extract an estimate of the parameter 
param_long <- param_summary %>%
	select(jdate, param, median) %>%
	rename(value = median) %>%
	mutate(source = "Model")

### Extract several confidence intervals in a format useful for plotting and storage
#param_temp <- param_summary %>%
#	select(jdate, param, q_25, q_75) %>%
#	rename(lower = q_25, upper = q_75) %>%
#	mutate(ci = "50%")

param_ci <- param_summary %>%
	select(jdate, param, q_2_5, q_97_5) %>%
	rename(lower = q_2_5, upper = q_97_5) %>%
	mutate(ci = "95%") %>%
#	bind_rows(param_temp)
	mutate(source = "Model") %>%
	mutate(ci = "Estimate")

#rm(param_temp)


###########################################################################
###  Calculate 95% credible interval
###########################################################################
mean_ci <- param_est$marginal$mean  %>%
	mutate(lower = mean - qnorm(0.975) * sigma_mean, upper = mean + qnorm(0.975)*sigma_mean) %>%
	group_by(jdate) %>%
	summarise(lower = median(lower, na.rm=TRUE), upper = median(upper, na.rm=TRUE)) %>%
	mutate(lower = exp(lower), upper = exp(upper)) %>%
	ungroup() %>%
	mutate(param = "mean")

shape_ci <- param_est$marginal$disp  %>%
	mutate(lower = disp - qnorm(0.975) * sigma_disp, upper = disp + qnorm(0.975)*sigma_disp) %>%
	group_by(jdate) %>%
	summarise(lower = median(lower, na.rm=TRUE), upper = median(upper, na.rm=TRUE)) %>%
	mutate(lower = disp_to_shape(lower), upper = disp_to_shape(upper)) %>%
	ungroup() %>%
	mutate(param = "shape")

theta_ci <- param_est$marginal$theta  %>%
	mutate(lower = theta - qnorm(0.975) * sigma_theta, upper = theta + qnorm(0.975)*sigma_theta) %>%
	group_by(jdate) %>%
	summarise(lower = median(lower, na.rm=TRUE), upper = median(upper, na.rm=TRUE)) %>%
	mutate(lower = logodds_to_p(lower), upper = logodds_to_p(upper)) %>%
	ungroup() %>%
	mutate(param = "theta")

ci_ribbon <- mean_ci %>%
	bind_rows(shape_ci)  %>%
	bind_rows(theta_ci) %>%
	mutate(source = "Model") %>%
	mutate(ci = "Error")

###########################################################################
###  Process true values for plotting
###########################################################################
### Process the true values
true_param_long <- true_param_stat %>%
	pivot_longer(
		cols = c(-jdate),
		names_to = "param",
		values_to = "value") %>%
	mutate(source = "True")

###########################################################################
###  Process MLE for plotting
###########################################################################
### Process the MLE 
mle_summary <- mle_synth_stat$draws %>%
	select(-draw, -type) %>%
	pivot_longer(
		cols = c(-jdate),
		names_to = "param",
		values_to = "value") %>%
	group_by(jdate, param) %>%
	summarise(median = median(value, na.rm=TRUE), q_25 = quantile(value, 0.25), q_75 = quantile(value, 0.75), q_2_5 = quantile(value, 0.025), q_97_5 = quantile(value, 0.975)) %>%
	ungroup()

mle_long  <- mle_summary %>%
	select(jdate, param, median) %>%
	rename(value = median) 

theta_temp <- mle_synth_stat$estimate %>%
	filter(ci == "estimate" & param == "theta") %>%
	select(-ci)

mle_long <- mle_long %>%
	bind_rows(theta_temp) %>%
	mutate(source = "MLE")

rm(theta_temp)

### Extract several confidence intervals in a format useful for plotting and storage
#mle_temp <- mle_summary %>%
#	select(jdate, param, q_25, q_75) %>%
#	rename(lower = q_25, upper = q_75) %>%
#	mutate(ci = "50%")

mle_ci <- mle_summary %>%
	select(jdate, param, q_2_5, q_97_5) %>%
	rename(lower = q_2_5, upper = q_97_5) %>%
#	bind_rows(mle_temp) %>%
	mutate(source = "MLE") %>%
	mutate(ci = "Error") %>%
	filter(param == "mean" | param == "shape" | param == "theta") %>%
	mutate(param = factor(param, levels = c("mean", "shape", "theta"), labels = c("Mean", "Shape", "Theta"))) 	%>%
	mutate(source = factor(source, levels = c("True", "MLE", "Model")))

#rm(mle_temp)


###########################################################################
###  Plot fit
###########################################################################
### Quick plot of Initial Values
month_breaks <- c(yday(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month")), 365)
month_labels <- c(as.character(month(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month"), label=TRUE)), "Jan")


plot_line <- param_long %>%
	bind_rows(mle_long) %>%
	bind_rows(true_param_long) %>%
	filter(param == "mean" | param == "shape" | param == "theta") %>%
	mutate(param = factor(param, levels = c("mean", "shape", "theta"), labels = c("Mean", "Shape", "Theta"))) 	%>%
	mutate(source = factor(source, levels = c("True", "MLE", "Model")))
	#mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta"))) 	


plot_ribbon <- ci_ribbon %>%
	bind_rows(param_ci)  %>%
	filter(param == "mean" | param == "shape" | param == "theta") %>%
	mutate(param = factor(param, levels = c("mean", "shape", "theta"), labels = c("Mean", "Shape", "Theta"))) 	
	#mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp", "theta"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion", "Theta"))) 	


p <- ggplot(plot_line, aes(x=jdate)) %>%
	+ geom_line(data = mle_ci, aes(y=lower), colour = "#a6cee3", alpha = 0.8, size = 0.27) %>%
	+ geom_line(data=mle_ci, aes(y=upper), colour = "#a6cee3", alpha = 0.8, size = 0.27) %>%
	+ geom_ribbon(data = plot_ribbon, aes(ymin=lower, ymax=upper, fill = ci), alpha = 0.6) %>%
	#+ geom_ribbon(data = mle_ci, aes(ymin=lower, ymax=upper, colour = source), fill = NULL, alpha = 0.5) %>%
	+ geom_line(aes(y=value, colour = source), alpha = 0.8) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_colour_manual(name = "Estimate", values = c("grey20", "#1f78b4", "#e41a1c")) %>%
	#+ scale_linetype_manual(name = "Estimate", values=c("solid", "twodash", "solid")) %>%
	+ scale_fill_manual(name = "95% CI", values = c("grey70", "grey30")) %>%
	+ scale_y_continuous(name="Mean") %>%
	+ coord_cartesian(xlim=c(1,365)) %>%
	+ facet_grid(param ~ ., scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p


### Save Figure
ggsave(file.path(write_figures_path, "stat_model_vs_mle.png"), p, width =4.5, height = 5.5, dpi = 300)
ggsave(file.path(write_figures_path, "stat_model_vs_mle.pdf"), p, width =4.5, height = 5.5)
ggsave(file.path(write_figures_path, "stat_model_vs_mle.svg"), p, width =4.5, height = 5.5)



### Try with facets
plot_line_facet <- plot_line %>%
	filter(source == "MLE" | source == "True") %>%
	mutate(facet = "MLE") %>%
	bind_rows(plot_line %>% filter(source == "Model" | source == "True") %>%
	mutate(facet = "Model"))

plot_ribbon_facet <- plot_ribbon %>% 
	mutate(facet = "Model") %>%
	bind_rows(mle_ci %>% mutate(facet="MLE"))

p <- ggplot(plot_line_facet, aes(x=jdate)) %>%
	+ geom_ribbon(data = plot_ribbon_facet, aes(ymin=lower, ymax=upper, fill = ci), alpha = 0.3) %>%
	#+ geom_ribbon(data = mle_ci, aes(ymin=lower, ymax=upper, colour = source), fill = NULL, alpha = 0.5) %>%
	+ geom_line(aes(y=value, colour = source)) %>%
	#+ geom_line(data = mle_ci, aes(y=lower, colour = source), alpha = 0.5) %>%
	#+ geom_line(data=mle_ci, aes(y=upper, colour = source), alpha = 0.5) %>%
	+ scale_x_continuous(name = "Julian Date", breaks=month_breaks, expand = c(0,0), sec.axis = sec_axis(~ . + 0, breaks = month_breaks, labels = month_labels)) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_colour_manual(name = "Estimate", values = c("grey20", "#377eb8", "#e41a1c")) %>%
	#+ scale_linetype_manual(name = "Estimate", values=c("solid", "twodash", "solid")) %>%
	+ scale_fill_manual(name = "95% CI", values = c("grey70", "grey30")) %>%
	+ scale_y_continuous(name="Mean") %>%
	+ coord_cartesian(xlim=c(1,365)) %>%
	+ facet_grid(param ~ facet, scales="free_y") %>%
	+ theme_bw(8) %>%
	+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


### Save Figure
ggsave(file.path(write_figures_path, "stat_model_vs_mle_facet.png"), p, width =6.5, height = 5.5, dpi = 300)
ggsave(file.path(write_figures_path, "stat_model_vs_mle_facet.pdf"), p, width =6.5, height = 5.5)
ggsave(file.path(write_figures_path, "stat_model_vs_mle_facet.svg"), p, width =6.5, height = 5.5)





###########################################################################
###  Save model results
###########################################################################

save(fit_2, b_full_rate, b_full_shape, file = "../output/mymodel.rda")



