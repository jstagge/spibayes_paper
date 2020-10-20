
###########################################################################
###  Check results
###########################################################################





plot_line <- param_est$estimate$gamma %>% 
	select(jdate, mean, shape, iteration) %>% 
	mutate(line = "Estimate") 


%>%
	bind_rows( true_param_stat %>% mutate(line = "True")) %>%
	mutate(line = factor(line, levels = c("True", "Estimate")))

mle_synth_stat$estimate %>% filter(param == "mean" & ci == "estimate")

head(param_est$gamma)

### Quick check
p <- ggplot(param_est$gamma, aes(x=jdate, y = mean)) %>%
	+ geom_line(aes(group = iteration), alpha=0.1) %>%
	+ geom_line(data = true_param_stat, aes(y=mean), colour = "red") %>%
	#+ geom_line(data = init_est$gamma, aes(y=mean), colour = "blue") %>%
	+ geom_line(data = mle_synth_stat$estimate %>% filter(param == "mean" & ci == "estimate"), aes(y=value), colour = "green") %>%
	+ theme_classic()%>%
	+ scale_x_continuous(name = "Julian Date", breaks=round(seq(1,365,length.out=20))) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_y_continuous(name="Mean") %>%
	+ coord_cartesian(xlim=c(1,365))
### Plot
p











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










### Quick plot to test
p <- ggplot(filter(param_est, draw=="X1"), aes(x=jdate, y=year, fill=mean)) %>%
	+ geom_tile() %>%
	+ scale_fill_viridis(name = "Mean")%>%
	+ theme_bw()%>%
	+ scale_x_continuous(name = "Julian Date", breaks=round(seq(1,365,length.out=20))) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_y_continuous(name="Year", expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1,365), ylim=c(1900,2020))
### Plot
p

### Save the plot
ggsave("./output/tensor_plot.png", p, width =12.5, height = 8, dpi = 300)

### Quick plot of all draws
ggplot(param_est, aes(x=jdate)) + geom_line(aes(group=draw, y=mean), colour="grey40", alpha=0.2) + theme_classic() + facet_grid(year ~ . )

#geom_line(data= mle_fit, aes(y=mean), colour="red") + 

ggplot(param_est %>% filter(jdate == unique(param_est$jdate)[seq(1,74, by =7)]), aes(x=year)) + geom_line(aes(group=draw, y=mean), colour="grey40", alpha=0.2) + theme_classic() + facet_grid(jdate ~ . )

### Convert to long format to plot all together
param_est_long <- param_long(param_est)
param_est_long <- param_est_long %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion")))
head(param_est_long)

### Plot all draws and all parameters
p <- ggplot(filter(param_est_long, year == 1990), aes(x=jdate)) + geom_ribbon(data = mle_plot, aes(ymin = lower_ci, ymax = upper_ci), fill="grey60", alpha=0.2) + geom_line(aes(group=draw, y=value), colour="grey40", alpha=0.2) + theme_classic() + facet_grid( param ~ . , scales = "free_y") + geom_line(data= mle_plot, aes(y=estimate), colour="red")
p

### Save the plot
ggsave("./output/tensor_plot3.png", p, width =12.5, height = 8, dpi = 300)

### Plot all draws and all parameters
p <- ggplot(filter(param_est_long, year == 1970), aes(x=jdate)) + geom_ribbon(data = mle_plot, aes(ymin = lower_ci, ymax = upper_ci), fill="grey60", alpha=0.2) + geom_line(aes(group=draw, y=value), colour="grey40", alpha=0.2) + theme_classic() + facet_grid( param ~ . , scales = "free_y") + geom_line(data= mle_plot, aes(y=estimate), colour="red")
p


### Calculate some summary statistics (quantiles) to plot the spline results
param_summ <- param_summary(param_est)
param_summ <- param_summ %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion")))
head(param_summ)

plot_df <- filter(param_summ, year == 1990)
### Create a plot comparing the true value, spline with uncertainty bounds, and the MLE estimate
p <- ggplot(plot_df, aes(x=jdate))  %>%
	+ geom_ribbon(data = mle_plot, aes(ymin = lower_ci, ymax = upper_ci), fill="grey60", alpha=0.2) %>%
	+ geom_line(data = mle_plot, aes(y=estimate), colour="#66c2a5") %>%
	+ geom_ribbon(data = plot_df, aes(ymin = perc_50_lower, ymax = perc_95_upper), fill="grey70", alpha = 0.5)  %>%
	+ geom_line(aes(y = median), colour="#fc8d62", size=1) %>% 
	+ geom_line(data = filter(param_summ, year == 1960), aes(y = median), colour="black", size=1) %>% 
	+ theme_classic() %>%
	+ scale_x_continuous(name = "Julian Date", breaks=round(seq(1,365,length.out=20)), expand = c(0, 0)) %>%
	+ scale_y_continuous(name="Parameter Estimate") %>%
	+ coord_cartesian(xlim=c(1,365)) %>%
    + facet_grid(param ~., scales="free_y")

### Plot
p

### Save the plot
ggsave("./output/tensor_plot2.png", p, width =12.5, height = 8, dpi = 300)


### Quick plot to test
p <- ggplot(param_summ %>% filter(param == "Mean"), aes(x=jdate, y=year, fill=median)) %>%
	+ geom_tile() %>%
	+ scale_fill_viridis()%>%
	+ theme_bw()%>%
	+ scale_x_continuous(name = "Julian Date", breaks=round(seq(1,365,length.out=20))) %>% ### This seems to break it, putting white lines , expand = c(0, 0)
	+ scale_y_continuous(name="Year", expand=c(0,0)) %>%
	+ coord_cartesian(xlim=c(1,365), ylim=c(1900,2020)) 
### Plot
p

### Save the plot
ggsave("./output/tensor_plot4.png", p, width =12.5, height = 8, dpi = 300)



new_data <- expand.grid(jdate = round(seq(1,336,length.out = 12)), year = seq(1950,2010, 1))

### Estimate parameters for all draws using the demo (single year) basis
param <- extract_params(model_fit = model_tens, basis = basis_full, newdata = new_data)

param_summ <- param_summary(param$param_est)
head(param_summ)

plot_df <- param_summ %>% filter(param == "shape")
mle_plot <- mle_plot %>%
	right_join(expand.grid(jdate = unique(param_summ$jdate), year = c(1950,2020))) %>%
	filter( param != "theta")

p <- ggplot(param_summ, aes(x=year))  %>%
#	+ geom_line(data = mle_plot, aes( y=value), colour="#66c2a5") %>%
	+ geom_ribbon(data = param_summ, aes(ymin = perc_95_lower, ymax = perc_95_upper), fill="grey70", alpha = 0.5)  %>%
	+ geom_ribbon(data = param_summ, aes(ymin = perc_50_lower, ymax = perc_50_upper), fill="grey50", alpha = 0.5) %>% 
	+ geom_line(aes(y = median), colour="#fc8d62", size=1) %>% 
	+ theme_classic() %>%
	+ scale_x_continuous(name = "Year", breaks=round(seq(1900,2020,10)), expand = c(0, 0)) %>%
	+ scale_y_continuous(name="Parameter Estimate") %>%
	#+ coord_cartesian(xlim=c(1,365)) %>%
    + facet_grid(param ~ jdate, scales = "free_y")

### Plot
p

### Save the plot
ggsave("./output/tensor_plot5.png", p, width =12.5, height = 8, dpi = 300)





new_data <- expand.grid(jdate = seq(1,366,60), year = seq(1939,2020, 1))

### Estimate parameters for all draws using the demo (single year) basis
param <- extract_params(model_fit = model_tens, basis = basis_full, newdata = new_data)

param_est <- param$param_est

### Calculate some summary statistics (quantiles) to plot the spline results
param_summ <- param_summary(param_est)
param_summ <- param_summ %>%
	mutate(param = factor(param, levels = c("mean", "scale", "rate", "shape", "disp"), labels = c("Mean", "Scale", "Rate", "Shape", "Dispersion")))
head(param_summ)

### Create a plot comparing the true value, spline with uncertainty bounds, and the MLE estimate
p <- ggplot(param_summ, aes(x=year))  %>%
#	+ geom_ribbon(data = mle_plot, aes(ymin = lower_ci, ymax = upper_ci), fill="grey60", alpha=0.2) %>%
#	+ geom_line(data = mle_plot, aes(y=estimate), colour="#66c2a5") %>%
#	+ geom_ribbon(data = plot_df, aes(ymin = perc_50_lower, ymax = perc_95_upper), fill="grey70", alpha = 0.5)  %>%
	+ geom_line(aes(y = median, group = jdate), size=1) %>% 
#	+ geom_line(data = filter(param_summ, year == 1960), aes(y = median), colour="black", size=1) %>% 
	+ theme_classic() %>%
	+ scale_x_continuous(name = "Year", breaks=round(seq(1940,2020,by=10)), expand = c(0, 0)) %>%
	+ scale_y_continuous(name="Parameter Estimate") %>%
	+ coord_cartesian(xlim=c(1930,2020)) %>%
    + facet_grid(param ~jdate, scales="free_y") %>%
	+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

### Plot
p












































mean_est <- exp(extract(fit_2, "mean_param")$mean_param)
ya <- apply(mean_est, 2, median)

plot(true_param$jdate, true_param$mean, type="l")
lines(fitting_df$jdate, apply(mean_est, 2, median), col="red")
#lines(fitting_df$jdate, apply(mean_est, 2, quantile, 0.975), col="blue")
#lines(fitting_df$jdate, apply(mean_est, 2, quantile, 0.025), col="blue")
lines(fitting_df$jdate, mean_est[1,], col="grey")

scale_est <- exp(extract(fit_2, "scale_param")$scale_param)
ya <- apply(scale_est, 2, median)

plot(true_param$jdate, true_param$scale, type="l")
lines(fitting_df$jdate, apply(scale_est, 2, median), col="red")
lines(fitting_df$jdate, apply(scale_est, 2, quantile, 0.975), col="blue")
lines(fitting_df$jdate, apply(scale_est, 2, quantile, 0.025), col="blue")



plot(density(exp(extract(fit_2, "b_0_mean")$b_0_mean)))
abline(v = 5, col="red")

plot(density(exp(extract(fit_2, "b_0_scale")$b_0_scale)))
abline(v = 5/0.5, col="red")

plot(fit_2, plotfun = "trace", pars = "lambda_mean", inc_warmup = TRUE)
plot(fit_2, plotfun = "trace", pars = "lambda_scale", inc_warmup = TRUE)


plot(density(exp(extract(fit_2, "lambda_mean")$lambda_mean)))
lines(density(exp(rnorm(5000,log(22.8),5))), col="blue")
abline(v = mle_mean_gam$sp /  mle_mean_gam$smooth[[1]]$S.scale, col="red")
abline(v = init_vals[[1]]$lambda_mean, col="green")



plot(density(extract(fit_2, "lambda_mean")$lambda_mean))
lines(density(rgamma(5000,1.5, rate=0.06579)), col="blue")
abline(v = mle_mean_gam$sp /  mle_mean_gam$smooth[[1]]$S.scale, col="red")
abline(v = init_vals[[1]]$lambda_mean, col="green")


plot(density(extract(fit_2, "lambda_scale")$lambda_scale))
lines(density(rgamma(5000,10,0.002)), col="blue")
abline(v = mle_scale_gam$sp /  mle_scale_gam$smooth[[1]]$S.scale, col="red")





###########################################################################
###  Run with fixed lambda
###########################################################################
### Create the data to send to stan model
data_fitting <- list(N = length(fitting_df$precip), 
	basis_dim = basis_dim, 
	y=fitting_df$precip, 
	X = as.matrix(X_reparam),  
	S = as.matrix(mle_mean_gam$smooth[[1]]$S[[1]]), 
	b_0_mean_prior=b_0_mean_prior, 
	b_0_scale_prior=b_0_scale_prior, 
	lambda_mean=2000, 
	lambda_scale=4000)

str(data_fitting)

#init_vals <- list(list(), list(), list())
init_vals <- list(list(b_0_mean = b_0_mean_prior[1], 
	b_0_scale = b_0_scale_prior[1], 
	b_mean=b_mean_init, 
	b_scale=b_scale_init))
	
init_vals

### Fit the model
fit_2 <- stan(file = "./stan_models/03-seasonal_spline/e-gammals_fixedlambda.stan", 
				data = data_fitting, 
				init = init_vals,
				iter = 500, chains = 1, verbose = TRUE)













precision_mat <-   mle_mean_gam$smooth[[1]]$S[[1]] * (mle_mean_gam$sp[[1]])

precision_mat <-   mle_mean_gam$smooth[[1]]$S[[1]] * (mle_mean_gam$sp[[1]] /mle_mean_gam$smooth[[1]]$S.scale)

sigma <- solve(precision_mat)
yoop <- mvrnorm(n=1000, rep(0, 28), Sigma = sigma)

plot(yoop[,1], yoop[,2])

plot(density(yoop[,1]))
sd(yoop[,1])


gammals_fit


plot(coef(gammals_fit)[2:29], type="b")
lines(coef(gammals_fit)[31:60], col="red")

plot(coef(gammals_fit)[31:60], type="b")


coef(gammals_fit)
plot(coef(gammals_fit))

c(gammals_fit$sp)

gammals_fit$smooth[[1]]$S
gammals_fit$smooth[[2]]$S

spline_reparam[[1]]$S

gammals_fit$smooth[[1]]$S.scale
spline_reparam[[1]]$S.scale

require(MASS)

sigma <- gammals_fit$sp[[1]] * gammals_fit$smooth[[1]]$S.scale *  gammals_fit$smooth[[1]]$S[[1]]

yoop <- mvrnorm(n=1000, rep(0, 28), Sigma = 1/sigma)

 plot(yoop[,1], yoop[,2])


var(mvrnorm(n=1000, rep(0, 2), Sigma))
var(mvrnorm(n=1000, rep(0, 2), Sigma, empirical = TRUE))


gammals_fit$smooth[[1]]$S.scale * 

precision_mat <-   gammals_fit$smooth[[1]]$S[[1]] * gammals_fit$sp[[2]]
sigma <- solve(precision_mat)
yoop <- mvrnorm(n=1000, rep(0, 28), Sigma = sigma)

 plot(yoop[,1], yoop[,2])








###########################################################################
###  Calculate prior and initial values from MLE fit
###########################################################################

### Fit using mgcv directly
ctrl <- list(nthreads=4)

gammals_fit <- gam(list(precip  
		~ s(jdate, bs=c("cc"), k = c(n_knots_jdate) ) ,
		~ s(jdate, bs=c("cc"), k = c(n_knots_jdate) ) 
	),
	data=fitting_df, 
	knots = list(jdate=knots_jdate), 
	family=gammals, 
	select=FALSE, 
	method="REML",
	control= ctrl)

### Build the basis function again, either with or without the reparameterization
#spline_orig <- smoothCon(s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data=fitting_df, knots = list(jdate=knots_jdate), absorb.cons=FALSE, null.space.penalty = TRUE, scale.penalty=TRUE)
#spline_reparam <- smoothCon(s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data=fitting_df, knots = list(jdate=knots_jdate), absorb.cons=TRUE, null.space.penalty = FALSE, scale.penalty=TRUE)


###########################################################################
###  Calculate prior and initial values from MLE fit
###########################################################################

### Create the prior for the mean intercept
b_0_mean_prior <- c(summary(gammals_fit)$p.table[1,1], summary(gammals_fit)$p.table[1,2])

### Create a vector for intializing the mean
b_mean_init <- c(coef(gammals_fit)[2:c(n_knots_jdate-1)])
lambda_mean_init <- c(gammals_fit$sp)[[1]] / gammals_fit$smooth[[1]]$S.scale

### Create the prior for the scale intercept
b_0_scale_prior <- c(summary(gammals_fit)$p.table[2,1])
b_0_scale_prior[2] <- log(b_0_scale_prior[1] + summary(gammals_fit)$p.table[2,2]) - log(b_0_scale_prior[1])
b_0_scale_prior[1] <- log(b_0_scale_prior[1])

### Create a vector for intializing the scale
b_scale_init <- c(coef(gammals_fit)[c(n_knots_jdate+1):length(coef(gammals_fit))])
lambda_scale_init <- c(gammals_fit$sp)[[2]] / gammals_fit$smooth[[1]]$S.scale













jagam <- jagam(precip  ~ s(jdate, bs=c("cc"), k = c(n_knots_jdate)  ),
	data=fitting_df, 
	knots = list(jdate=knots_jdate), 
	family=gaussian(link = "log"), 
    file="test.jag")

jagam$pregam$S
jagam$jags.ini$lambda

precision_mat <- jagam$pregam$S[[1]] * jagam$jags.ini$lambda
sigma <- solve(precision_mat)

yoop <- mvrnorm(n=1000, rep(0, 28), Sigma = sigma)

 plot(yoop[,1], yoop[,2])


multi_normal_prec






### Build the basis function again, either with or without the reparameterization
spline_orig <- smoothCon(s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data=fitting_df, knots = list(jdate=knots_jdate), null.space.penalty = FALSE)
spline_reparam <- smoothCon(s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data=fitting_df, knots = list(jdate=knots_jdate), absorb.cons=TRUE, null.space.penalty = FALSE)

#spline_orig <- smoothCon(s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data=fitting_df, knots = list(jdate=knots_jdate), null.space.penalty = TRUE)
#spline_reparam <- smoothCon(s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data=fitting_df, knots = list(jdate=knots_jdate), absorb.cons=TRUE, null.space.penalty = TRUE)

### Extract the matrices for basis and penalty term
X_orig <- spline_orig[[1]]$X
s_orig <- spline_orig[[1]]$S
s_reparam <- spline_reparam[[1]]$S[[1]]

### Reparameterize both using the QR decomposition following Wood 
### Where Z is the Q matrix without the first column, used to reparameterize
C <- rep(1, nrow(X_orig)) %*% X_orig
qrc <- qr(t(C))
Z <- qr.Q(qrc,complete=TRUE)[,(nrow(C)+1):ncol(C)]

### Calculate reparameterized matrices for basis and penalty
X_reparam <- X_orig%*%Z
#X_reparam2 <- spline_reparam[[1]]$X

head(X_orig)
head(X_reparam)
head(s_orig)
head(s_reparam)
#b <- gam(log(rate_mle) ~ s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data=fitting_df, knots = list(jdate=knots_jdate), select=TRUE)

### Create the data to be used for model
basis_dim <- dim(X_reparam)[2]








### Fit the model
fit_2 <- stan(file = "./stan_models/03-seasonal_spline/e-gamma_hurdle_spline_loc_e.stan", 
				data = data_fitting, 
				init = init_vals,
				iter = 800, chains = 1, verbose = FALSE)
				
				
				
				
				
				
				
				
				
				
				




###########################################################################
###  Calculate prior and initial values from MLE fit
###########################################################################
### Fit a gam using the same basis 
mle_mean_gam <- gam(log(mean) ~ s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data=mle_fit, knots = list(jdate=knots_jdate), select=FALSE)
summary(mle_mean_gam)
plot(mle_mean_gam)

### Create the prior for the mean intercept
b_0_mean_prior <- c(summary(mle_mean_gam)$p.table[1], summary(mle_mean_gam)$p.table[2])

### Create a vector for intializing the mean
b_mean_init <- c(coef(mle_mean_gam)[2:length(coef(mle_mean_gam))])

lambda_mean_init <- c(mle_mean_gam$sp)

### Doublecheck
b_init_test <- X_reparam %*% b_mean_init +  b_0_mean_prior[1]
plot(fitting_df$jdate, exp(b_init_test), type="l")
rm(b_init_test)


### Fit a gam using the same basis 
mle_scale_gam <- gam(log(scale) ~ s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data=mle_fit, knots = list(jdate=knots_jdate), select=FALSE)
summary(mle_scale_gam)
plot(mle_scale_gam)

### Create the prior for the scale intercept
b_0_scale_prior <- c(summary(mle_scale_gam)$p.table[1], summary(mle_scale_gam)$p.table[2])

### Create a vector for intializing the scale
b_scale_init <- c(coef(mle_scale_gam)[2:length(coef(mle_scale_gam))])

lambda_scale_init <- c(mle_scale_gam$sp)

### Doublecheck
b_init_test <- X_reparam %*% b_scale_init +  b_0_scale_prior[1]
plot(fitting_df$jdate, exp(b_init_test), type="l")
rm(b_init_test)



###########################################################################
###  Run the model
###########################################################################

### Create the data to send to stan model
data_fitting <- list(N = length(fitting_df$precip), basis_dim = basis_dim, y=fitting_df$precip, X = X_reparam,  S = s_reparam, b_0_mean_prior=b_0_mean_prior, b_0_scale_prior=b_0_scale_prior, lambda_mean_init=lambda_mean_init, lambda_scale_init=lambda_scale_init)

str(data_fitting)

#init_vals <- list(list(), list(), list())
init_vals <- list(list(b_0_mean = b_0_mean_prior[1], b_0_scale = b_0_scale_prior[1], b_mean=b_mean_init, b_scale=b_scale_init, lambda_mean_init=lambda_mean_init, lambda_scale_init=lambda_scale_init))
init_vals

### Fit the model
fit_3 <- stan(file = "./stan_models/03-seasonal_spline/e-gamma_mean_scaleersion.stan", 
				data = data_fitting, 
				init = init_vals,
				iter = 200, chains = 1, verbose = FALSE)

### So fast





### Fit a gam using the same basis 
ya <- gammals(precip ~ s(jdate, bs=c("cc"), k = c(n_knots_jdate)),  ~ s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data=fitting_df, knots = list(jdate=knots_jdate), select=TRUE))
summary(ya)
plot(ya)


### Fit using mgcv directly
ctrl <- list(nthreads=4)

fit_3 <- gam(list(precip  
		~ s(jdate, bs=c("cc"), k = c(n_knots_jdate) ) ,
		~ s(jdate, bs=c("cc"), k = c(n_knots_jdate) ) 
	),
	data=fitting_df, 
	knots = list(jdate=knots_jdate), 
	family=gammals, 
	select=TRUE, 
	method="REML",
	control= ctrl)

### You can't put the scale.penalty option in here

c(coef(fit_3)[2:length(coef(fit_3))])

c(fit_3$sp)



summary(fit_3)
plot(fit_3)				
		
smoothCon		

maS <- norm(S) / norm(X, type = "I")^2				### Scaling factor for S

sapply(fit_3$smooth, "[[", "S.scale") / fit_3$sp

S_gammals <- fit_3$smooth[[1]]$S
S_gammals <- S_gammals[[1]]

x_mat <- model.matrix(fit_3)
X_gammals <- x_mat[,seq(2,29)]

maS <- norm(S_gammals) / norm(X_gammals, type = "I")^2


ya <- smoothCon(s(jdate, bs=c("cc"), k = c(n_knots_jdate)), data = fitting_df, knots = list(jdate=knots_jdate), absorb.cons=TRUE, scale.penalty=FALSE)
ya[[1]]$S


ya <- smoothCon(s(jdate, bs=c("cc"), k = c(n_knots)), data = fitting_df, knots = list(jdate=knots_jdate), absorb.cons=TRUE, scale.penalty=TRUE)
ya[[1]]$S	


	
###########################################################################
###  Check results
###########################################################################

print(fit_2, pars = c("b_0_mean", "b_0_scale"))

#
exp(1.59)
exp(2.28)
### True values are 5 and 5/0.5 or 10

### Check the trace plots to confirm the chains converge
plot(fit_2, plotfun = "trace", pars = "b_0_mean", inc_warmup = TRUE)

plot(density(exp(extract(fit_2, "b_0_mean")$b_0_mean)))
abline(v = 5, col="red")

plot(density(exp(extract(fit_2, "b_0_scale")$b_0_scale)))
abline(v = 5/0.5, col="red")

plot(fit_2, plotfun = "trace", pars = "lambda_mean", inc_warmup = TRUE)
plot(fit_2, plotfun = "trace", pars = "lambda_scale", inc_warmup = TRUE)

### Check chains without warmup
plot(fit_2, plotfun = "trace", pars = "b_mean", inc_warmup = TRUE)
plot(fit_2, plotfun = "trace", pars = "b_scale", inc_warmup = TRUE)

### Check the distributions of beta values
plot(fit_2, show_density = TRUE, ci_level = 0.5, pars = "b_mean", fill_color = "lightblue") + theme_classic()
plot(fit_2, show_density = TRUE, ci_level = 0.5, pars = "b_scale", fill_color = "lightblue") + theme_classic()

mean_est <- extract(fit_2, "mean_est")$mean_est
ya <- apply(mean_est, 2, median)

plot(true_param$jdate, true_param$mean, type="l")
lines(fitting_df$jdate, apply(mean_est, 2, median), col="red")
lines(fitting_df$jdate, apply(mean_est, 2, quantile, 0.975), col="blue")
lines(fitting_df$jdate, apply(mean_est, 2, quantile, 0.025), col="blue")


scale_est <- extract(fit_2, "scale_est")$scale_est
ya <- apply(scale_est, 2, median)

plot(true_param$jdate, true_param$scale, type="l")
lines(fitting_df$jdate, ya, col="red")


shape_est <- extract(fit_2, "shape_est")$shape_est
ya <- apply(shape_est, 2, median)

plot(true_param$jdate, true_param$shape, type="l")
lines(fitting_df$jdate, ya, col="red")

lambda_mean <- extract(fit_2, "lambda_mean")$lambda_mean
 plot(density(lambda_mean))
 
lambda_scale <- extract(fit_2, "lambda_scale")$lambda_scale
 plot(density(lambda_scale))

###########################################################################
###  Check results
###########################################################################
#print(fit_2)

### Create the full basis matrix by adding the intercept column
X_full_reparam <- cbind(rep(1,dim(X_reparam)[1]), X_reparam)
demo_basis_reparam <- cbind(rep(1,dim(demo_basis_reparam)[1]), demo_basis_reparam)

### Extract the spline coefficients and intercept for mean
b_mean <- extract(fit_2, "b_mean")$b_mean
b_0_mean <- extract(fit_2, "b_0_mean")$b_0_mean
### Combine the intercept and spline coefficients into a single matrix
b_full_mean <- cbind(matrix(b_0_mean, dim(b_mean)[1], 1), b_mean)

### Extract the spline coefficients and intercept for rate
b_scale <- extract(fit_2, "b_scale")$b_scale
b_0_scale <- extract(fit_2, "b_0_scale")$b_0_scale
### Combine the intercept and spline coefficients into a single matrix
b_full_scale <- cbind(matrix(b_0_scale, dim(b_scale)[1], 1), b_scale)

### Calculate the estimate of rate based on the jdate 1 to 365 dataframe
mean_est_jdate <- exp(demo_basis_reparam %*% t(b_full_mean))
mean_est <- data.frame(jdate_demo, rate=mean_est_jdate)

### Gather the results into a long format
mean_est_long <- mean_est %>%
	gather("draw", "mean", -jdate)

### Quick plot of all draws
ggplot(mean_est_long, aes(x=jdate, y=mean, group=draw)) + geom_line(colour="grey40", alpha=0.2) + theme_classic()

### Calculate some summary statistics (quantiles) to plot the spline results
mean_est_summary <- mean_est_long %>%
	group_by(jdate) %>%
	summarise(median = median(mean), perc_95_lower = quantile(mean, 0.025), perc_95_upper = quantile(mean, 0.975), perc_50_upper = quantile(mean, 0.25), perc_50_lower = quantile(mean, 0.75))

### Create a plot comparing the true value, spline with uncertainty bounds, and the MLE estimate
p <- ggplot(mean_est_summary, aes(x=jdate))  %>%
	+ geom_line(data = mle_fit, aes(y=mean), colour="#66c2a5") %>%
	+ geom_ribbon(data = mean_est_summary, aes(ymin = perc_95_lower, ymax = perc_95_upper), fill="grey70", alpha = 0.5)  %>%
	+ geom_ribbon(data = mean_est_summary, aes(ymin = perc_50_lower, ymax = perc_50_upper), fill="grey50", alpha = 0.5) %>% 
	+ geom_line(data = true_param, aes(y=mean), size=1) %>%
	+  geom_line(aes(y = median), colour="#fc8d62", size=1) %>% 
	+  theme_classic() %>%
	+ scale_x_continuous(name = "Julian Date", breaks=round(seq(1,365,length.out=20)), expand = c(0, 0)) %>%
	+ scale_y_continuous(name="Mean Parameter") %>%
	+ coord_cartesian(xlim=c(1,365))
### Plot
p

### Save the plot
ggsave("../output/model2_rate_spline_v_mle.png", p, width =12.5, height = 8, dpi = 300)


#### Repeat all of this for shape parameter
### Calculate the estimate of shape based on the jdate 1 to 365 dataframe
shape_est_jdate <- exp(demo_basis_reparam %*% t(b_full_shape))
shape_est <- data.frame(jdate_demo, rate=shape_est_jdate)

### Gather the results into a long format
shape_est_long <- shape_est %>%
	gather("draw", "shape", -jdate)

### Calculate summary statistics
shape_est_summary <- shape_est_long %>%
	group_by(jdate) %>%
	summarise(median = median(shape), perc_95_lower = quantile(shape, 0.025), perc_95_upper = quantile(shape, 0.975), perc_50_upper = quantile(shape, 0.25), perc_50_lower = quantile(shape, 0.75))

### Quick plot of all draws
ggplot(shape_est_long, aes(x=jdate, y=shape, group=draw)) + geom_line(colour="grey40", alpha=0.2) + theme_classic()

### Create a plot comparing the true value, spline with uncertainty bounds, and the MLE estimate
p <- ggplot(shape_est_summary, aes(x=jdate))  %>%
	+ geom_line(data = mle_fit, aes(y=shape), colour="#66c2a5") %>%
	+ geom_ribbon(data = shape_est_summary, aes(ymin = perc_95_lower, ymax = perc_95_upper), fill="grey70", alpha = 0.5)  %>%
	+ geom_ribbon(data = shape_est_summary, aes(ymin = perc_50_lower, ymax = perc_50_upper), fill="grey50", alpha = 0.5) %>% 
	+ geom_line(data = true_param, aes(y=shape), size=1) %>%
	+  geom_line(aes(y = median), colour="#fc8d62", size=1) %>% 
	+  theme_classic() %>%
	+ scale_x_continuous(name = "Julian Date", breaks=round(seq(1,365,length.out=20)), expand = c(0, 0)) %>%
	+ scale_y_continuous(name="Shape Parameter") %>%
	+ coord_cartesian(xlim=c(1,365))
### Plot
p
### Save the plot
ggsave("../output/model2_shape_spline_v_mle.png", p, width =12.5, height = 8, dpi = 300)



###########################################################################
###  Convert coeficients back to original basis (before reparam)
###########################################################################
### Need to go row by row of betas
b_mean_orig <- Z %*% b_mean[1,]

plot(jdate_demo$jdate, exp(demo_basis_x %*% b_mean_orig + b_0_mean[1]), type="l"); for(j in seq(1,n_knots-1)){lines(jdate_demo$jdate, exp(demo_basis_x[,j] * b_mean_orig[j] + b_0_mean[1]), col= rainbow(n_knots)[j])}

 

