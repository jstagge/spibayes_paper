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
###  Install and Load Devtools for custom packages
###########################################################################
install.packages("fs", repos = "https://cran.case.edu/")
install.packages("usethis", repos = "https://cran.case.edu/")
install.packages("vctrs", repos = "https://cran.case.edu/")
install.packages("lifecycle", repos = "https://cran.case.edu/")

install.packages("backports", repos = "https://cran.case.edu/")
install.packages("rprojroot", repos = "https://cran.case.edu/")
install.packages("devtools", repos = "https://cran.case.edu/")
install.packages("broom", repos = "https://cran.case.edu/")
install.packages("ggplot2", repos = "https://cran.case.edu/")
install.packages("Rcpp", repos = "https://cran.case.edu/")
#require(devtools)

###########################################################################
###  Install packages from R repo
###########################################################################

install.packages("dplyr", repos = "https://cran.case.edu/")
install.packages("tibble", repos = "https://cran.case.edu/")
install.packages("tidyverse", repos = "https://cran.case.edu/")
install.packages("here", repos = "https://cran.case.edu/")

### To access GHCND
#install.packages("isdparser", repos = "https://cloud.r-project.org")
#install.packages("rnoaa", repos = "https://cloud.r-project.org")

### To save in SVG
install.packages("svglite", repos = "https://cloud.r-project.org")
install.packages("viridis", repos = "https://cloud.r-project.org")

### Packages for spi
install.packages("fitdistrplus", repos = "https://cloud.r-project.org")
install.packages("lubridate", repos = "https://cloud.r-project.org")

install.packages("mgcv", repos = "https://cloud.r-project.org")
### There was an error message with this version
#install_version("mgcv", version = "1.8-31", repos = "http://cran.us.r-project.org")
#install.packages("dplyr",  repos = "https://cloud.r-project.org")

install.packages("ggthemes", repos = "https://cloud.r-project.org")
install.packages("R6", repos = "https://cloud.r-project.org")

#install.packages("rgdal", repos = "https://cloud.r-project.org")
#install.packages("ggrepel", repos = "https://cloud.r-project.org")
#install.packages("sf", repos = "https://cloud.r-project.org")
#install.packages("maps", repos = "https://cloud.r-project.org")
install.packages("furrr", repos = "https://cloud.r-project.org")
install.packages("zoo", repos = "https://cloud.r-project.org")
#install.packages("rnaturalearth", repos = "https://cloud.r-project.org")

install.packages("vroom", repos = "https://cloud.r-project.org")
install.packages("cmdstanr", repos = "https://mc-stan.org/r-packages")

require(cmdstanr)
install_cmdstan()
cmdstan_path()
cmdstan_version()

install.packages("rstan", repos = "https://cloud.r-project.org")
install.packages("StanHeaders", repos = "https://cloud.r-project.org")

###########################################################################
###  Install custom spibayes package
###########################################################################
### For custom MLE functions
#devtools::install_github("staggelab/spibayes",  ref= "add_variance")

install.packages("/users/PAS1627/stagge11/code/spibayes", 
repos = NULL, 
type = "source")

