# TITLE:            Configuration File
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATE:             March 2024 
# OVERVIEW:         File to specify environmental variables used across multiple scripts (e.g., dir). This file will be sourced in each code file to avoid the necessity of repeatedly specifying directories, dates, etc. 

# Specify data directory
dir='/mnt/research/ibeem/variability_testing/'
# dir='/mnt/research/ibeem/variability/'

# Specify model run date (for step 6 and beyond)
gl_run_date <- '2024-04-14'
# gl_run_date <- '2023-10-17'


# # Download necessary packages if not already downloaded
# packages <- c("tidyverse", "sf", "data.table", "moments", "terra", "ape",
#               "phytools", "MCMCvis", "ape", "geiger", "ggplot2",
#               "rnaturalearth", "viridis", "plotfunctions")
# 
# package_download <- function(x){
#   if ((x %in% installed.packages()[,1]) == FALSE){
#     print(paste0("Installing package: ", x))
#     install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
#   }
# }
# 
# lapply(packages, package_download)
# 
# # Separately check cmdstanr (not available on CRAN)
# if (("cmdstanr" %in% installed.packages()[,1]) == FALSE){
#   print("Installing package: cmdstanr")
#   install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#   cmdstanr::install_cmdstan(cores = 2)
# }