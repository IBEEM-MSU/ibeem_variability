# TITLE:            Configuration File
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATE:             March 2024 
# OVERVIEW:         File to specify environmental variables used across multiple scripts (e.g., dir). This file will be sourced in each code file to avoid the necessity of repeatedly specifying directories, dates, etc. 

# Specify data directory 
dir='/mnt/research/ibeem/variability_testing/'

# Download necessary packages if not already downloaded 
packages <- c("tidyverse", "sf", "data.table", "moments", "terra", "ape", 
              "phytools", "cmdstanr", "MCMCvis", "ape", "geiger", "ggplot2", 
              "rnaturalearth", "colorplaner", "viridis", "plotfunctions")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if ((x %in% installed.packages()[,1]) == FALSE){
      install.packages(x, dependencies = TRUE)
    }
  }
)