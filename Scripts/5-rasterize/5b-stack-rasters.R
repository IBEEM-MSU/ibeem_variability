#############################
# Stack rasters
#
#############################


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)


# Specify top-level directory -------------------------------------------------------

dir <- '/mnt/research/ibeem/variability/'
# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# read in bird data -------------------------------------------------------

# #read in life history and other traits
# #use to filter species of interest (e.g., no seabirds)
# bird_df <- read.csv(paste0(dir, 'data/L2/main-bird-data.csv')) %>%
#   dplyr::filter(!is.na(GenLength))


# read in rasters ---------------------------------------------------------

#list files of interest
lf <- list.files(paste0(dir, 'data/L2/range-raster/'), full.names = TRUE)
lf_gl <- grep('GenLength', lf, value = TRUE)
lf_dh <- grep('delta_haldane', lf, value = TRUE)

#read in stacked
gl_stack <- terra::rast(lf_gl)
dh_stack <- terra::rast(lf_dh)

#calculate median
med_gl <- terra::app(gl_stack, fun = function(x) median(x, na.rm = TRUE))
med_dh <- terra::app(dh_stack, fun = function(x) median(x, na.rm = TRUE))

#calculate number of species in each grid cell
n_sp <- terra::app(gl_stack, fun = function(x) sum(!is.na(x)))

#plot
plot(med_gl)
plot(n_sp)

