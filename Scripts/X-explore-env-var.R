################
# Explore env variability metrics
#
# NOTES:
# - visualize differences across world
################


# Specify dir --------------------------------------------------

XXXX


# load packages -----------------------------------------------------------

library(tidyverse)
library(data.table)
library(moments)


# read in data -------------------------------------------------

env_out <- read.csv('/mnt/research/ibeem/L2/climate/era5/Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv')
# env_out <- read.csv('/mnt/research/ibeem/L2/climate/era5/Env-var-6_7_8.csv')


# explore -----------------------------------------------------------------

# temp_out2 <- dplyr::filter(env_out, 
#                            lon > -170, lon < -50,
#                            lat > 15, lat < 75)


#convert to raster to visualize
dplyr::select(env_out, lon, lat, slope) %>%
  terra::rast() %>%
  plot(main = 'slope')

dplyr::select(env_out, lon, lat, sd_resid) %>%
  terra::rast() %>%
  plot(main = 'sd resid')

dplyr::mutate(env_out, sl_sd = slope / sd_resid) %>%
  dplyr::select(lon, lat, sl_sd) %>%
  terra::rast() %>%
  plot(main = 'slope/sd')


# plot --------------------------------------------------------------------


