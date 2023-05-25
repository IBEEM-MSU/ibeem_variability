################
# Explore env variability metrics
#
# NOTES:
# - apply land mask
################


# Specify dir --------------------------------------------------

dir <- '~/Downloads/era5/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)


# read in data -------------------------------------------------

env_var <- read.csv(paste0(dir, 'Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv'))
env_GAM_var <- read.csv(paste0(dir, 'Env-var-GAM-1_2_3_4_5_6_7_8_9_10_11_12.csv'))


# stats -------------------------------------------------------------------

#degrees C per year
dplyr::filter(env_var, var == 'temp') %>%
  summarize(median(slope))
dplyr::filter(env_var, var == 'temp') %>%
  summarize(median(kurt))
dplyr::filter(env_var, var == 'temp') %>%
  summarize(median(skew))
dplyr::filter(env_var, var == 'temp') %>%
  summarize(median(spectral_beta))
dplyr::filter(env_var, var == 'temp') %>%
  summarize(median(rho_l1))


# plots -----------------------------------------------------------------

# NA_out <- dplyr::filter(env_var, 
#                            lon > -170, lon < -50,
#                            lat > 15, lat < 75)

#convert to raster to visualize
dplyr::filter(env_var, var == 'temp') %>%
  dplyr::select(lon, lat, slope) %>%
  terra::rast() %>%
  plot(main = 'TEMP - slope')

dplyr::filter(env_var, var == 'temp') %>%
  dplyr::select(lon, lat, sd_resid) %>%
  terra::rast() %>%
  plot(main = 'TEMP - sd resid')

dplyr::filter(env_var, var == 'temp') %>%
  dplyr::mutate(sl_sd = slope / sd_resid) %>%
  dplyr::select(lon, lat, sl_sd) %>%
  terra::rast() %>%
  plot(main = 'TEMP - slope/sd')

dplyr::filter(env_var, var == 'temp') %>%
  dplyr::select(lon, lat, kurt) %>%
  terra::rast() %>%
  plot(main = 'TEMP - kurt')

dplyr::filter(env_var, var == 'temp') %>%
  dplyr::select(lon, lat, skew) %>%
  terra::rast() %>%
  plot(main = 'TEMP - skew')

#need to figure out why there are some negative values...
dplyr::filter(env_var, var == 'temp') %>%
  dplyr::select(lon, lat, spectral_beta) %>%
  terra::rast() %>%
  plot(main = 'TEMP - spectral exponent')

dplyr::filter(env_var, var == 'temp') %>%
  dplyr::mutate(rho2_l1 = rho_l1^2) %>%
  dplyr::select(lon, lat, rho2_l1) %>%
  terra::rast() %>%
  plot(main = 'TEMP - rho^2_l1')

#which areas are highly predictable (P; high temporal autocorrelation) on short time scales (S) and have low intrinsic variability (IV)?
#high P, short S, low IV = short LH
#low P, long S, high IV = long LH

#high P/short S = faster LH
#high IV = slower LH

#is P at lag = gen time the same across all species?

