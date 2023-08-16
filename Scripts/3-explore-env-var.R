################
# Explore env variability metrics
#
# NOTES:
# - apply land mask
################


# Specify dir --------------------------------------------------

#path CY machine
dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)


# read in data -------------------------------------------------

#environmental variability
env_var <- read.csv(paste0(dir, 'Data/L2/climate/era5/Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv'))
env_var_GAM <- read.csv(paste0(dir, 'Data/L2/climate/era5/Env-var-GAM-1_2_3_4_5_6_7_8_9_10_11_12.csv'))
#READ IN FROM L1
env_season <- read.csv(paste0(dir, 'Data/L1/climate/era5/Env-seasonality-1_2_3_4_5_6_7_8_9_10_11_12.csv'))

#landmask
lm <- sf::st_read(paste0(dir, 'Data/L0/landmask/ne_10m_land.shp')) %>%
  dplyr::filter(featurecla == 'Land')


# rasterize and mask env var data -------------------------------------------------------

#multiband raster (each band different env metric)
#var = temp or precip
#function
mbr_fun <- function(input, VAR)
{
  #rasterize
  trast <- dplyr::filter(input, var == VAR) %>%
    dplyr::select(-cell_id, -var) %>%
    terra::rast(crs = "epsg:4326") %>%
    #mask out water
    terra::mask(terra::vect(lm))
  
  #get relative slope
  sl_resid <- trast[['slope']] / trast[['sd_resid']]
  names(sl_resid) <- 'rel_slope'
  
  #add band
  trast2 <- c(trast, sl_resid)
  
  return(trast2)
}


#run function - putting some cells outside of lat/lon bounds for some reason
ev_temp <- mbr_fun(input = env_var, VAR = 'temp')
ev_precip <- mbr_fun(input = env_var, VAR = 'precip')
ev_temp_GAM <- mbr_fun(input = env_var_GAM, VAR = 'temp')
ev_precip_GAM <- mbr_fun(input = env_var_GAM, VAR = 'precip')

#bands
names(ev_temp)

#check
# plot(ev_temp[[1]])


# stats -------------------------------------------------------------------

#median over globe
terra::global(ev_temp[['slope']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['sd_resid']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['kurt']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['skew']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['spectral_beta']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['rho_l1']], fun = function(x) median(x, na.rm = TRUE))
terra::global(ev_temp[['rel_slope']], fun = function(x) median(x, na.rm = TRUE))


# plots -----------------------------------------------------------------

# NA_out <- dplyr::filter(env_var, 
#                            lon > -170, lon < -50,
#                            lat > 15, lat < 75)

#function to plot lm and gam
gc_fun <- function(rast, rast_gam, var)
{
  par(mfrow = c(2,1))
  
  #get range of values  
  rng <- range(terra::global(c(rast[[var]], rast_gam[[var]]), 
                             fun = function(x) range(x, na.rm = TRUE)))
  plot(rast[[var]], 
       main = var,
       range = rng)
  plot(rast_gam[[var]], 
       main = paste0(var, ' - GAM'),
       range = rng)
}

#run fun
#slope doesn't vary since calcated with lm for both
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'slope')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'sd_resid')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'kurt')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'skew')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'spectral_beta')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'rho_l1')
gc_fun(rast = ev_temp, 
       rast_gam = ev_temp_GAM,
       var = 'rel_slope')


#which areas are highly predictable (P; high temporal autocorrelation) on short time scales (S) and have low intrinsic variability (IV)?
#high P, short S, low IV = short LH
#low P, long S, high IV = long LH

#high P/short S = faster LH
#high IV = slower LH

#is P at lag = gen time the same across all species?

