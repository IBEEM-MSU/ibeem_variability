################
# Explore env variability metrics
#
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
env_main <- read.csv(paste0(dir, 'Data/L2/climate/era5/Env-main.csv'))


# rasterize ------------------------------------------------------------------

#multiband raster (each band different env metric)
#var = temp or precip
#function
mbr_fun <- function(input, VAR)
{
  #rasterize
  trast <- dplyr::filter(input, var == VAR, valid == TRUE) %>%
    dplyr::select(lon, lat, mean, slope, sd_resid, sd_season, kurt, skew, 
                  spectral_beta, rho_l1, rel_slope, PC1, PC2, PC3) %>%
    terra::rast(crs = "epsg:4326")
  return(trast)
}


#run function - putting some cells outside of lat/lon bounds for some reason
ev_temp <- mbr_fun(input = env_main, VAR = 'temp')
ev_precip <- mbr_fun(input = env_main, VAR = 'precip')

#bands
names(ev_precip)


# plots -------------------------------------------------------------------

#temperature
plot(ev_temp[[1]], main = 'Mean')
plot(ev_temp[[3]], main = 'Inter-annual sd')
plot(ev_temp[[4]], main = 'Intra-annual sd')

plot(ev_temp[[5]], main = 'Kurtosis')
plot(ev_temp[[6]], main = 'Skew')
plot(ev_temp[[7]], main = 'Spectral exponent')
plot(ev_temp[[8]], main = 'Rho Lag 1')

plot(ev_temp[[2]], main = 'Slope')
plot(ev_temp[[9]], main = 'Relative slope')

plot(ev_temp[[10]], main = 'PC1 (+ = v mean, ^ inter, ^ intra)')
plot(ev_temp[[11]], main = 'PC2 (+ = ^ inter, v intra)')
plot(ev_temp[[12]], main = 'PC3 (+ = ^ mean, ^ intra)')


#Precip
plot(ev_precip[[1]], main = 'Mean')
plot(ev_precip[[3]], main = 'Inter-annual sd')
plot(ev_precip[[4]], main = 'Intra-annual sd')

plot(ev_precip[[5]], main = 'Kurtosis')
plot(ev_precip[[6]], main = 'Skew')
plot(ev_precip[[7]], main = 'Spectral exponent')
plot(ev_precip[[8]], main = 'Rho Lag 1')

plot(ev_precip[[2]], main = 'Slope')
plot(ev_precip[[9]], main = 'Relative slope')

#there are weird - is this bc precip data is less reliable?
plot(ev_precip[[10]], main = 'PC1 (+ = v mean, ^ inter, ^ intra)')
plot(ev_precip[[11]], main = 'PC2 (+ = ^ inter, v intra)')
plot(ev_precip[[12]], main = 'PC3 (+ = ^ mean, ^ intra)')


# stats -------------------------------------------------------------------

# #median over globe
# terra::global(ev_temp[['slope']], fun = function(x) median(x, na.rm = TRUE))
# terra::global(ev_temp[['sd_resid']], fun = function(x) median(x, na.rm = TRUE))
# terra::global(ev_temp[['kurt']], fun = function(x) median(x, na.rm = TRUE))
# terra::global(ev_temp[['skew']], fun = function(x) median(x, na.rm = TRUE))
# terra::global(ev_temp[['spectral_beta']], fun = function(x) median(x, na.rm = TRUE))
# terra::global(ev_temp[['rho_l1']], fun = function(x) median(x, na.rm = TRUE))
# terra::global(ev_temp[['rel_slope']], fun = function(x) median(x, na.rm = TRUE))


# compare lm and GAM -----------------------------------------------------------------

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

