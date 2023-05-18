################
# Calc env variability metrics
#
# NOTES:
# - determine how best to chunk up world for metric calculation
# - write results back to database
# - create tif files of all env metrics
# - worth thinking about the most appropriate lags for temporal autocorrelation (in terms of predictability)
# - can we fit regression models with t-distributed errors (estimating and not specificy degrees of freedom parameter) in R at a reasonable speed? (seems possibly -> https://stats.stackexchange.com/questions/117980/regression-with-t-distributed-errors-and-massrlm)
################


# specify dir -------------------------------------------------------------

#path for data on CY machine - remember trailing slash
temp_data_dir <- '~/Desktop/temp/'
precip_data_dir <- '~/Desktop/precip/'
#directory to save out intermediate file
sample_output_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/Sample_output/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)


# read in data -------------------------------------------------

#which months averaged over
months <- 1
temp_out <- read.csv(paste0(temp_data_dir, 'Temp-', 
                           paste0(as.character(months), collapse = '_'), '.csv'))


# Env metrics -------------------------------------------------------------

#metrics:
#-slope over time
#-variance (or sd) of residuals (i.e., variability of detrended time series)
#-kurtosis (i.e., heavy-tailed ness - or model time series using t-distribution and estimate degrees of freedom parameter)
#-temporal autocorrelation (of residuals - see Leung et al. 2020 Eco Letters for use of this metric for 'predictability')

#North America only
temp_out2 <- dplyr::filter(temp_out, 
                           lon > -170, lon < -50,
                           lat > 15, lat < 75)
#Entire globe
#temp_out2 <- temp_out

#add unique cell id
#convert to data table (much faster)
temp_out3 <- dplyr::group_by(temp_out2, lat, lon) %>%
  dplyr::mutate(cell_id = cur_group_id()) %>%
  dplyr::ungroup() %>%
  data.table::as.data.table()

#unique cells
uci <- unique(temp_out3$cell_id)

#data.frame to fill
env_out <- data.frame(cell_id = uci,
                      lon = NA,
                      lat = NA,
                      slope = NA, #slope of linear model
                      se_slope = NA, #standard error of slope
                      sd_resid = NA, #sd of residuals for linear model
                      kurt = NA, #kurtosis
                      skew = NA, #skewness
                      rho_l1 = NA, #temporal autocorrelation lag 1
                      rho_l2 = NA, #lag 2
                      rho_l3 = NA,
                      rho_l4 = NA,
                      rho_l5 = NA)

#loop through each cell to calc metrics
#takes a handful (on the order of 30 min) of minutes on laptop
for (i in 1:length(uci))
{
  #i <- 1
  print(paste0('Processing cell ', i, ' of ', length(uci)))
  
  te <- temp_out3[cell_id == uci[i],]
  
  #linear model fit
  fit <- summary(lm(temp ~ year, data = te))
  
  #residuals from model
  t_resid <- residuals(fit)
  
  #temporal autocorrelation
  ar_fit <- acf(t_resid, lag.max = 5, plot = FALSE)
  
  #fill df
  env_out$lon[i] <- te$lon[1]
  env_out$lat[i] <- te$lat[2]
  env_out$slope[i] <- fit$coefficients[2,1]
  env_out$se_slope[i] <- fit$coefficients[2,2]
  env_out$sd_resid[i] <- sd(t_resid)
  env_out$kurt[i] <- moments::kurtosis(t_resid)
  env_out$skew[i] <- moments::skewness(t_resid)
  env_out$rho_l1[i] <- ar_fit$acf[2]
  env_out$rho_l2[i] <- ar_fit$acf[3]
  env_out$rho_l3[i] <- ar_fit$acf[4]
  env_out$rho_l4[i] <- ar_fit$acf[5]
  env_out$rho_l5[i] <- ar_fit$acf[6]
}

#save out intermediate object to avoid having to reprocess
saveRDS(env_out, paste0(sample_output_dir, 'env_var_out.rds'))
#env_out <- readRDS(paste0(sample_output_dir, 'env_var_out.rds'))


# explore -----------------------------------------------------------------

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
