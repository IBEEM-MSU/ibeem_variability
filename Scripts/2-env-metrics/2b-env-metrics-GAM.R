################
# Calc env variability metrics
#
# args:
# - in file
# - out file
#
# NOTES:
# - what are the most appropriate lags for temporal autocorrelation (in terms of predictability)
# - can we fit regression models with t-distributed errors (estimating and not specificy degrees of freedom parameter) in R at a reasonable speed? (seems possibly -> https://stats.stackexchange.com/questions/117980/regression-with-t-distributed-errors-and-massrlm)
################

# get args fed to script --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args[1] = in file (L1 env data) - '/mnt/research/ibeem/L1/climate/era5/ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv'
#args[2] = out file (L2 env data) - '/mnt/research/ibeem/L2/climate/era5/Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv'


# load packages -----------------------------------------------------------

library(tidyverse)
library(data.table)
library(moments)


# read in data -------------------------------------------------

env_out <- read.csv(args[1])


# Env metrics -------------------------------------------------------------

#metrics:
#-slope over time (and se)
#-variance (or sd) of residuals (i.e., variability of detrended time series)
#-kurtosis (i.e., heavy-tailed ness - or model time series using t-distribution and estimate degrees of freedom parameter)
#-skew (i.e., skew of residuals)
#-temporal autocorrelation (of residuals - see Leung et al. 2020 Eco Letters for use of this metric for 'predictability')

#add unique cell id
#convert to data table (much faster)
env_out2 <- dplyr::group_by(env_out, lat, lon) %>%
  dplyr::mutate(cell_id = cur_group_id()) %>%
  dplyr::ungroup() %>%
  data.table::as.data.table()

#unique cells
uci <- unique(env_out2$cell_id)

#data.frame to fill
env_df <- data.frame(cell_id = rep(uci, each = 2),
                     var = rep(c('temp', 'precip'), length(uci)),
                     lon = NA,
                     lat = NA,
                     mean = NA,
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
time <- proc.time()
counter <- 1
for (i in 1:length(uci))
{
  #i <- 1
  print(paste0('Processing cell ', i, ' of ', length(uci)))
  
  te <- env_out2[cell_id == uci[i],]
  
  #linear model fit
  fit_temp <- summary(lm(temp ~ year, data = te))
  fit_precip <- summary(lm(precip ~ year, data = te))
  
  #residuals from model
  temp_resid <- residuals(fit_temp)
  precip_resid <- residuals(fit_precip)
  
  #temporal autocorrelation
  temp_acf <- acf(temp_resid, lag.max = 5, plot = FALSE)
  precip_acf <- acf(precip_resid, lag.max = 5, plot = FALSE)
  
  #fill df
  env_df$lon[counter:(counter + 1)] <- rep(te$lon[1], 2)
  env_df$lat[counter:(counter + 1)] <- rep(te$lat[1], 2)
  env_df$mean[counter:(counter + 1)] <- c(mean(te$temp), mean(te$precip))
  env_df$slope[counter:(counter + 1)] <- c(fit_temp$coefficients[2,1], fit_precip$coefficients[2,1])
  env_df$se_slope[counter:(counter + 1)] <- c(fit_temp$coefficients[2,2], fit_precip$coefficients[2,2])
  env_df$sd_resid[counter:(counter + 1)] <- c(sd(temp_resid), sd(precip_resid))
  env_df$kurt[counter:(counter + 1)] <- c(moments::kurtosis(temp_resid), moments::kurtosis(precip_resid))
  env_df$skew[counter:(counter + 1)] <- c(moments::skewness(temp_resid), moments::skewness(precip_resid))
  env_df$rho_l1[counter:(counter + 1)] <- c(temp_acf$acf[2], precip_acf$acf[2])
  env_df$rho_l2[counter:(counter + 1)] <- c(temp_acf$acf[3], precip_acf$acf[3])
  env_df$rho_l3[counter:(counter + 1)] <- c(temp_acf$acf[4], precip_acf$acf[4])
  env_df$rho_l4[counter:(counter + 1)] <- c(temp_acf$acf[5], precip_acf$acf[5])
  env_df$rho_l5[counter:(counter + 1)] <- c(temp_acf$acf[6], precip_acf$acf[6])
  
  counter <- counter + 2
}
print('Run time')
print(proc.time() - time)


# write to csv ------------------------------------------------------------

write.csv(env_df, args[2], row.names = FALSE)
