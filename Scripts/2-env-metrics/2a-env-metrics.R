################
# Calc env variability metrics
#
# args:
# - in file
# - out file
#
################

# get args fed to script --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
# args <- rep(NA, 2)
# args[1] <- '/mnt/research/ibeem/variability/data/L1/climate/era5/ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv'
# args[2] <- '/mnt/research/ibeem/variability/data/L2/climate/era5/Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv'


# load packages -----------------------------------------------------------

library(tidyverse)
library(data.table)
library(moments)
library(lomb)


# read in data -------------------------------------------------

env_out <- read.csv(args[1])


# Env metrics -------------------------------------------------------------

#metrics:
#-mean
#-slope over time (and se)
#-variance (or sd) of residuals (i.e., variability of detrended time series)
#-kurtosis (i.e., heavy-tailed ness - or model time series using t-distribution and estimate degrees of freedom parameter)
#-skew (i.e., skew of residuals)
#-spectral exponent ('color' of noise)
#-temporal autocorrelation (of residuals - see Leung et al. 2020 Eco Letters for use of this metric for 'predictability')

#add unique cell id
#convert to data table (much faster)
env_out2 <- dplyr::group_by(env_out, lat, lon) %>%
  dplyr::mutate(cell_id = cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cell_id, year) %>%
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
                     sd_year = NA, #sd of residuals for linear model
                     sd_season = NA, #mean yearly sd across months
                     kurt = NA, #kurtosis
                     skew = NA, #skewness
                     sp_color_year = NA, #spectral exponent
                     rho_l1 = NA, #temporal autocorrelation lag 1
                     rho_l2 = NA, #lag 2
                     rho_l3 = NA,
                     rho_l4 = NA,
                     rho_l5 = NA)

#loop through each cell to calc metrics
time <- proc.time()
counter <- 1
for (i in 1:length(uci))
{
  #i <- 1
  print(paste0('Processing cell ', i, ' of ', length(uci)))
  
  te <- env_out2[cell_id == uci[i],]
  
  #linear model fit (yearly env ~ year)
  fit_temp <- summary(lm(mean_temp ~ year, data = te))
  fit_precip <- summary(lm(mean_precip ~ year, data = te))
  
  #residuals from model
  temp_resid <- residuals(fit_temp)
  precip_resid <- residuals(fit_precip)
  
  #spectral analysis using Lomb-Scargle Periodogram
  #between freq 2/(n*dt) and 1/(2*dt), where dt = 1 and n = 72; 0.0278 to 0.5
  #following Marshall and Burgess 2015 Eco Letters
  #see also Vasseur and Yodzis 2004 Ecology
  #period = 1/freq; ~2 - 36 years
  ll_freq <- 0.0278
  ul_freq <- 0.5
  temp_spec <- lomb::lsp(temp_resid, 
                         from = ll_freq, to = ul_freq, 
                         type = 'frequency',
                         normalize =  'standard', 
                         plot = FALSE)
  #spectral exponent (1/f^beta)
  temp_spec_fit <- summary(lm(log10(temp_spec$power) ~ log10(temp_spec$scanned)))$coefficients[,1]
  #precip - only run if residuals exist (i.e., all precip values were not 0)
  if (sd(precip_resid) > 0)
  {
    precip_spec <- lomb::lsp(precip_resid, 
                             from = ll_freq, to = ul_freq, 
                             type = 'frequency', 
                             normalize = 'standard', 
                             plot = FALSE)
    precip_spec_fit <- summary(lm(log10(precip_spec$power) ~ log10(precip_spec$scanned)))$coefficients[,1]
  } else {
    precip_spec_fit <- rep(NA, 2)
  }
  
  #temporal autocorrelation
  temp_acf <- acf(temp_resid, lag.max = 5, plot = FALSE)
  precip_acf <- acf(precip_resid, lag.max = 5, plot = FALSE)
  
  #fill df
  env_df$lon[counter:(counter + 1)] <- rep(te$lon[1], 2)
  env_df$lat[counter:(counter + 1)] <- rep(te$lat[1], 2)
  env_df$mean[counter:(counter + 1)] <- c(mean(te$mean_temp), mean(te$mean_precip))
  env_df$slope[counter:(counter + 1)] <- c(fit_temp$coefficients[2,1], 
                                           fit_precip$coefficients[2,1])
  env_df$se_slope[counter:(counter + 1)] <- c(fit_temp$coefficients[2,2], 
                                              fit_precip$coefficients[2,2])
  env_df$sd_year[counter:(counter + 1)] <- c(sd(temp_resid), sd(precip_resid))
  env_df$sd_season[counter:(counter + 1)] <- c(mean(te$season_temp), 
                                               mean(te$season_precip))
  env_df$kurt[counter:(counter + 1)] <- c(moments::kurtosis(temp_resid), 
                                          moments::kurtosis(precip_resid))
  env_df$skew[counter:(counter + 1)] <- c(moments::skewness(temp_resid), 
                                          moments::skewness(precip_resid))
  env_df$sp_color_year[counter:(counter + 1)] <- c(temp_spec_fit[2] * -1, 
                                                   precip_spec_fit[2] * -1)
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
