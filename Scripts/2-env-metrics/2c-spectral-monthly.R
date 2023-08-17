################
# detrend and calc spectral exponent monthly time series

# need to calculate one chunk at a time (loop through chunked csv files)
# input/output are dir not files
################

# get args fed to script --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- rep(NA, 2)
#args[1] = '/mnt/research/ibeem/variability/data/L1/climate/era5/'
#args[2] = '/mnt/research/ibeem/variability/data/L2/climate/era5/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(data.table)
library(moments)
library(lomb)


# list files -------------------------------------------------

lf <- list.files(args[1], full.names = TRUE)
env_files <- grep('ERA5-monthly-time-series-chunk', lf, value = TRUE)


# loop through files ------------------------------------------------------

counter2 <- 1
for (j in 1:length(env_files))
{
  # j <- 1
  tfile <- read.csv(env_files[j])
  
  #env metrics
  #-spectral exponent ('color' of noise - 0 <= B <= 0.5 = white, 0.5 < B <= 1.5 = red, 1.5 < brown <= 2)
  
  #create unique time point
  ut <- unique(tfile[,c('year', 'month')]) %>%
    dplyr::arrange(year, month)
  ut_df <- data.frame(num = 1:NROW(ut), 
                      ut)
  
  tfile2 <- dplyr::group_by(tfile, lat, lon) %>%
    dplyr::mutate(cell_id = cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(ut_df, by = c('year', 'month')) %>%
    dplyr::arrange(cell_id, year, month) %>%
    data.table::as.data.table()
  
  head(dplyr::arrange(tfile, lat, lon), n = 100)
  
  #unique cells
  uci <- unique(tfile2$cell_id)
  
  #data.frame to fill
  t_df <- data.frame(cell_id = rep(uci, each = 2),
                       var = rep(c('temp', 'precip'), length(uci)),
                       lon = NA,
                       lat = NA,
                       spectral_alpha = NA, #intercept of log10(power) ~ log10(freq)
                       spectral_beta = NA) #spectral exponent
  
  #loop through each cell to calc metrics
  counter <- 1
  for (i in 1:length(uci))
  {
    #i <- 1
    print(paste0('Processing chunk ', j, ' of 10; cell ', i, ' of ', length(uci)))
    
    te <- tfile[cell_id == uci[i],]
    
    #linear model fit (yearly env ~ num)
    fit_temp <- summary(lm(mean_temp ~ num, data = te))
    fit_precip <- summary(lm(mean_precip ~ num, data = te))
    
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
    if (sum(precip_resid) > 0)
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
    
    #fill df
    env_df$lon[counter:(counter + 1)] <- rep(te$lon[1], 2)
    env_df$lat[counter:(counter + 1)] <- rep(te$lat[1], 2)
    env_df$spectral_alpha[counter:(counter + 1)] <- c(temp_spec_fit[1], 
                                                      precip_spec_fit[1])
    env_df$spectral_beta[counter:(counter + 1)] <- c(temp_spec_fit[2] * -1, 
                                                     precip_spec_fit[2] * -1)
    
    counter <- counter + 2
  }
  if (j == 1)
  {
    main_out <- data.frame(cell_id = rep(NA, NROW(env_df) * 10)
                           var = NA,
                           lon = NA,
                           lat = NA,
                           spectral_alpha = NA,
                           spectral_beta = NA)
  }
  
  main_out[counter2:(counter2 + NROW(env_df) -1),] <- env_df
  
  #advance counter
  counter2 <- counter2 + NROW(env_df)
  
}

print('Run time')
print(proc.time() - time)

# write to csv ------------------------------------------------------------

write.csv(main_out, paste0(args[2], 'Env-spectral-exp-monthly.csv'), row.names = FALSE)

