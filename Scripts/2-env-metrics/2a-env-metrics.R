# TITLE:            Calculate environmental variability metrics 
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Environmental data csv from script 1b
# DATA OUTPUT:      CSV with summary environmental metrics across all years
# DATE:             May 2023 
# OVERVIEW:         Summarize environmental variability within and between years for temperature and precipitation data. 


# load environmental variables ------------------------------------------------

source("./Scripts/0-config.R")

################
# args:
# - in file
# - out file
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


# read in data -------------------------------------------------

env_out <- read.csv(args[1])


# Env metrics -------------------------------------------------------------

#metrics:
#-mean
#-slope over time (and se)
#-variance (or sd) of residuals (i.e., variability of detrended time series)

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
                     cv_year = NA,
                     sd_season = NA, #mean yearly sd across months
                     cv_season = NA)

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
  
  #fill df
  env_df$lon[counter:(counter + 1)] <- rep(te$lon[1], 2)
  env_df$lat[counter:(counter + 1)] <- rep(te$lat[1], 2)
  env_df$mean[counter:(counter + 1)] <- c(mean(te$mean_temp), mean(te$mean_precip))
  env_df$slope[counter:(counter + 1)] <- c(fit_temp$coefficients[2,1], 
                                           fit_precip$coefficients[2,1])
  env_df$se_slope[counter:(counter + 1)] <- c(fit_temp$coefficients[2,2], 
                                              fit_precip$coefficients[2,2])
  env_df$sd_year[counter:(counter + 1)] <- c(sd(temp_resid), sd(precip_resid))
  env_df$cv_year[counter:(counter + 1)] <- c(sd(temp_resid) / mean(te$mean_temp), 
                                             sd(precip_resid) / mean(te$mean_precip))
  env_df$sd_season[counter:(counter + 1)] <- c(mean(te$season_temp), 
                                               mean(te$season_precip))
  env_df$cv_season[counter:(counter + 1)] <- c(mean(te$season_temp) / mean(te$mean_temp), 
                                               mean(te$season_precip) / mean(te$mean_precip))
  
  counter <- counter + 2
}
print('Run time')
print(proc.time() - time)


# write to csv ------------------------------------------------------------

write.csv(env_df, args[2], row.names = FALSE)
