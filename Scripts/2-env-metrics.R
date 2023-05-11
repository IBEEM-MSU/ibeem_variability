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

#list files for temp
fn_temp <- list.files(temp_data_dir, full.names = TRUE)

#likely need to query specified lat/lon sections of world to fit in memory and speed up processing
#for now, read in first 4 files (entire world for 4 years) and bind together to test env. metric calculation
df_temp <- lapply(fn_temp[1:4], read_csv) %>% 
  dplyr::bind_rows()


# process data ------------------------------------------------------------

#add cell id to df
df_temp2 <- dplyr::group_by(df_temp, lat, lon) %>%
  dplyr::mutate(cell_id = cur_group_id()) %>%
  dplyr::ungroup()

#Summer average (JJA) temp for each year for each cell
#takes several minutes on laptop
mn_JJA_temp <- dplyr::filter(df_temp2, month %in% c(6, 7, 8)) %>%
  #calc mean temp for each cell/year
  dplyr::group_by(cell_id, year) %>%
  dplyr::summarize(mn_temp = mean(temp)) %>%
  dplyr::ungroup() %>%
  #join with cell_id info
  dplyr::left_join(unique(df_temp2[,c('lon', 'lat', 'cell_id')]))

#North America only
mn_JJA_temp_NA <- dplyr::filter(mn_JJA_temp, 
                                lon > -170, lon < -50,
                                lat > 15, lat < 75)

#clean up to reduce memory pressure
rm(df_temp)
rm(df_temp2)
rm(mn_JJA_temp)
gc()

#save out intermediate object to avoid having to reprocess
saveRDS(mn_JJA_temp_NA, paste0(sample_output_dir, 'mn_JJA_temp_NA.rds'))
#mn_JJA_temp_NA <- readRDS(paste0(sample_output_dir, 'mn_JJA_temp_NA.rds'))


# plot --------------------------------------------------------------------

# ggplot(mn_JJA_temp_NA, aes(year, mn_temp, col = factor(cell_id))) + 
#   geom_line(size = 0.25, alpha = 0.3) + 
#   xlab('year') +
#   ylab('temp (C)') +
#   theme_bw() +
#   theme(legend.position = 'none')


# calc metrics ------------------------------------------------------------

#metrics:
#-slope over time
#-variance (or sd) of residuals (i.e., variability of detrended time series)
#-kurtosis (i.e., heavy-tailed ness - or model time series using t-distribution and estimate degrees of freedom parameter)
#-temporal autocorrelation (of residuals - see Leung et al. 2020 Eco Letters for use of this metric for 'predictability')

#unique cells
uci <- unique(mn_JJA_temp_NA$cell_id)

#data.frame to fill
out <- data.frame(cell_id = rep(NA, length(uci)),
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
  
  tt <- dplyr::filter(mn_JJA_temp_NA, cell_id == uci[i])
  
  #linear model fit
  fit <- summary(lm(mn_temp ~ year, data = tt))
  
  #residuals from model
  t_resid <- residuals(fit)
  
  #temporal autocorrelation
  ar_fit <- acf(t_resid, lag.max = 5, plot = FALSE)
  
  #fill df
  out$cell_id[i] <- uci[i]
  out$lon[i] <- tt$lon[1]
  out$lat[i] <- tt$lat[2]
  out$slope[i] <- fit$coefficients[2,1]
  out$se_slope[i] <- fit$coefficients[2,2]
  out$sd_resid[i] <- sd(t_resid)
  out$kurt[i] <- moments::kurtosis(t_resid)
  out$skew[i] <- moments::skewness(t_resid)
  out$rho_l1[i] <- ar_fit$acf[2]
  out$rho_l2[i] <- ar_fit$acf[3]
  out$rho_l3[i] <- ar_fit$acf[4]
  out$rho_l4[i] <- ar_fit$acf[5]
  out$rho_l5[i] <- ar_fit$acf[6]
}

#save out intermediate object to avoid having to reprocess
saveRDS(out, paste0(sample_output_dir, 'env_var_out.rds'))
#out <- readRDS(paste0(sample_output_dir, 'env_var_out.rds'))


# explore -----------------------------------------------------------------

#convert to raster to visualize
dplyr::select(out, lon, lat, slope) %>%
  terra::rast() %>%
  plot(main = 'slope')

dplyr::select(out, lon, lat, sd_resid) %>%
  terra::rast() %>%
  plot(main = 'sd resid')

dplyr::mutate(out, sl_sd = slope / sd_resid) %>%
dplyr::select(lon, lat, sl_sd) %>%
  terra::rast() %>%
  plot(main = 'slope/sd')
