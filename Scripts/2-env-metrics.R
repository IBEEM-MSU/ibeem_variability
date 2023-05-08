################
# Calc env variability metrics
#
# NOTES:
# - determine how best to chunk up world for metric calculation
# - write results back to database
# - create tif files of all env metrics
# - can we fit regresion models with t-distributed errors (estimating and not specificy degrees of freedom parameter) in R at a resonable speed? (seems yes -> https://stats.stackexchange.com/questions/117980/regression-with-t-distributed-errors-and-massrlm)
################


# specify dir -------------------------------------------------------------

#path for data on CY machine - remember trailing slash
temp_data_dir <- '~/Downloads/T2m/'
precip_data_dir <- '~/Downloads/PREC/'


# load packages -----------------------------------------------------------

library(tidyverse)


# read in data -------------------------------------------------

#list files for temp
fn <- list.files(temp_data_dir, full.names = TRUE)

#likely need to query specified lat/lon sections of world to fit in memory and speed up processing
#for now, read in first 5 files (entire world for 5 years) and bind together to test env. metric calculation
df <- lapply(fn[1:5],read_csv) %>% 
  dplyr::bind_rows()


# process data ------------------------------------------------------------

#Summer average (JJA) temp for each year for each cell
mn_JJA_temp <- dplyr::filter(df, month %in% c(6, 7, 8)) %>%
  #add cell id
  dplyr::group_by(lat, lon) %>%
  dplyr::mutate(cell_id = cur_group_id()) %>%
  dplyr::ungroup() %>%
  #calc mean temp for each cell/year
  dplyr::group_by(cell_id, year) %>%
  dplyr::summarize(mn_temp = mean(temp))


# calc metrics ------------------------------------------------------------

#metrics:
#-slope over time
#-variance (or sd) of residuals (i.e., variability of detrended time series)
#-kurtosis (heavy-tailed ness - order model time series using t-distribution)
#-temporal autocorrelation (of residuals) - AR1 model (see Leung et al. 2020 Eco Letters for use of this metric for 'predictability)

#unique cells
uci <- unique(mn_JJA_temp$cell_id)

#data.frame to fill
out <- data.frame(cell_id = rep(NA, length(uci)).
                  lon = NA,
                  lat = NA,
                  slope = NA, #slope of linear model
                  se_slope = NA, #standard error of slope
                  sd_resid = NA, #sd of residuals for linear model
                  kurt = NA, #kurtosis
                  skew = NA, #skewness
                  rho = NA) #temporal autocorrelation

#loop through each cell to calc metrics
for (i in 1:length(uci))
{
  #i <- 1
  tt <- dplyr::filter(mn_JJA_temp, cell_id == uci[i])
  
  #linear model fit
  fit <- summary(lm(temp ~ year, data = tt))
  #residuals from model
  tresid <- residuals(fit)
  
  #AR1 model
  ar_fit <- ar(tresid)
  
  #fill df
  out$cell_id[i] <- uci[i]
  out$lon[i] <- tt$lon[1]
  out$lat[i] <- tt$lon[2]
  out$slope[i] <- fit$coefficients[1,2]
  out$se_slope[i] <- fit$coefficients[2,2]
  out$sd_resid[i] <- sd(t_resid)
  out$kurt[i] <- moments::kurtosis(t_resid)
  out$skew[i] <- moments::skewness(t_resid)
  out$rho[i] <- XXXX
}


# explore -----------------------------------------------------------------

#convert to raster to visualize
out

# #Check - turn df into rast
# library(terra)
# pp <- dplyr::filter(tmp_df2, month == 2) %>%
#   dplyr::select(-year, -month)
# plot(terra::rast(pp))
