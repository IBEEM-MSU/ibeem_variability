# TITLE:            Environmental variation over time
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Environmental data csv from script 1b
# DATA OUTPUT:      Trend in env var
# DATE:             July 2024
# OVERVIEW:         Calculate trends in interannual env var using moving window approach


# load environment variables ------------------------------------------------

source("./Scripts/0-config.R")


# load packages -----------------------------------------------------------

library(tidyverse)
library(data.table)
library(zoo)
library(terra)


# read in data -------------------------------------------------

env_out <- read.csv(paste0(dir, 'data/L1/climate/era5/ERA5-1_2_3_4_5_6_7_8_9_10_11_12.csv')) %>%
  dplyr::group_by(lat, lon) %>%
  dplyr::mutate(cell_id = cur_group_id()) %>%
  dplyr::ungroup()

#landmask
lm <- sf::st_read(paste0(dir, 'data/L0/landmask/ne_10m_land.shp')) %>%
  dplyr::filter(featurecla == 'Land')

#apply land mask
#just need one mask for both vars
frast <- dplyr::select(env_out, lon, lat, cell_id) %>%
  terra::rast(crs = "epsg:4326") %>%
  #mask out water
  terra::mask(terra::vect(lm), touches = TRUE)

#coords (lat/lon)
t_crds <- terra::crds(frast)
#values
t_vals <- data.frame(terra::values(frast)) %>%
  dplyr::filter(!is.na(cell_id))

#masked df (valid col)
#combine, filter above -60 lat, add valid col
t_cval <- cbind(t_crds, t_vals) %>%
  dplyr::rename(lat = y, lon = x) %>%
  dplyr::filter(lat > -60) %>%
  dplyr::mutate(valid = TRUE)

#merge env data
# remove masked pixels 
env_mrg <- dplyr::left_join(env_out, t_cval, by = c('cell_id', 'lon', 'lat'))

na_idx <- which(is.na(env_mrg$valid))
#set all ocean and land < 60 S lat to FALSE for field valid
env_mrg$valid[na_idx] <- FALSE

env_mrg2 <- dplyr::filter(env_mrg, valid == TRUE) %>%
  dplyr::arrange(cell_id, year) %>%
  dplyr::select(lat, lon, year, 
                mean_temp, mean_precip,
                season_temp, season_precip) %>% 
  dplyr::group_by(lat, lon) %>%
  dplyr::mutate(cell_id = cur_group_id()) %>%
  dplyr::ungroup() %>%
  data.table::as.data.table()


# calc residuals of linear model over time --------------------------------

#for each cell
t_resid <- env_mrg2[, 
                    residuals(lm(mean_temp ~ year, .SD)), by = cell_id]
p_resid <- env_mrg2[, 
                    residuals(lm(mean_precip ~ year, .SD)), by = cell_id]


# rolling sd/cv for each cell -------------------------------------------------

#for each cell
W <- 10 #window size
t_rsd <- t_resid[, zoo::rollapplyr(V1, W, 
                                   function(x) sd(x, na.rm = TRUE)), 
                 by = cell_id]
p_rsd <- p_resid[, zoo::rollapplyr(V1, W, 
                                   function(x) sd(x, na.rm = TRUE) / 
                                     mean(x, na.rm = TRUE)), 
                 by = cell_id]

#add time
t_rsd2 <- dplyr::group_by(t_rsd, cell_id) %>%
  dplyr::mutate(t = 1:(72 - W + 1)) %>%
  dplyr::ungroup() %>%
  data.table::as.data.table()
p_rsd2 <- dplyr::group_by(p_rsd, cell_id) %>%
  dplyr::mutate(t = 1:(72 - W + 1)) %>%
  dplyr::ungroup() %>%
  data.table::as.data.table()


# fit lm for rolling sd/cv ---------------------------------------------------

t_fit <- t_rsd2[, 
                {
                  fit <- lm(V1 ~ t, .SD)
                  list(int = fit$coefficients[1], sl = fit$coefficients[2])
                },
                  by = cell_id]

#remove outliers (3 MAD) and NaN
p_rsd3 <- dplyr::filter(p_rsd2, is.finite(V1)) %>%
  data.table::as.data.table() %>%
  dplyr::mutate(mad3 = 3 * mad(V1)) %>%
  dplyr::filter(V1 < (median(V1) + mad3) & V1 > (median(V1) - mad3))

p_fit <- p_rsd3[, 
                {
                  fit <- lm(V1 ~ t, .SD)
                  list(int = fit$coefficients[1], sl = fit$coefficients[2])
                },
                by = cell_id]


# plots inter ------------------------------------------------------------

#histograms of slopes
hist(t_fit$sl, breaks = 50,
     main = 'Slope: tau_inter_T',
     xlab = 'Slope')
med_t <- median(t_fit$sl)
abline(v = med_t, col = 'red', lwd = 3)
# legend('topleft',
#        legend = paste0('Median = ', round(med_t, 4)),
#        bty = 'n')

hist(p_fit$sl, breaks = 50,
     main = 'Slope: tau_inter_P',
     xlab = 'Slope')
med_p <- median(p_fit$sl)
abline(v = med_p, col = 'red', lwd = 3)
# legend('topleft',
#        legend = paste0('Median = ', round(med_p, 4)),
#        bty = 'n')


# save out inter -------------------------------------------------------------

#residuals from lm
saveRDS(t_resid, paste0(dir, 'data/L2/climate/era5/Temp-resid.rds'))
saveRDS(p_resid, paste0(dir, 'data/L2/climate/era5/Precip-resid.rds'))

#rolling sd
saveRDS(t_rsd2, paste0(dir, 'data/L2/climate/era5/Temp-roll-sd.rds'))
saveRDS(p_rsd3, paste0(dir, 'data/L2/climate/era5/Precip-roll-sd.rds'))

#model fit
saveRDS(t_fit, paste0(dir, 'data/L2/climate/era5/Temp-sd-time.rds'))
saveRDS(p_fit, paste0(dir, 'data/L2/climate/era5/Precip-sd-time.rds'))


# trend in intra ----------------------------------------------------------

ts_fit <- env_mrg2[,{
  fit <- lm(season_temp ~ year, .SD)
  list(int = fit$coefficients[1], sl = fit$coefficients[2])
},
by = cell_id]

ps_fit <- env_mrg2[,{
  fit <- lm(season_precip ~ year, .SD)
  list(int = fit$coefficients[1], sl = fit$coefficients[2])
},
by = cell_id]


# plots intra ------------------------------------------------------------

#histograms of slopes
hist(ts_fit$sl, breaks = 50,
     main = 'Slope: tau_intra_T',
     xlab = 'Slope')
med_ts <- median(ts_fit$sl)
abline(v = med_ts, col = 'red', lwd = 3)

hist(ps_fit$sl, breaks = 50,
     main = 'Slope: tau_intra_P',
     xlab = 'Slope')
med_ps <- median(ps_fit$sl)
abline(v = med_ps, col = 'red', lwd = 3)


# save out intra -------------------------------------------------------------

#model fit
saveRDS(ts_fit, paste0(dir, 'data/L2/climate/era5/Temp-intra-sd-time.rds'))
saveRDS(ps_fit, paste0(dir, 'data/L2/climate/era5/Precip-intra-sd-time.rds'))
