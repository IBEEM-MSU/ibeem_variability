################
# Mask out land and everything below 60 S lat for env var data and calc PCA
#
################


# Specify dir --------------------------------------------------

#path CY machine
# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
dir <- '/mnt/research/ibeem/variability/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)


# read in data -------------------------------------------------

#environmental variability - add relative slope
env_var <- read.csv(paste0(dir, 'data/L2/climate/era5/Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv')) %>%
  dplyr::mutate(rel_slope = slope / sd_year)

#monthly spectral exponent
lf <- list.files(paste0(dir, 'data/L1/climate/era5/'), full.names = TRUE)
ef <- grep('spectral', lf, value = TRUE)
#read all files in and bind - change name to match env_var
se <- do.call(rbind, lapply(ef, read.csv)) %>%
  #NOTE this cell_id is not the same as env_var cell_id (so remove it)
  dplyr::select(-cell_id) %>%
  reshape2::melt(id.vars = c('lon', 'lat'),
                 value.name = 'spectral_exp', variable.name = 'var')
se$var <- as.character(se$var)
se$var[which(se$var == 'spectral_beta_temp')] <- 'temp'
se$var[which(se$var == 'spectral_beta_precip')] <- 'precip'

#landmask
lm <- sf::st_read(paste0(dir, 'data/L0/landmask/ne_10m_land.shp')) %>%
  dplyr::filter(featurecla == 'Land')

#read in DHI rast
dhi_rast <- terra::rast(paste0(dir, 'data/L1/DHI/DHI_rast.tif'))


# apply land mask -------------------------------------------------------------

#just need one mask for both vars
frast <- dplyr::filter(env_var, var == 'temp') %>%
  dplyr::select(lon, lat, cell_id) %>%
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
  #dplyr::filter(lat > -66.5, lat < 66.5) %>%
  dplyr::mutate(valid = TRUE)


# process DHI -------------------------------------------------------------

#resample DHI raster to env.dat.rast
dhi_rast_rs <- terra::resample(dhi_rast, frast, method = 'average')

#df format
dhi_df <- data.frame(terra::crds(dhi_rast_rs), 
                     terra::values(dhi_rast_rs, na.rm = TRUE)) %>%
  dplyr::rename(lat = y, lon = x)

#plot
# dplyr::select(dhi_df, lon, lat, dhi_cum_mean) %>% 
#   terra::rast() %>% 
#   plot()

#where cumulative NDVI is very low, insert NA into CV year
dhi_df$dhi_cv_year[which(dhi_df$dhi_cum_mean < 0.01)] <- NA

# dplyr::select(dhi_df, lon, lat, dhi_cv_year) %>%
#   terra::rast() %>%
#   plot()


# merge env data ----------------------------------------------------------

#merge masked env var with spectral exp
env_mrg <- dplyr::left_join(env_var, se, 
                          by = c('var', 'lon', 'lat')) %>%
  #land mask
  dplyr::left_join(t_cval, by = c('cell_id', 'lon', 'lat')) %>%
  #DHI
  dplyr::left_join(dhi_df, by = c('lon', 'lat')) %>%
  dplyr::rename(sp_color_month = spectral_exp)
na_idx <- which(is.na(env_mrg$valid))
#set all ocean and land < 60 S lat to FALSE for field valid
env_mrg$valid[na_idx] <- FALSE

# saveRDS(env_mrg, '~/tt.rds')

# plt <- dplyr::filter(env_mrg, var == 'temp', valid == TRUE) %>%
#   dplyr::select(lon, lat, sp_color_year) %>%
#   terra::rast(crs = "epsg:4326")
# plot(plt)


# merge and write out -----------------------------------------------------

#temp, precip, and DHI together
tt_temp <- dplyr::filter(env_mrg, var == 'temp') %>%
  dplyr::select(-dhi_cum_mean, -dhi_cv_year, -dhi_cv_season)
names(tt_temp) <- c('cell_id', 'var', 'lon', 'lat', 
                    paste0('temp_', names(tt_temp)[-c(1:4)]))
tt_precip <- dplyr::filter(env_mrg, var == 'precip') %>%
  dplyr::select(-c(var, valid, cell_id, lon, lat,
                   -dhi_cum_mean, -dhi_cv_year, -dhi_cv_season))
names(tt_precip) <- paste0('precip_', names(tt_precip))
tt_dhi <- dplyr::filter(env_mrg, var == 'temp') %>%
  dplyr::select(dhi_cum_mean, dhi_cv_year, dhi_cv_season)

#merge
tt_mrg <- cbind(tt_temp, tt_precip, tt_dhi) %>%
  dplyr::rename(valid = temp_valid) %>%
  dplyr::select(-var) %>%
  #reorder
  dplyr::select(cell_id, lon, lat, 
                temp_mean, precip_mean, 
                temp_sd_year, precip_sd_year,
                temp_cv_year, precip_cv_year,
                temp_sd_season, precip_sd_season,
                temp_cv_season, precip_cv_season,
                temp_slope, precip_slope,
                temp_se_slope, precip_se_slope,
                temp_rel_slope, precip_rel_slope,
                temp_rng_season, precip_rng_season,
                temp_kurt, precip_kurt,
                temp_skew, precip_skew,
                temp_sp_color_year, precip_sp_color_year,
                temp_sp_color_month, precip_sp_color_month,
                temp_rho_l1, precip_rho_l1,
                temp_rho_l2, precip_rho_l2,
                temp_rho_l3, precip_rho_l3,
                temp_rho_l4, precip_rho_l4,
                temp_rho_l5, precip_rho_l5,
                dhi_cum_mean, dhi_cv_year, dhi_cv_season,
                valid)

#write to csv
write.csv(tt_mrg, paste0(dir, 'data/L2/climate/era5/Env-main.csv'), 
          row.names = FALSE)
