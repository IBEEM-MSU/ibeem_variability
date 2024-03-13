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

#landmask
lm <- sf::st_read(paste0(dir, 'data/L0/landmask/ne_10m_land.shp')) %>%
  dplyr::filter(featurecla == 'Land')


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
  dplyr::mutate(valid = TRUE)


# merge env data ----------------------------------------------------------

#merge masked env var with spectral exp
env_mrg <- dplyr::left_join(env_var, se, 
                          by = c('var', 'lon', 'lat')) %>%
  #land mask
  dplyr::left_join(t_cval, by = c('cell_id', 'lon', 'lat')) %>%
  #DHI
  dplyr::left_join(dhi_df, by = c('lon', 'lat')) %>%
na_idx <- which(is.na(env_mrg$valid))
#set all ocean and land < 60 S lat to FALSE for field valid
env_mrg$valid[na_idx] <- FALSE

# saveRDS(env_mrg, '~/tt.rds')


# merge and write out -----------------------------------------------------

#temp, precip, and DHI together
tt_temp <- dplyr::filter(env_mrg, var == 'temp')
names(tt_temp) <- c('cell_id', 'var', 'lon', 'lat', 
                    paste0('temp_', names(tt_temp)[-c(1:4)]))
tt_precip <- dplyr::filter(env_mrg, var == 'precip') %>%
  dplyr::select(-c(var, valid, cell_id, lon, lat))
names(tt_precip) <- paste0('precip_', names(tt_precip))

#merge
tt_mrg <- cbind(tt_temp, tt_precip) %>%
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
                valid)

#write to csv
write.csv(tt_mrg, paste0(dir, 'data/L2/climate/era5/Env-main.csv'), 
          row.names = FALSE)
