# TITLE:            Create raster stack  
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Peter Williams, Jeff Dozer, Adriana Uscanga, Lala Kounta, Kelly Kapsar, Phoebe Zarnetske, Pat Bills
# DATA INPUT:       Bird breeding range rasters (6a), env data (2b), model output (5)
# DATA OUTPUT:      Raster stack with summarized gl, dT, and dP values
# DATE:             October 2023 
# OVERVIEW:         Creates raster stack with mean, median, and sd of generation length, delta T, and delta P across all species.

# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)


# Specify top-level directory -------------------------------------------------------

dir <- '/mnt/research/ibeem/variability/'
# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
gl_run_date <- '2023-10-17'


# read in env data --------------------------------------------------------

#read in data to make sure grids are the same
env.dat <- read.csv(paste0(dir, 'data/L2/climate/era5/Env-main.csv')) %>%
  #only 'valid' cells (those over land and > -60S lat)
  dplyr::filter(valid == TRUE) %>%
  dplyr::select(-valid)

#rasterize env data for extraction (to ensure that ranges that don't intersect cell centers are captured)
env.dat.rast <- dplyr::select(env.dat, lon, lat,
                              grep('temp_mean', colnames(env.dat), value = TRUE)) %>%
  terra::rast(crs = "epsg:4326")


# read in bird data -------------------------------------------------------

#from results
bird_df <- readRDS(paste0(dir, 'Results/bird-gl-phylo-vint-', gl_run_date, 
                          '/bird-gl-phylo-vint-data-', gl_run_date, '.rds'))$pro_data
ids <- bird_df$ID

#copy rasters to local
# scp -r 'ccy@rsync.hpcc.msu.edu:/mnt/home/ccy/ibeem/variability/data/L2/range-raster' /Users/caseyyoungflesh/Desktop/
# #filter rasters by valid species (no seabirds/migrants)
# for (i in 1:length(ids))
# {
#   #i <- 1
#   print(paste0(i, ' of ', length(ids)))
#   system(paste0('cp ~/Desktop/range-raster/', ids[i], '-breeding-GenLength.tif ~/Desktop/range-raster-filtered/'))
# }


# get file names ---------------------------------------------------------

#list files of interest
lf <- list.files(paste0(dir, 'data/L2/range-raster/'), full.names = TRUE)

#species ids of rasters
id_ras <- sapply(str_split(sapply(str_split(lf, '/'), tail, 1), 
                           '-'), head, 1)

#get relevant file names
#grep match out of memory...
mrg <- dplyr::left_join(data.frame(id = ids),
                        data.frame(idx = 1:length(id_ras), id = as.numeric(id_ras)), 
                        by = 'id') %>%
  dplyr::arrange(id)

#index for relevant files
lf2 <- lf[mrg$idx]

#files separately
lf_gl <- grep('GenLength', lf2, value = TRUE)
lf_dT <- grep('delta_T', lf2, value = TRUE)
lf_dP <- grep('delta_P', lf2, value = TRUE)


# read in and process -----------------------------------------------------

#stack
gl_stack <- terra::rast(lf_gl)
dT_stack <- terra::rast(lf_dT)
dP_stack <- terra::rast(lf_dP)

# #apply land mask
# gl_stack2 <- terra::mask(tt, env.dat.rast)
# dh_stack2 <- terra::mask(dh_stack, env.dat.rast)

#calculate median - close to mean of logged values for gen length VVV
#https://www.wikiwand.com/en/Log-normal_distribution
med_gl <- terra::app(gl_stack, fun = function(x) median(x, na.rm = TRUE))
med_dT <- terra::app(dT_stack, fun = function(x) median(x, na.rm = TRUE))
med_dP <- terra::app(dP_stack, fun = function(x) median(x, na.rm = TRUE))

names(med_gl) <- 'median_gl'
names(med_dT) <- 'median_dT'
names(med_dP) <- 'median_dP'

#calculate sd
sd_gl <- terra::app(gl_stack, fun = function(x) sd(x, na.rm = TRUE))
sd_dT <- terra::app(dT_stack, fun = function(x) sd(x, na.rm = TRUE))
sd_dP <- terra::app(dP_stack, fun = function(x) sd(x, na.rm = TRUE))

names(sd_gl) <- 'sd_gl'
names(sd_dT) <- 'sd_dT'
names(sd_dP) <- 'sd_dP'

#calculate mean
mn_gl <- terra::app(gl_stack, fun = function(x) mean(x, na.rm = TRUE))
mn_dT <- terra::app(dT_stack, fun = function(x) mean(x, na.rm = TRUE))
mn_dP <- terra::app(dP_stack, fun = function(x) mean(x, na.rm = TRUE))

names(mn_gl) <- 'mean_gl'
names(mn_dT) <- 'mean_dT'
names(mn_dP) <- 'mean_dP'

#calculate number of species in each grid cell
n_sp <- terra::app(gl_stack, fun = function(x) sum(!is.na(x)))
names(n_sp) <- 'n_sp'

#plot
# plot(med_gl)
# plot(n_sp)

#combine into one raster
mrg_ras <- c(med_gl, sd_gl, med_dT, sd_dT, med_dP, sd_dP, mn_gl, mn_dT, mn_dP, n_sp)

#mask out non-land
mrg_ras2 <- terra::mask(mrg_ras, env.dat.rast, 
                  inverse = FALSE)


# terra::writeRaster(med_gl,
#                    filename = paste0(dir, 'data/L3/med-test.tif'),
#                    overwrite = TRUE)
# tt <- terra::rast('~/Downloads/med-test.tif')
# plot(tt)


# save out tifs -----------------------------------------------------------

#species in Sahara (checked with QGIS):
# 1107 - Scissor-tailed kite
# 2018 - Nile Valley sunbird
# 9268 - brown-necked raven
# 10487 - common ostrich -> locations in Sahara where this is the only resident species
terra::writeRaster(mrg_ras2,
                   filename = paste0(dir, 'data/L3/raster-gl-dT-dP-nsp.tif'),
                   overwrite = TRUE)

