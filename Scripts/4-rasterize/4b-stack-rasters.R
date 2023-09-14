#############################
# Stack rasters
#
#############################


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)


# Specify top-level directory -------------------------------------------------------

dir <- '/mnt/research/ibeem/variability/'
# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# read in bird data -------------------------------------------------------

# #read in life history and other traits
# #use to filter species of interest (only resident, terrestrial birds)
'%ni%' <- Negate('%in%')

or_excl <- c('Sphenisciformes', #penguins 
             'Procellariiformes', #tubenoses
             'Pelecaniformes', #pelicans
             'Suliformes', #gannets/boobies
             'Phaethontiformes', #tropicbirds
             'Charadriiformes', #skuas, gulls, terns, skimmers, auks
             'Anseriformes', #waterfowl
             'Ciconiiformes', #storks
             'Gaviiformes', #aquatic birds (loons and divers)
             'Gruiformes', #cranes, rails - Family Psophiidae are not waterbirds, but there are few members (~6 species)
             'Phoenicopteriformes', #flamingos and relatives
             'Podicipediformes') #grebes

bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data.csv')) %>%
  dplyr::arrange(Accepted_name) %>%
  dplyr::filter(Order %ni% or_excl,
                Migration == 1)

ids <- bird_df$ID


# read in env data for landmask -----------------------------------------------

# # Climate data (takes a couple of minutes to load
# # only valid cells (land and N of -60S lat)
# env.dat <- read.csv(paste0(dir, 'data/L2/climate/era5/Env-main.csv')) %>%
#   #only 'valid' cells (those over land and > -60S lat)
#   dplyr::filter(valid == TRUE) %>%
#   dplyr::select(-valid)
# 
# #rasterize env data for land mask
# env.dat.rast <- dplyr::select(env.dat, lon, lat,
#                               grep('temp_mean', colnames(env.dat), value = TRUE)) %>%
#   terra::rast(crs = "epsg:4326")


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

#two files separately
lf_gl <- grep('GenLength', lf2, value = TRUE)
lf_dh <- grep('delta_haldane', lf2, value = TRUE)


# read in and process -----------------------------------------------------

#stack
gl_stack <- terra::rast(lf_gl)
dh_stack <- terra::rast(lf_dh)
names(gl_stack) <- unique(mrg$id)
names(dh_stack) <- unique(mrg$id)

# #apply land mask
# gl_stack2 <- terra::mask(tt, env.dat.rast)
# dh_stack2 <- terra::mask(dh_stack, env.dat.rast)

gl_stack2 <- gl_stack
dh_stack2 <- dh_stack

#calculate median - close to mean of logged values for gen length VVV
#https://www.wikiwand.com/en/Log-normal_distribution
med_gl <- terra::app(gl_stack2, fun = function(x) median(x, na.rm = TRUE))
med_dh <- terra::app(dh_stack2, fun = function(x) median(x, na.rm = TRUE))
names(med_gl) <- 'median_gl'
names(med_dh) <- 'median_dh'

# #calculate mean
# mn_gl <- terra::app(gl_stack2, fun = function(x) mean(x, na.rm = TRUE))
# mn_dh <- terra::app(dh_stack2, fun = function(x) mean(x, na.rm = TRUE))
# names(med_gl) <- 'mean_gl'
# names(med_dh) <- 'mean_dh'

#calculate sd
sd_gl <- terra::app(gl_stack2, fun = function(x) sd(x, na.rm = TRUE))
sd_dh <- terra::app(dh_stack2, fun = function(x) sd(x, na.rm = TRUE))
names(sd_gl) <- 'sd_gl'
names(sd_dh) <- 'sd_dh'

#calculate number of species in each grid cell
n_sp <- terra::app(gl_stack2, fun = function(x) sum(!is.na(x)))
names(n_sp) <- 'n_sp'

#plot
# plot(med_gl)
# plot(n_sp)

#combine into one raster
mrg_ras <- c(med_gl, sd_gl, med_dh, sd_dh, n_sp)


# save out tifs -----------------------------------------------------------

terra::writeRaster(mrg_ras,
                   filename = paste0(dir, 'data/L3/raster-gl-dh-nsp.tif'),
                   overwrite = TRUE)
