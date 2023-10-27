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

# #read in life history and other traits
# #use to filter species of interest (only resident, terrestrial birds)
'%ni%' <- Negate('%in%')

or_excl <- c('Sphenisciformes', #penguins 
             'Procellariiformes', #tubenoses
             'Pelecaniformes', #pelicans
             'Suliformes', #gannets/boobies
             'Phaethontiformes', #tropicbirds
             'Charadriiformes')#, #skuas, gulls, terns, skimmers, auks
#'Anseriformes', #waterfowl
#'Ciconiiformes', #storks
#'Gaviiformes', #aquatic birds (loons and divers)
#'Gruiformes', #cranes, rails - Family Psophiidae are not waterbirds, but there are few members (~6 species)
#'Phoenicopteriformes', #flamingos and relatives
#'Podicipediformes') #grebes

'%ni%' <- Negate('%in%')
bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data-birdtreeX.csv')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::filter(Order %ni% or_excl,
                Migration == 1) %>%
  #filter out precip outlier
  dplyr::filter(precip_cv_season < 2.5) %>%
  dplyr::mutate(lMass = log(Mass),
                lGL = log(GenLength),
                lAb = log(Modeled_age_first_breeding),
                # lAb = log(Measured_age_first_breeding),
                lMl = log(Modeled_max_longevity),
                species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  #drop duplicated species (for now)
  dplyr::group_by(Birdtree_name) %>%
  dplyr::slice_head() %>%
  dplyr::ungroup() %>%
  #ex = number of years before temp exceeds 2 sd
  #AND
  #delta_t = how much temp will change in 1 generation (in sds)
  #AND
  #delta_haldane = how much change (in sd) per generation
  #AND
  #n_gen = how many gens before temp will exceed 2 sd
  dplyr::mutate(ex = 2 * temp_sd_year / temp_slope,
                ex_season = 2 * temp_sd_season / temp_slope,
                delta_t = temp_slope / temp_sd_year * GenLength,
                delta_t_season = temp_slope / temp_sd_year * GenLength,
                #(degrees / year) * (sd / degrees) * (year / gen) = sd / gen
                delta_haldane = (temp_slope / temp_sd_year) * GenLength,
                delta_haldane_season = (temp_slope / temp_sd_season) * GenLength,
                n_gen = ex / GenLength)

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
lf_dh <- grep('delta_haldane', lf2, value = TRUE)
lf_ab <- grep('-Ab', lf2, value = TRUE)
lf_ml <- grep('-Ml', lf2, value = TRUE)
lf_cs <- grep('-Cs', lf2, value = TRUE)
lf_s <- grep('-S', lf2, value = TRUE)
lf_pc1 <- grep('-LH_PC1', lf2, value = TRUE)
lf_pc2 <- grep('-LH_PC2', lf2, value = TRUE)
lf_pc3 <- grep('-LH_PC3', lf2, value = TRUE)


# read in and process -----------------------------------------------------

#stack
gl_stack <- terra::rast(lf_gl)
dh_stack <- terra::rast(lf_dh)
ab_stack <- terra::rast(lf_ab)
ml_stack <- terra::rast(lf_ml)
cs_stack <- terra::rast(lf_cs)
s_stack <- terra::rast(lf_s)
pc1_stack <- terra::rast(lf_pc1)
pc2_stack <- terra::rast(lf_pc2)
pc3_stack <- terra::rast(lf_pc3)

#currently number of names is several species too long (result of duplicate species after merging -- just a couple)
# names(gl_stack) <- unique(mrg$id)
# names(dh_stack) <- unique(mrg$id)
# names(ab_stack) <- unique(mrg$id)
# names(ml_stack) <- unique(mrg$id)
# names(cs_stack) <- unique(mrg$id)
# names(s_stack) <- unique(mrg$id)

# #apply land mask
# gl_stack2 <- terra::mask(tt, env.dat.rast)
# dh_stack2 <- terra::mask(dh_stack, env.dat.rast)

#calculate median - close to mean of logged values for gen length VVV
#https://www.wikiwand.com/en/Log-normal_distribution
med_gl <- terra::app(gl_stack, fun = function(x) median(x, na.rm = TRUE))
med_dh <- terra::app(dh_stack, fun = function(x) median(x, na.rm = TRUE))
med_ab <- terra::app(ab_stack, fun = function(x) median(x, na.rm = TRUE))
med_ml <- terra::app(ml_stack, fun = function(x) median(x, na.rm = TRUE))
med_cs <- terra::app(cs_stack, fun = function(x) median(x, na.rm = TRUE))
med_s <- terra::app(s_stack, fun = function(x) median(x, na.rm = TRUE))
med_pc1 <- terra::app(pc1_stack, fun = function(x) median(x, na.rm = TRUE))
med_pc2 <- terra::app(pc2_stack, fun = function(x) median(x, na.rm = TRUE))
med_pc3 <- terra::app(pc3_stack, fun = function(x) median(x, na.rm = TRUE))
names(med_gl) <- 'median_gl'
names(med_dh) <- 'median_dh'
names(med_ab) <- 'median_ab'
names(med_ml) <- 'median_ml'
names(med_cs) <- 'median_cs'
names(med_s) <- 'median_s'
names(med_pc1) <- 'median_pc1'
names(med_pc2) <- 'median_pc2'
names(med_pc3) <- 'median_pc3'

# #calculate mean
# mn_gl <- terra::app(gl_stack, fun = function(x) mean(x, na.rm = TRUE))
# mn_dh <- terra::app(dh_stack, fun = function(x) mean(x, na.rm = TRUE))
# names(med_gl) <- 'mean_gl'
# names(med_dh) <- 'mean_dh'

#calculate sd
sd_gl <- terra::app(gl_stack, fun = function(x) sd(x, na.rm = TRUE))
sd_dh <- terra::app(dh_stack, fun = function(x) sd(x, na.rm = TRUE))
sd_ab <- terra::app(ab_stack, fun = function(x) sd(x, na.rm = TRUE))
sd_ml <- terra::app(ml_stack, fun = function(x) sd(x, na.rm = TRUE))
sd_cs <- terra::app(cs_stack, fun = function(x) sd(x, na.rm = TRUE))
sd_s <- terra::app(s_stack, fun = function(x) sd(x, na.rm = TRUE))
sd_pc1 <- terra::app(pc1_stack, fun = function(x) sd(x, na.rm = TRUE))
sd_pc2 <- terra::app(pc2_stack, fun = function(x) sd(x, na.rm = TRUE))
sd_pc3 <- terra::app(pc3_stack, fun = function(x) sd(x, na.rm = TRUE))
names(sd_gl) <- 'sd_gl'
names(sd_dh) <- 'sd_dh'
names(sd_ab) <- 'sd_ab'
names(sd_ml) <- 'sd_ml'
names(sd_cs) <- 'sd_cs'
names(sd_s) <- 'sd_s'
names(sd_pc1) <- 'sd_pc1'
names(sd_pc2) <- 'sd_pc2'
names(sd_pc3) <- 'sd_pc3'

#calculate number of species in each grid cell
n_sp <- terra::app(gl_stack, fun = function(x) sum(!is.na(x)))
names(n_sp) <- 'n_sp'

#plot
# plot(med_gl)
# plot(n_sp)

#combine into one raster
mrg_ras <- c(med_gl, sd_gl, med_dh, sd_dh, med_ab, sd_ab, 
             med_ml, sd_ml, med_cs, sd_cs, med_s, sd_s, 
             med_pc1, sd_pc1, med_pc2, sd_pc2, med_pc3, sd_pc4, 
             n_sp)

#mask out non-land
mrg_ras2 <- terra::mask(mrg_ras, env.dat.rast, 
                  inverse = FALSE)


# save out tifs -----------------------------------------------------------

#species in Sahara (checked with QGIS):
# 1107 - Scissor-tailed kite
# 2018 - Nile Valley sunbird
# 9268 - brown-necked raven
# 10487 - common ostrich -> locations in Sahara where this is the only resident species
terra::writeRaster(mrg_ras2,
                   filename = paste0(dir, 'data/L3/raster-LH-nsp.tif'),
                   overwrite = TRUE)

