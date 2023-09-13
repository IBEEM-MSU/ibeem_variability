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
# #use to filter species of interest (e.g., no seabirds)
'%ni%' <- Negate('%in%')
bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data.csv')) %>%
  dplyr::arrange(Accepted_name) %>%
  dplyr::filter(!is.na(GenLength),
                Habitat %ni% c('Marine', 'Coastal'),
                Primary.Lifestyle %ni% ('Aquatic'),
                Migration == 1)
  
nms <- unique(gsub(' ', '_', unique(bird_df$Accepted_name)))


# get file names ---------------------------------------------------------

#list files of interest
lf <- list.files(paste0(dir, 'data/L2/range-raster/'), full.names = TRUE)

#species names of rasters
sp_ras <- sapply(str_split(sapply(str_split(lf, '/'), tail, 1), 
                           '-'), head, 1)

#get relevant file names
#grep match out of memory...
mrg <- dplyr::left_join(data.frame(idx = 1:length(sp_ras), name = sp_ras), 
                        data.frame(name = nms), by = 'name')

#index for relevant files
lf2 <- lf[mrg$idx]

#two files separately
lf_gl <- grep('GenLength', lf2, value = TRUE)
lf_dh <- grep('delta_haldane', lf2, value = TRUE)


# read in and process -----------------------------------------------------

gl_stack <- terra::rast(lf_gl)
dh_stack <- terra::rast(lf_dh)

#calculate median
med_gl <- terra::app(gl_stack, fun = function(x) median(x, na.rm = TRUE))
med_dh <- terra::app(dh_stack, fun = function(x) median(x, na.rm = TRUE))
names(med_gl) <- 'median_gl'
names(med_dh) <- 'median_dh'

#calculate number of species in each grid cell
n_sp <- terra::app(gl_stack, fun = function(x) sum(!is.na(x)))
names(n_sp) <- 'n_sp'

#plot
# plot(med_gl)
# plot(n_sp)


# save out tifs -----------------------------------------------------------

terra::writeRaster(med_gl,
                   filename = paste0(dir, 'data/L3/raster-gl.tif'),
                   overwrite = TRUE)

terra::writeRaster(med_dh,
                   filename = paste0(dir, 'data/L3/raster-dh.tif'),
                   overwrite = TRUE)

terra::writeRaster(n_sp,
                   filename = paste0(dir, 'data/L3/raster-nsp.tif'),
                   overwrite = TRUE)
