# TITLE:            Create rasterized range maps for individual species 
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Bird breeding ranges (1d), env data (2b), model output (5)
# DATA OUTPUT:      Rasterized maps of variables across bird ranges
# DATE:             October 2023 
# OVERVIEW:         Creates rasters for generation length, delta T, and delta P for each bird species



rm(list = ls())

# load environmental variables ------------------------------------------------

source("./Scripts/0-config.R")

# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)

# Get the current file to process -----------------------------------------

file.name <- commandArgs(trailingOnly = TRUE)
# Testing
# file.name <- 'BTIDsPiece-9.rda'
if(length(file.name) == 0) base::stop('Need to give the file name to process')


# Read in data ------------------------------------------------------------

# Loads in current set of ids (a vector called ids)
load(paste0(dir, 'data/L1/range/bird-breeding/', file.name))

#read in data to make sure grids are the same
env.dat <- read.csv(paste0(dir, 'data/L2/climate/era5/Env-main.csv')) %>%
  #only 'valid' cells (those over land and > -60S lat)
  dplyr::filter(valid == TRUE) %>%
  dplyr::select(-valid)

#rasterize env data for extraction (to ensure that ranges that don't intersect cell centers are captured)
env.dat.rast <- dplyr::select(env.dat, lon, lat,
                              grep('temp_mean', colnames(env.dat), value = TRUE)) %>%
  terra::rast(crs = "epsg:4326")

#read in life history and other traits
#from results
bird_df <- readRDS(paste0(dir, 'Results/bird-gl-phylo-vint-', gl_run_date, 
               '/bird-gl-phylo-vint-data-', gl_run_date, '.rds'))$pro_data


# loop through ranges -----------------------------------------------------

#from 3b...R
counter <- 1
for (i in 1:length(ids))
{
  #i <- 1
  print(paste0("Currently on species ", i, " out of ", length(ids)))
  curr.sp <- ids[i]
  
  # if(file.exists(paste0(dir, 'data/L1/range/bird-breeding/', 
  #                       curr.sp, '-breeding.shp'))){
  curr.range <- sf::st_read(paste0(dir, 'data/L1/range/bird-breeding/birdtree-', 
                                   curr.sp, '-breeding.shp'),
                              quiet = TRUE)
  # }else(print(paste0(dir, 'data/L1/range/bird-breeding/', 
  #                    curr.sp, '-breeding.shp does not exist.')))
  
  #filter trait for species of interest (use ID)
  tdf <- dplyr::filter(bird_df, ID == curr.sp)
  
  #only proceed if species was part of final df
  if (NROW(tdf) > 0)
  {
    #if more than one polygon (resident + breeding ranges), merge them
    if (NROW(curr.range) > 1)
    {
      #to address 'duplciate vertex' errors:
      # https://github.com/r-spatial/sf/issues/1762
      #to address lines across globe for species that cross date line
     # https://gis.stackexchange.com/questions/462335/st-union-in-r-creates-artifacts-for-global-data
      tcrs <- sf::st_crs(curr.range) # save crs for later
      sf::st_crs(curr.range) <- NA  # set crs to missing
      curr.range2 <- try(sf::st_union(curr.range)) #join
      sf::st_crs(curr.range2) <- tcrs #reassign crs
      
      if (inherits(curr.range2, "try-error"))
      {
        tcrs <- sf::st_crs(curr.range)
        sf::st_crs(curr.range)
        curr.range2 <- try(sf::st_union(sf::st_make_valid(curr.range)))
        sf::st_crs(curr.range2) <- tcrs
        
        if (inherits(curr.range2, "try-error"))
        {
          # to avoid some 'duplciate vertex' errors:
          # https://github.com/r-spatial/sf/issues/1762
          sf::sf_use_s2(FALSE)
          curr.range2 <- sf::st_union(curr.range)
          sf::sf_use_s2(TRUE)
        }
      }
    } else {
      curr.range2 <- curr.range
    }
  
    #rasterize range with same crs as temp raster
    curr.range.rast <- terra::rasterize(terra::vect(curr.range2), 
                                        env.dat.rast, 
                                        touches = TRUE)
    #two layers
    curr.range.rast2 <- c(curr.range.rast, curr.range.rast,
                          curr.range.rast)
    
    #change layer names
    names(curr.range.rast2) <- c('GenLength', 'delta_T', 'delta_P')
    
    curr.range.rast2[['GenLength']][curr.range.rast2['GenLength'] == 1] <- exp(tdf$lGL)
    curr.range.rast2[['delta_T']][curr.range.rast2['delta_T'] == 1] <- tdf$temp_delta
    curr.range.rast2[['delta_P']][curr.range.rast2['delta_P'] == 1] <- tdf$precip_delta
    
    # save out as single-species tifs (one per metric)
    terra::writeRaster(curr.range.rast2[['GenLength']], 
                       filename = paste0(dir, 'data/L2/range-raster/', 
                                         gsub(' ', '_', curr.sp),
                                         '-breeding-GenLength.tif'),
                       overwrite = TRUE)
    
    terra::writeRaster(curr.range.rast2[['delta_T']],
                       filename = paste0(dir, 'data/L2/range-raster/', 
                                         gsub(' ', '_', curr.sp),
                                       '-breeding-delta_T.tif'),
                       overwrite = TRUE)
    
    terra::writeRaster(curr.range.rast2[['delta_P']],
                       filename = paste0(dir, 'data/L2/range-raster/', 
                                         gsub(' ', '_', curr.sp),
                                         '-breeding-delta_P.tif'),
                       overwrite = TRUE)
  }
}

