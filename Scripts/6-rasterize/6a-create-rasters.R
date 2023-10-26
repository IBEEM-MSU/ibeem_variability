#############################
# Create rasterized range maps for birds
#
# Rasters for: 
# - Gl (gen length)
# - delta haldane
# - Ab (age first breeding)
# - Cs (clutch size)
# - Ml (maximum longevity)
# - S (survival)
#############################


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)


# Specify top-level directory -------------------------------------------------------

dir <- '/mnt/research/ibeem/variability/'
# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# Get the current file to process -----------------------------------------

file.name <- commandArgs(trailingOnly = TRUE)
# Testing
# file.name <- 'BLIDsPiece-9.rda'
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


# loop through ranges -----------------------------------------------------

#from 3b...R
counter <- 1
for (i in 1:length(ids))
{
  #i <- 42
  print(paste0("Currently on species ", i, " out of ", length(ids)))
  curr.sp <- ids[i]  
  curr.range <- sf::st_read(paste0(dir, 'data/L1/range/bird-breeding/', 
                                   curr.sp, '-breeding.shp'),
                            quiet = TRUE)
  
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
                          curr.range.rast, curr.range.rast,
                          curr.range.rast, curr.range.rast)
    
    #change layer names
    names(curr.range.rast2) <- c('GenLength', 'delta_haldane',
                                 'Ab', 'Ml', 'Cs', 'S')
    
    curr.range.rast2[['GenLength']][curr.range.rast2['GenLength'] == 1] <-  tdf$GenLength
    curr.range.rast2[['delta_haldane']][curr.range.rast2['delta_haldane'] == 1] <- tdf$delta_haldane
    curr.range.rast2[['Ab']][curr.range.rast2['Ab'] == 1] <- tdf$Modeled_age_first_breeding
    curr.range.rast2[['Ml']][curr.range.rast2['Ml'] == 1] <- tdf$Modeled_max_longevity
    curr.range.rast2[['Cs']][curr.range.rast2['Cs'] == 1] <- tdf$Measured_clutch_size
    curr.range.rast2[['S']][curr.range.rast2['S'] == 1] <- tdf$Modeled_survival
    
    # save out as single-species tifs (one per metric)
    terra::writeRaster(curr.range.rast2[['GenLength']], 
                       filename = paste0(dir, 'data/L2/range-raster/', 
                                         gsub(' ', '_', curr.sp),
                                         '-breeding-GenLength.tif'),
                       overwrite = TRUE)
    
    terra::writeRaster(curr.range.rast2[['delta_haldane']],
                       filename = paste0(dir, 'data/L2/range-raster/', 
                                         gsub(' ', '_', curr.sp),
                                       '-breeding-delta_haldane.tif'),
                       overwrite = TRUE)
    
    terra::writeRaster(curr.range.rast2[['Ab']],
                       filename = paste0(dir, 'data/L2/range-raster/', 
                                         gsub(' ', '_', curr.sp),
                                         '-breeding-Ab.tif'),
                       overwrite = TRUE)
    
    terra::writeRaster(curr.range.rast2[['Ml']],
                       filename = paste0(dir, 'data/L2/range-raster/', 
                                         gsub(' ', '_', curr.sp),
                                         '-breeding-Ml.tif'),
                       overwrite = TRUE)
    
    terra::writeRaster(curr.range.rast2[['Cs']],
                       filename = paste0(dir, 'data/L2/range-raster/', 
                                         gsub(' ', '_', curr.sp),
                                         '-breeding-Cs.tif'),
                       overwrite = TRUE)
    
    terra::writeRaster(curr.range.rast2[['S']],
                       filename = paste0(dir, 'data/L2/range-raster/', 
                                         gsub(' ', '_', curr.sp),
                                         '-breeding-S.tif'),
                       overwrite = TRUE)
  }
}

