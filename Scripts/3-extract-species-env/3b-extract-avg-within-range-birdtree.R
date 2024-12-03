# TITLE:            Extract environmental data for bird ranges
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Range shapefiles from from script 1d, Env variability from 2b
# DATA OUTPUT:      CSV files containing environmental variables across bird ranges for subsets of species 
# DATE:             September 2023 
# OVERVIEW:         Extract environmental data for bird ranges for subsets of species 


rm(list = ls())


# load environment variables ------------------------------------------------

source('./Scripts/0-config.R')


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)


# Get the current file to process -----------------------------------------

# file.name <- commandArgs(trailingOnly = TRUE)
# Testing
file.name <- 'BTIDsPiece-1.rda'
# if(length(file.name) == 0) base::stop('Need to give the file name to process')


# Read in data ------------------------------------------------------------

# Loads in current set of ids (a vector called ids)
load(paste0(dir, 'data/L1/range/bird-breeding/', file.name))

# Climate data (takes a couple of minutes to load)
# only valid cells (land and N of 60S lat)
env.dat <- read.csv(paste0(dir, 'data/L2/climate/era5/Env-main.csv')) %>%
  #only 'valid' cells (those over land and > -60S lat)
  dplyr::filter(valid == TRUE) %>%
  dplyr::select(-valid)

#rasterize env data for extraction (so ensure that ranges that don't intersect cell centers are captured)
env.dat.rast <- dplyr::select(env.dat, lon, lat,
                               grep('temp', colnames(env.dat), value = TRUE),
                               grep('precip', colnames(env.dat), value = TRUE)) %>%
  terra::rast(crs = "epsg:4326")

#just mean temp, precip for each cell
tpd.dat.rast <- env.dat.rast[[c('temp_mean', 'precip_mean')]]


# summarize ---------------------------------------------------------------

env.dat.rast[[c('temp_sd_year', 'temp_sd_season',
             'precip_cv_year', 'precip_cv_season')]] %>%
  terra::global(median, na.rm = TRUE)


# Loop through all the species in the current file ------------------------

# Make data frame to hold everything
cn <- c('ID', names(env.dat.rast), 
        paste0(names(env.dat.rast), '_sd_space'),
        'temp_rng_space', 'precip_rng_space',
        'cen_lon', 'cen_lat', 'range_size_km2')
env.out <- data.frame(matrix(NA, nrow = length(ids), ncol = length(cn)))
colnames(env.out) <- cn                      

counter <- 1
for (i in 1:length(ids))
{
  #i <- 1
  print(paste0("Currently on species ", i, " out of ", length(ids)))
  curr.sp <- ids[i]
  curr.range <- sf::st_read(paste0(dir, 'data/L1/range/bird-breeding/birdtree-', 
                                   curr.sp, '-breeding.shp'),
                            quiet = TRUE)
  
  #if more than one polygon (resident + breeding ranges), merge them
  if (NROW(curr.range) > 1)
  {
    #to address 'duplciate vertex' errors:
    # https://github.com/r-spatial/sf/issues/1762
    #to address lines across globe for species that cross date line
    # https://gis.stackexchange.com/questions/462335/st-union-in-r-creates-artifacts-for-global-data
    tcrs <- sf::st_crs(curr.range) # save crs for later
    sf::st_crs(curr.range) <- NA  # set crs to missing
    curr.range2 <- try(sf::st_union(curr.range))
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
  
  # exclude Inf for those cells that divide by 0 for CV
  # keep only finite values within 3 MAD of median if range covers more than 5 cells
  filt_fun <- function(x)
  {
    x2 <- x[is.finite(x)]
    if (length(x2) > 5)
    {
      x3 <- x2[(x2 < (median(x2) + 3 * mad(x2))) & (x2 > (median(x2) - 3 * mad(x2)))]
    } else {
      x3 <- x2
    }
    return(x3)
  }
  
  sd_filt_fun <- function(x)
  {
    x2 <- filt_fun(x)
    if (length(x2) > 1)
    {
      x3 <- sd(x2)
    } else {
      x3 <- 0
    }
    return(x3)
  }
  
  # Extract mean env value across range from raster
  env.mn.vals <- terra::extract(env.dat.rast,
                 terra::vect(curr.range2),
                 touches = TRUE,
                 fun = function(x) mean(filt_fun(x)))

  #sd of env values across space
  env.sd.vals <- terra::extract(env.dat.rast,
                            terra::vect(curr.range2),
                            touches = TRUE,
                            fun = function(x) sd_filt_fun(x))
  
  #range of mean values across space
  rng.vals <- terra::extract(tpd.dat.rast,
                             terra::vect(curr.range2),
                             touches = TRUE,
                             fun = function(x) diff(range(filt_fun(x))))
  
  # reproject to laea (equal area)
  curr.range.tr <- sf::st_transform(curr.range2, crs = "+proj=laea")
  # get centroid - ignore warning
  cen_ll <- sf::st_centroid(curr.range.tr) %>%
    sf::st_transform(4326) %>%
    sf::st_coordinates() %>%
    as.numeric()
  # get range size - km^2
  rsize_km2 <- as.numeric(round(sf::st_area(curr.range.tr) / 1000^2, 0))
  
  #fill df - current id, env vals, centroid, range size
  env.out[counter,] <- c(curr.sp, env.mn.vals[-1], 
                         env.sd.vals[-1], rng.vals[-1],
                         cen_ll, rsize_km2)
  
  #advance counter
  counter <- counter + 1
  
  #clean up memory
  rm(curr.range, curr.range2, curr.range.tr, env.mn.vals, env.sd.vals, rng.vals)
  gc()
}

#write out file
save(env.out, file = paste0(dir, 'data/L2/range-env-pieces-birdtree/new-summarized-data-piece-BT-', 
				stringr::str_extract(file.name, '\\d+')))
