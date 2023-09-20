# 3b: extract env metrics across species ranges

rm(list = ls())

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
# file.name <- 'BTIDsPiece-14.rda'
if(length(file.name) == 0) base::stop('Need to give the file name to process')


# Read in data ------------------------------------------------------------

# Loads in current set of ids (a vector called ids)
load(paste0(dir, 'data/L1/range/bird-breeding/', file.name))

# Climate data (takes a couple of minutes to load)
# only valid cells (land and N of -60S lat)
env.dat <- read.csv(paste0(dir, 'data/L2/climate/era5/Env-main.csv')) %>%
  #only 'valid' cells (those over land and > -60S lat)
  dplyr::filter(valid == TRUE) %>%
  dplyr::select(-valid)

#rasterize env data for extraction (so ensure that ranges that don't intersect cell centers are captured)
env.dat.rast <- dplyr::select(env.dat, lon, lat,
                               grep('temp', colnames(env.dat), value = TRUE),
                               grep('precip', colnames(env.dat), value = TRUE),
                              grep('dhi', colnames(env.dat), value = TRUE)) %>%
  terra::rast(crs = "epsg:4326")

#just mean temp precip, and cum DHI for each cell
tpd.dat.rast <- env.dat.rast[[c('temp_mean', 'precip_mean', 'dhi_cum_mean')]]


# Loop through all the species in the current file ------------------------

# Make data frame to hold everything
cn <- c('ID', names(env.dat.rast), 
        'temp_rng_space', 'precip_rng_space', 'dhi_cum_rng_space',
        'temp_sd_space', 'precip_sd_space', 'dhi_cum_sd_space',
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
  
  # Extract median env value across range from raster
  env.vals <- terra::extract(env.dat.rast,
                 terra::vect(curr.range2),
                 touches = TRUE,
                 fun = function(x) median(x, na.rm = TRUE))
  
  #range of mean values across space
  rng.vals <- terra::extract(tpd.dat.rast,
                             terra::vect(curr.range2),
                             touches = TRUE,
                             fun = function(x) diff(range(x, na.rm = TRUE)))
  
  #sd of mean values across space
  sd.vals <- terra::extract(tpd.dat.rast,
                             terra::vect(curr.range2),
                             touches = TRUE,
                             fun = function(x) sd(x, na.rm = TRUE))
  
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
  env.out[counter,] <- c(curr.sp, env.vals[-1], 
                         rng.vals[-1], sd.vals[-1],
                         cen_ll, rsize_km2)
  
  #advance counter
  counter <- counter + 1
  
  #clean up memory
  rm(curr.range, curr.range2, curr.range.tr, env.vals, rng.vals, sd.vals)
  gc()
}

#write out file
save(env.out, file = paste0(dir, 'data/L2/range-env-pieces-birdtree/summarized-data-piece-BT-', 
				stringr::str_extract(file.name, '\\d+')))
