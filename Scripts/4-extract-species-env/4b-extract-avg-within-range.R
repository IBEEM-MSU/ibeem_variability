# 4b-example-species-extract.R: a script to start thinking about how we'll want
#                              to calculate metrics within each species range
#                              and output all the associated information for 
#                              running models of the form trait ~ climate
rm(list = ls())
library(tidyverse)
library(sf)
library(terra)

# Specify directory -------------------------------------------------------
BL.dir <- '/mnt/research/ibeem/variability/data/L1/range/bird-breeding/'
env.dir <- '/mnt/research/ibeem/variability/data/L2/climate/era5/'
out.dir <- '/mnt/research/ibeem/variability/data/L2/range-env-pieces/'

# Get the current file to process -----------------------------------------
file.name <- commandArgs(trailingOnly = TRUE)
# Testing
# file.name <- 'BLIDsPiece-14.rda'
if(length(file.name) == 0) base::stop('Need to give the file name to process')


# Read in data ------------------------------------------------------------
# Loads in current set of ids (a vector called ids)
load(paste0(BL.dir, file.name))
# Climate data (takes a couple of minutes to load
# NOTE: can change this to the GAM stuff if that's what we want.

#only valid cells (land and N of -60S lat)
env.dat <- read.csv(paste0(env.dir, 'Env-main.csv')) %>%
  #only 'valid' cells (those over land and > -60S lat)
  dplyr::filter(valid == TRUE) %>%
  dplyr::select(-valid)

#rasterize env data for extraction (so ensure that ranges that don't intersect cell centers are captured)
env.dat.rast <- dplyr::select(env.dat, lon, lat,
                               grep('temp', colnames(env.dat), value = TRUE),
                               grep('precip', colnames(env.dat), value = TRUE),
                              grep('mys', colnames(env.dat), value = TRUE)) %>%
  terra::rast(crs = "epsg:4326")


# Loop through all the species in the current file ------------------------

# Make data frame to hold everything
cn <- c('ID', names(env.dat.rast), 'cen_lon', 'cen_lat', 'range_size_km2')
env.out <- data.frame(matrix(NA, nrow = length(ids), ncol = length(cn)))
colnames(env.out) <- cn                      

counter <- 1
for (i in 1:length(ids)) {
  #i <- 1
  print(paste0("Currently on species ", i, " out of ", length(ids)))
  curr.sp <- ids[i]  
  curr.range <- sf::st_read(paste0(BL.dir, curr.sp, '-breeding.shp'))
  
  # Extract mean env value across range from raster
  env.vals <- terra::extract(env.dat.rast,
                 terra::vect(curr.range),
                 touches = TRUE,
                 fun = function(x) mean(x, na.rm = TRUE))

  # reproject to laea (equal area)
  curr.range.tr <- sf::st_transform(curr.range, crs = "+proj=laea")
  # get centroid - ignore warning
  cen_ll <- sf::st_centroid(curr.range.tr) %>%
    sf::st_transform(4326) %>%
    sf::st_coordinates() %>%
    as.numeric()
  # get range size - km^2
  rsize_km2 <- as.numeric(round(sf::st_area(curr.range.tr) / 1000^2, 0))
  
  
  # #TO DO: RASTERIZE SHP FILE
  # #create empty rast
  # er <- temp.dat.rast[[1]]
  # er[] <- NA
  # #rasterize range with same crs as temp raster
  # curr.range.rast <- terra::rasterize(vect(curr.range), tt)
  #
  # REPLACE RASTER VALUES (1'S) WITH GENERATION LENGTH (WILL NEED TO READ IN GEN LENGTH DATA FROM 4C...)
  # save out as single-species tifs
  
  
  #fill df - current id, env vals, centroid, range size
  env.out[counter,] <- c(curr.sp, env.vals[-1], cen_ll, rsize_km2)
  
  #advance counter
  counter <- counter + 1
  
  #clean up memory
  rm(curr.range, curr.range.tr, env.vals)
  gc()
}

#write out file
save(env.out, file = paste0(out.dir, 'summarized-data-piece-', 
				stringr::str_extract(file.name, '\\d+')))
