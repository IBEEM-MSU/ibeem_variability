#stack bird range rasters (with gen time) and make plots

#rasterize range map


# read in data ------------------------------------------------------------

#read in data to make sure grids are the same
env.dat <- read.csv(paste0(dir, 'data/L2/climate/era5/Env-main.csv')) %>%
  #only 'valid' cells (those over land and > -60S lat)
  dplyr::filter(valid == TRUE) %>%
  dplyr::select(-valid)

#rasterize env data for extraction (so ensure that ranges that don't intersect cell centers are captured)
env.dat.rast <- dplyr::select(env.dat, lon, lat,
                              grep('temp', colnames(env.dat), value = TRUE)) %>%
  terra::rast(crs = "epsg:4326")


# loop through ranges -----------------------------------------------------

#from 3b...R
counter <- 1
for (i in 1:length(ids))
{
  #i <- 1
  print(paste0("Currently on species ", i, " out of ", length(ids)))
  curr.sp <- ids[i]  
  curr.range <- sf::st_read(paste0(dir, 'data/L1/range/bird-breeding/', 
                                   curr.sp, '-breeding.shp'),
                            quiet = TRUE)
  
  #if more than one polygon (resident + breeding ranges), merge them
  if (NROW(curr.range) > 1)
  {
    #to address 'duplciate vertex' errors:
    # https://github.com/r-spatial/sf/issues/1762
    curr.range2 <- try(sf::st_union(curr.range))
    
    if (inherits(curr.range2, "try-error"))
    {
      curr.range2 <- try(sf::st_union(sf::st_make_valid(curr.range)))
      
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

  ######################################
  # #TO DO: RASTERIZE SHP FILES
  # #create empty rast
  # er <- env.dat.rast[[1]]
  # er[] <- NA
  # #rasterize range with same crs as temp raster
  # 
  # ADD RASTER VALUES (1'S) WITH GENERATION LENGTH (WILL NEED TO READ IN GEN LENGTH DATA FROM 4C...)
  # tt represents values of interest (e.g., gen length)
  # curr.range.rast <- terra::rasterize(vect(curr.range), tt)
  #
  # save out as single-species tifs
  # terra::writeRaster()
  #
  # will stack individual .tifs in another script to get median gen length, etc. 
  ######################################
}
