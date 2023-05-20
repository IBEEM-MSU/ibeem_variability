################
# Convert netcdf ERA data into df, averaged of specified months
#
# args:
# - in dir
# - out dir
# - months to avg over (separated by comma, no space)
################


# get args fed to script --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args[1] = in directory (that contains both temp and precip data) - '/mnt/research/ibeem/L0/climate/era5/'
#args[2] = out directory - where to write file - '/mnt/research/ibeem/L1/climate/era5/'
#args[3] = months to average over - 1,2,3,4,5,6,7,8,9,10,11,12

#convert to numeric
MONTHS <- as.numeric(strsplit(args[3], ',')[[1]])

# load packages -----------------------------------------------------------

library(tidyverse)
library(ncdf4)


# function ----------------------------------------------------------------

#function to process temp and precip data
#returns one data.frame for all years and specified cells
#with metric averaged over specified months 

# ARGS:
# startvallat = starting value for lattitude cell (grid is 721 cells tall)
# startvallon = starting value for longitude cell (grid is 1440 cells wide)
# lenlat = number of cells tall for each chunk 
# lenlon = number of cells wide for each chunk
# months = vector of numbered months to average over (e.g., 1 is Jan) - must be consecutive months
proc_fun <- function(startvallat = 500, 
                     startvallon = 500, 
                     lenlon = 10, 
                     lenlat = 10,
                     months = 1)
{
  #list full file paths - corresponding to var of interest
  files_temp <- grep('moda', list.files(paste0(args[1], '/T2m/'), full.names = TRUE), value = TRUE)
  files_precip <- grep('moda', list.files(paste0(args[1], '/PREC/'), full.names = TRUE), value = TRUE)

  #loop through files
  counter <- 1 #where to slot in data into df
  for (i in 1:length(files_temp)) #assuming num of temp files are the same as num precip
  {
    #i <- 1
    print(paste0('Processing year: ', i, ' of ', length(files_temp)))
    
    #open netcdf file  
    tt_temp <- ncdf4::nc_open(files_temp[i])
    tt_precip <- ncdf4::nc_open(files_precip[i])
    
    # Extract data
    # Read lat, lon, and time for each observation
    lon_360 <- ncdf4::ncvar_get(tt_temp, "longitude",
                                start = startvallon, 
                                count = lenlon)
    lat <- ncdf4::ncvar_get(tt_temp, "latitude", verbose = F,
                            start = startvallat, 
                            count = lenlat)
    time <- ncdf4::ncvar_get(tt_temp, "time")[months]
    
    #convert lon to -180 to 180
    lon <- lon_360
    lon[which(lon_360 > 180)] <- (360 - lon_360[which(lon_360 > 180)]) * -1
    
    #extract temp/precip data
    #Replace missing vals in array with NA
    #temp
    var_array_temp <- ncdf4::ncvar_get(tt_temp, "VAR_2T", 
                                       start = c(startvallon, startvallat, min(months)), 
                                       count = c(lenlon, lenlat, length(months))) # 3dim array
    fillvalue_temp <- ncdf4::ncatt_get(tt_temp, "VAR_2T","_FillValue")
    var_array_temp[var_array_temp == fillvalue_temp$value] <- NA
    
    #convert from K to C
    var_array_temp2 <- var_array_temp - 273.15
    
    #precip
    var_array_precip <- ncdf4::ncvar_get(tt_precip, "TCRW", 
                                  start = c(startvallon, startvallat, min(months)),
                                  count = c(lenlon, lenlat, length(months))) # 3dim array
    fillvalue_precip <- ncdf4::ncatt_get(tt_precip, "TCRW","_FillValue")
    var_array_precip[var_array_precip == fillvalue_precip$value] <- NA
    var_array_precip2 <- var_array_precip
    
    # Close netcdf file
    ncdf4::nc_close(tt_temp)
    ncdf4::nc_close(tt_precip)
    
    #covert dates to ymd
    ymd_dates <- lubridate::ymd("1900-01-01") + lubridate::hours(time)
    
    ######
    # #Check - plot rast
    # library(terra)
    # temp_slice <- var_array_temp2[,,1]
    # plot(terra::rast(t(temp_slice)))
    ######
    
    #collapse months into one measure if there is more than one month
    #collapse date into single measure
    if (length(months) > 1)
    {
      var_array_temp3 <- apply(var_array_temp2, c(1, 2), mean)
      var_array_precip3 <- apply(var_array_precip2, c(1, 2), mean)
    } else {
      var_array_temp3 <- var_array_temp2
      var_array_precip3 <- var_array_precip2
    }
    
    # reshape temp array into vector
    var_array_temp3_long <- round(as.vector(var_array_temp3), 2)
    var_array_precip3_long <- round(as.vector(var_array_precip3), 2)
    
    # create a dataframe
    llt <- expand.grid(lon, lat, lubridate::year(ymd_dates)[1])
    tmp_df <- data.frame(llt, var_array_temp3_long, var_array_precip3_long)
    
    colnames(tmp_df) <- c('lon', 'lat', 'year', 'temp', 'precip')
    
    ######
    # #Check - turn df into rast
    # library(terra)
    # pp <- dplyr::filter(tmp_df) %>%
    #   dplyr::select(-year, -precip)
    # plot(terra::rast(pp))
    ######
    
    #create empty df if first file
    if (i == 1)
    {
      out <- data.frame(lon = rep(NA, NROW(tmp_df) * length(files_temp)),
                        lat = NA,
                        year = NA,
                        temp = NA,
                        precip = NA)
    }

    #fill df
    out[counter:(counter + NROW(tmp_df) - 1),] <- tmp_df
    #advance counter
    counter <- counter + NROW(tmp_df)
    
    # #clean up memory
    # rm(llt)
    # rm(var_array)
    # rm(var_array2)
    # rm(var_array2_long)
    # rm(tmp_df)
    # rm(tmp_df2)
    # gc()
  }
  
  return(out)
}


# run function ------------------------------------------------------------

#writes out data.frame with lat, lon, year, month, temp and precip
env_out <- proc_fun(startvallat = 1, 
                     startvallon = 1,
                     lenlon = 1440, #entire world
                     lenlat = 721, #entire world
                     months = MONTHS)


# check -------------------------------------------------------------------

# library(terra)
# pp <- dplyr::filter(env_out, year == 1979) %>%
#   dplyr::select(-year, -precip)
# plot(terra::rast(pp))


# write to csv ------------------------------------------------------------

#ERA5-MONTHS.csv
write.csv(env_out, paste0(args[2], 'ERA5-', 
                           paste0(MONTHS, collapse = '_'), '.csv'), 
          row.names = FALSE)