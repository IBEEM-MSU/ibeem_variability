################
# Convert netcdf ERA data into df
# calculate env variabilty metrics
#
# NOTES:
# -could mask water (has values)
################


# specify dir -------------------------------------------------------------

#path for data on CY machine - remember trailing slash
temp_data_dir <- '~/Downloads/T2m/'
precip_data_dir <- '~/Downloads/PREC/'

#where to write files to
out_dir_temp <- '~/Desktop/temp/'
out_dir_precip <- '~/Desktop/precip/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(ncdf4)
library(data.table)


# function ----------------------------------------------------------------

#function to process temp OR precip data
#returns one data.frame for all years and specified cells
#with metric averaged over specified months 

# ARGS:
# startvallat = starting value for lattitude cell (grid is 721 cells tall)
# startvallon = starting value for longitude cell (grid is 1440 cells wide)
# lenlat = number of cells tall for each chunk 
# lenlon = number of cells wide for each chunk
# months = vector of numbered months to average over (e.g., 1 is Jan)
proc_fun <- function(type = 'TEMP', 
                     startvallat = 500, 
                     startvallon = 500, 
                     lenlon = 10, 
                     lenlat = 10,
                     months = 1)
{
  #list full file paths - corresponding to var of interest
  if (type == 'TEMP')
  {
    files <- grep('moda', list.files(temp_data_dir, full.names = TRUE), value = TRUE)
  }
  if (type == 'PREC')
  {
    files <- grep('moda', list.files(precip_data_dir, full.names = TRUE), value = TRUE)
  }
  if (type != 'TEMP' & type != 'PREC')
  {
    stop('Arg "type" must be either "TEMP" or "PREC"')
  }
  
  
  #loop through files
  counter <- 1 #where to slot in data into df
  for (i in 1:length(files))
  {
    #i <- 1
    print(paste0('Processing ' , type, ' file: ', i, ' of ', length(files)))
    
    #open netcdf file  
    tt <- ncdf4::nc_open(files[i])
    
    # Extract data
    # Read lat, lon, and time for each observation
    lon_360 <- ncdf4::ncvar_get(tt, "longitude",
                                start = startvallon, 
                                count = lenlon)
    lat <- ncdf4::ncvar_get(tt, "latitude", verbose = F,
                            start = startvallat, 
                            count = lenlat)
    time <- ncdf4::ncvar_get(tt, "time")[months]
    
    #convert lon to -180 to 180
    lon <- lon_360
    lon[which(lon_360 > 180)] <- (360 - lon_360[which(lon_360 > 180)]) * -1
    
    #extract temp/precip data
    #Replace missing vals in array with NA
    if (type == 'TEMP')
    {
      var_array <- ncdf4::ncvar_get(tt, "VAR_2T", 
                                    start = c(startvallon, startvallat, min(months)), 
                                    count = c(lenlon, lenlat, max(months))) # 3dim array
      # var_array <- ncdf4::ncvar_get(tt, "VAR_2T")
      fillvalue <- ncdf4::ncatt_get(tt, "VAR_2T","_FillValue")
      var_array[var_array == fillvalue$value] <- NA
      
      #convert from K to C
      var_array2 <- var_array - 273.15
    }
    if (type == 'PREC')
    {
      var_array <- ncdf4::ncvar_get(tt, "TCRW", 
                                    start = c(startvallon, startvallat, min(months)),
                                    count = c(lenlon, lenlat, max(months))) # 3dim array
      #var_array <- ncdf4::ncvar_get(tt, "TCRW")
      fillvalue <- ncdf4::ncatt_get(tt, "TCRW","_FillValue")
      var_array[var_array == fillvalue$value] <- NA
      var_array2 <- var_array
    }
    
    # Close netcdf file
    ncdf4::nc_close(tt)
    
    #covert dates to ymd
    ymd_dates <- lubridate::ymd("1900-01-01") + lubridate::hours(time)
    
    ######
    # #Check - plot rast
    # library(terra)
    # temp_slice <- var_array2[,,1]
    # plot(terra::rast(t(temp_slice)))
    ######
    
    #collapse months into one measure if there is more than one month
    #collapse date into single measure
    if (length(dim(var_array2)) > 2)
    {
      var_array3 <- apply(var_array2, c(1, 2), mean)
    } else {
      var_array3 <- var_array2
    }
    
    # reshape temp array into vector
    var_array3_long <- round(as.vector(var_array3), 2)
    
    # create a dataframe
    llt <- expand.grid(lon, lat, lubridate::year(ymd_dates))
    tmp_df <- data.frame(llt, var_array3_long)
    
    if (type == 'TEMP')
    {
      colnames(tmp_df) <- c('lon', 'lat', 'year', 'temp')
    }
    if (type == 'PREC')
    {
      colnames(tmp_df) <- c('lon', 'lat', 'year', 'precip')
    }
    
    ######
    # #Check - turn df into rast
    # library(terra)
    # pp <- dplyr::filter(tmp_df) %>%
    #   dplyr::select(-year)
    # plot(terra::rast(pp))
    ######
    
    #create empty df if first file
    if (i == 1)
    {
      if (type == 'TEMP')
      {
        out <- data.frame(lon = rep(NA, NROW(tmp_df) * length(files)),
                          lat = NA,
                          year = NA,
                          temp = NA)
      }
      if (type == 'PREC')
      {
        out <- data.frame(lon = rep(NA, NROW(tmp_df) * length(files)),
                          lat = NA,
                          year = NA,
                          precip = NA)
      }
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

#writes out data.frame with lat, lon, year, month, and env variable (temp or precip)
temp_out <- proc_fun(type = 'TEMP', 
                     startvallat = 1, 
                     startvallon = 1,
                     lenlon = 1440, #entire world
                     lenlat = 721, #entire world
                     months = 1) # JAN

#precip_out <- proc_fun(type = 'PREC')


# check -------------------------------------------------------------------

# library(terra)
# pp <- dplyr::filter(temp_out, year == 1979) %>%
#   dplyr::select(-year)
# plot(terra::rast(pp))


# write to csv ------------------------------------------------------------

#VAR-MONTHS.csv
write.csv(temp_out, paste0(out_dir_temp, 'Temp-', 
                           paste0(as.character(months), collapse = '_'), '.csv'), 
          row.names = FALSE)
write.csv(precip_out, paste0(out_dir_precip, 'Precip-', 
                             paste0(as.character(months), collapse = '_'), '.csv'), 
          row.names = FALSE)
