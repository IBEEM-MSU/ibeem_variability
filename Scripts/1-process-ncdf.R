################
# Convert netcdf ERA data into df
# outputs csv for each year, then combines csvs
#
# NOTES:
# -could mask water (has values)
# -could change long values (0 to 360 now)
# -add processing for Precip
################


# specify dir -------------------------------------------------------------

#path for data on CY machine - remember trailing slash
data_dir <- '~/Downloads/T2m/'
out_dir <- '~/Desktop/temp/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(ncdf4)


# notes -------------------------------------------------------------------

#why are some files named e5 and some e5p?

#list only files of interest
files <- grep('moda', list.files(data_dir, full.names = TRUE), value = TRUE)


# process ------------------------------------------------------------

#loop through files
# counter <- 1 #for use in combining all years into single df
for (i in 1:length(files))
{
  #i <- 1
  print(paste0('Processing file: ', i, ' of ', length(files)))
  
  #open netcdf file  
  tt <- ncdf4::nc_open(files[i])
  
  # Extract data
  # Read lat, lon, and time for each observation
  lon <- ncdf4::ncvar_get(tt, "longitude")
  lat <- ncdf4::ncvar_get(tt, "latitude", verbose = F)
  time <- ncdf4::ncvar_get(tt, "time")
  
  #temp data
  temp_array <- ncdf4::ncvar_get(tt, "VAR_2T") # 3dim array
  #dims: lon, lat, time
  #dim(temp_array)
  
  #Replace missing vals in array with NA
  fillvalue <- ncdf4::ncatt_get(tt, "VAR_2T","_FillValue")
  temp_array[temp_array == fillvalue$value] <- NA
  
  # Close netcdf file
  ncdf4::nc_close(tt)
  
  # process data
  #convert from K to C
  temp_array_C <- temp_array - 273.15
  
  #covert dates to ymd
  ymd_dates <- lubridate::ymd("1900-01-01") + lubridate::hours(time)
  
  ######
  # #Check
  # library(terra)
  # temp_slice <- temp_array_C[,,1]
  # plot(terra::rast(t(temp_slice)))
  ######
  
  #adapted from: https://pjbartlein.github.io/REarthSysSci/netCDF.html#introduction
  # reshape temp array into vector
  temp_array_C_long <- round(as.vector(temp_array_C), 2)
  
  # create a dataframe
  llt <- as.matrix(expand.grid(lon, lat, ymd_dates))
  tmp_df <- data.frame(cbind(llt, temp_array_C_long))
  colnames(tmp_df) <- c('lon', 'lat', 'date', 'temp')
  
  #rearrange rows, add year and month
  tmp_df2 <- dplyr::arrange(tmp_df, lon, lat, date) %>%
    dplyr::mutate(year = lubridate::year(date),
                  month = lubridate::month(date)) %>%
    dplyr::select(-date)

  #clean up memory
  rm(llt)
  rm(temp_array)
  rm(temp_array_C)
  rm(temp_array_C_long)
  rm(tmp_df)
  gc()
  
  ######
  # #Check - turn df into rast
  # #library(terra)
  # pp <- dplyr::filter(tmp_df2, year == 1979, month == 1) %>%
  #   dplyr::select(-year, -month)
  # plot(terra::rast(pp))
  ######
  
  #for one df all files (hits mem limit on laptop)
  # #create empty df if first file
  # if (i == 1)
  # {
  #   out <- data.frame(lon = rep(NA, NROW(tmp_df2) * length(files)),
  #                     lat = NA,
  #                     year = NA,
  #                     month = NA,
  #                     temp = NA)
  # }
  # 
  # #fill df
  # out[counter:(counter + NROW(tmp_df2)) - 1] <- tmp_df2
  # #advance counter
  # counter <- counter + NROW(tmp_df2)
  
  #write out on csv per year
  write.csv(tmp_df2, paste0(out_dir, 'T2m-', tmp_df2$year[1], '.csv'))
  
  rm(tmp_df2)
  gc()
}


# combine csvs ------------------------------------------------------------

fn <- list.files(out_dir, full.names = TRUE)
tt <- read.csv(fn[1])
head(tt)
