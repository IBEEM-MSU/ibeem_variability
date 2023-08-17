################
# Convert netcdf ERA data into df, monthly temp/precip values for each cell
# Use this to calculate spectral exponents in other script
#
# args:
# - in dir
# - out dir
# - chunk number
################


# get args fed to script --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- rep(NA, 3)
#args[1] = '/mnt/research/ibeem/variability/data/L0/climate/era5/'
#args[2] = '/mnt/research/ibeem/variability/data/L1/climate/era5/'
#args[3] = 1


# load packages -----------------------------------------------------------

library(tidyverse)
library(ncdf4)


# function ----------------------------------------------------------------

#function to process temp and precip data
#returns one data.frame for all years and specified cells
#with metric averaged over specified months 

# ARGS:
# startvallat = starting value for latitude cell (grid is 721 cells tall)
# startvallon = starting value for longitude cell (grid is 1440 cells wide)
# lenlat = number of cells tall for each chunk 
# lenlon = number of cells wide for each chunk
# months = vector of numbered months to average over (e.g., 1 is Jan) - must be consecutive months
proc_fun <- function(startvallat = 500, 
                     startvallon = 500, 
                     lenlon = 10, 
                     lenlat = 10,
                     months = 1:12)
{
  #list full file paths - corresponding to var of interest
  files_temp <- grep('moda', list.files(paste0(args[1], '/T2m/'), full.names = TRUE), value = TRUE)
  files_precip <- grep('moda', list.files(paste0(args[1], '/PREC/'), full.names = TRUE), value = TRUE)
  
  #loop through files
  counter <- 1 #where to slot in data into df
  for (i in 1:length(files_temp)) #assuming num of temp files are the same as num precip
  {
    #i <- 1
    
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
    
    #squar root transform precip
    var_array_precip2 <- sqrt(var_array_precip)
    
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
    
    # reshape temp array into vector
    var_array_temp3_long <- round(as.vector(var_array_temp2), 2)
    var_array_precip3_long <- round(as.vector(var_array_precip2), 2)
    
    
    # create a dataframe
    llt <- expand.grid(lon, 
                       lat, 
                       lubridate::year(ymd_dates)[1], 
                       lubridate::month(ymd_dates))
    tmp_df <- data.frame(llt, var_array_temp3_long, var_array_precip3_long)
    
    colnames(tmp_df) <- c('lon', 'lat', 'year', 'month', 'temp', 'precip')

    #create empty df if first file
    if (i == 1)
    {
      out <- data.frame(lon = rep(NA, NROW(tmp_df) * length(files_temp)),
                        lat = NA,
                        year = NA,
                        month = NA,
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

#401 chunks
#21 different jobs
#20 chunks per job, but 1 job with only 1 chunk

#each job runs through a limited set of lats and all lons
lenlat <- 721 %/% 20 #entire world / number chunks
lenlon <- 1440 / 20 #entire world / number chunks
lat_start <- (as.numeric(args[3]) * lenlat) - lenlat + 1
lon_start <- 1

#in last lat chunk remainder of len / # chunks
if (as.numeric(args[3]) == 21)
{
  lenlat <- 1
}


chunk_counter <- 1
counter2 <- 1
for (k in 1:20) #lon
{
  # k <- 1
  #writes out data.frame with lat, lon, year, month, temp and precip
  env_out <- proc_fun(startvallon = lon_start,
                      startvallat = lat_start, 
                      lenlon = lenlon,
                      lenlat = lenlat)
  
  if (k == 1)
  {
    #create unique time point
    ut <- unique(env_out[,c('year', 'month')]) %>%
      dplyr::arrange(year, month)
    ut_df <- data.frame(num = 1:NROW(ut), 
                        ut)
  }
  
  env_out2 <- dplyr::group_by(env_out, lat, lon) %>%
    dplyr::mutate(cell_id = cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(ut_df, by = c('year', 'month')) %>%
    dplyr::arrange(cell_id, year, month) %>%
    data.table::as.data.table()
  
  #unique cells
  uci <- unique(env_out2$cell_id)
  
  #data.frame to fill
  t_df <- data.frame(cell_id = uci,
                     lon = NA,
                     lat = NA,
                     spectral_beta_temp = NA, #spectral exponent
                     spectral_beta_precip = NA) #spectral exponent
  
  #loop through each cell to calc metrics
  counter <- 1
  for (i in 1:length(uci))
  {
    #i <- 1
    print(paste0('Processing chunk ', chunk_counter, ' of 20; cell ', 
                 i, ' of ', length(uci)))
    
    #just one cell
    te <- env_out2[cell_id == uci[i],]
    
    #linear model fit (yearly env ~ num)
    fit_temp <- summary(lm(temp ~ num, data = te))
    fit_precip <- summary(lm(precip ~ num, data = te))
    
    #residuals from model
    temp_resid <- residuals(fit_temp)
    precip_resid <- residuals(fit_precip)
    
    #spectral analysis using Lomb-Scargle Periodogram
    #between freq 2/(n*dt) and 1/(2*dt), where dt = 1 and n = 72; 0.0278 to 0.5
    #following Marshall and Burgess 2015 Eco Letters
    #see also Vasseur and Yodzis 2004 Ecology
    #period = 1/freq; ~2 - 36 years
    ll_freq <- 0.0278
    ul_freq <- 0.5
    temp_spec <- lomb::lsp(temp_resid, 
                           from = ll_freq, to = ul_freq, 
                           type = 'frequency',
                           normalize =  'standard', 
                           plot = FALSE)
    #spectral exponent (1/f^beta)
    temp_spec_fit <- summary(lm(log10(temp_spec$power) ~ log10(temp_spec$scanned)))$coefficients[,1]
    #precip - only run if residuals exist (i.e., all precip values were not 0)
    if (sum(precip_resid) > 0)
    {
      precip_spec <- lomb::lsp(precip_resid, 
                               from = ll_freq, to = ul_freq, 
                               type = 'frequency', 
                               normalize = 'standard', 
                               plot = FALSE)
      precip_spec_fit <- summary(lm(log10(precip_spec$power) ~ log10(precip_spec$scanned)))$coefficients[,1]
    } else {
      precip_spec_fit <- rep(NA, 2)
    }
    
    #fill df
    t_df$lon[counter] <- te$lon[1]
    t_df$lat[counter] <- te$lat[1]
    t_df$spectral_beta_temp[counter] <- temp_spec_fit[2] * -1
    t_df$spectral_beta_precip[counter] <- precip_spec_fit[2] * -1
    
    counter <- counter + 1
  }
  if (chunk_counter == 1)
  {
    main_out <- data.frame(cell_id = rep(NA, NROW(t_df) * 20),
                           lon = NA,
                           lat = NA,
                           spectral_beta_temp = NA,
                           spectral_beta_precip = NA)
  }
  
  main_out[counter2:(counter2 + NROW(t_df) - 1),] <- t_df
  
  #advance counters
  chunk_counter <- chunk_counter + 1
  counter2 <- counter2 + NROW(t_df)
  
  #advance lon marker
  lon_start <- lon_start + lenlon
}


# write to csv ------------------------------------------------------------

write.csv(main_out, paste0(args[2], 'Env-spectral-exp-monthly-', args[3], '.csv'), 
          row.names = FALSE)
