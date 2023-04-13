################################################################################
# TITLE: NetCDF Processing Example Code
#
# PURPOSE: This script provides example code for opening, aggregating, and 
# processing netcdf files in R. 
#
# AUTHORS: Kelly Kapsar, 
# CREATED: 2023-04-13
# LAST UPDATED ON: 
#
################################################################################


################ EXAMPLE 1: WIND SPEED DATA ################ 
# File is located on IBEEM google drive > L0 > example_netcdf_files
# Original data source: https://resources.marine.copernicus.eu/?option=com_csw&view=details&product_id=WIND_GLO_WIND_L4_NRT_OBSERVATIONS_012_004

# Open netcdf file 
wind <- nc_open("..//CERSAT-GLO-BLENDED_WIND_L4-V6-OBS_FULL_TIME_SERIE_1626911972037.nc")
# Save metadata to a text file
{
  sink('../CERSAT-GLO-BLENDED_WIND_L4-V6-OBS_FULL_TIME_SERIE_16269119720374.txt')
  print(wind)
  sink()
}

# Read lat lon and time for each observation
lon <- ncvar_get(wind, "lon")
lat <- ncvar_get(wind, "lat", verbose = F)
t <- ncvar_get(wind, "time")

head(lon)

# Read in data from the wind variable and verify the dimensions of the array
wind.array <- ncvar_get(wind, "wind_speed") # 3dim array
dim(wind.array)

# Identify fill value and replace with NA
fillvalue <- ncatt_get(wind, "wind_speed","_FillValue")
fillvalue

wind.array[wind.array == fillvalue$value] <- NA

# Close netcdf file
nc_close(wind)

# Isolate and plot a random time step to check
wind.slice <- wind.array[,,1]

dim(wind.slice) #2dim

wind.r <- raster(t(wind.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat),ymx=max(lat),
                 # Found projection on the website
                 crs=CRS("+proj=longlat +datum=WGS84 +no_defs")) %>%
  flip(direction="y")%>%
  raster::projectRaster(crs=prj)

wind.df <- as.data.frame(wind.r, xy=TRUE) %>% drop_na()
colnames(wind.df) <- c("x","y","windspeed")

plot(wind.r)


# Map of study area with wind data
ggplot() +
  geom_sf(data=basemap.crop, fill="gray", color="black", lwd=0.5) +
  geom_raster(data=wind.df, aes(x=x, y=y, fill=windspeed), alpha=0.9) +
  scale_fill_gradient2() +
  geom_sf(data=study, fill=NA, color="red")


test <- aperm(wind.array, c(2,1,3))
# Make a raster brick of all values 
wind_brick <- brick(test, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
                    crs=CRS("+proj=longlat +datum=WGS84 +no_defs")) %>% 
  flip(direction="y") %>%
  raster::projectRaster(crs=prj) 

# Convert date from seconds since 01/01/1970 to yyyy-mm-dd format
t2 <- as.POSIXct("1900-01-01 00:00")+as.difftime(t,units="hours")
t2 <- format(t2, "%G-W%V")

# Name raster layers after the date that they portray
names(wind_brick) <-t2
wind_brick <- raster::setZ(wind_brick, t2)

# Save output file
# writeRaster(wind_brick, "../Data_Processed/wind_AOOS_cropped_4336.tif")

# Fun animation of the raster brick
# animate(wind_brick, pause=0.5, n=1)

# Calculate mean monthly wind rasters 
wind_week <- zApply(wind_brick, t2, fun = mean)

# animate(wind_week, pause=0.5, n=1)

# plot monthly average wind data for November, 2018
windweek.df <- as.data.frame(wind_week[[1]], xy=TRUE)  %>% drop_na()
colnames(windweek.df) <- c("x","y","windspeed")

ggplot() +
  geom_sf(data=basemap.crop, fill="gray", color="black", lwd=0.5) +
  geom_raster(data=windweek.df, aes(x=x, y=y, fill=windspeed), alpha=0.9) +
  scale_fill_gradient2() +
  geom_sf(data=study, fill=NA, color="red")

# animate(wind_week, pause=0.5, n=1)

writeRaster(wind_week, "../Data_Processed/wind_weekly.tif", options="INTERLEAVE=BAND", overwrite=T)
saveRDS(wind_week, "../Data_Processed/wind_weekly.rds")
wind <- wind_week
save(wind, file="../Data_Processed/wind_weekly.rda")

# Save output file
# writeRaster(wind_week, "../Data_Processed/wind_weekly.nc",
# overwrite=TRUE, format="CDF",
# varname="wind_speed", varunit="m/s",
# longname="Wind Speed -- raster brick to netCDF",
# xname="lon",   yname="lat",zname="time",
# zunit="numeric")