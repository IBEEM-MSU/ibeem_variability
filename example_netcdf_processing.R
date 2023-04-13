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

# Load libraries 
library(ncdf4)

################ EXAMPLE 1: WIND SPEED DATA ################ 
# File is located on IBEEM google drive > L0 > example_netcdf_files
# Original data source: https://resources.marine.copernicus.eu/?option=com_csw&view=details&product_id=WIND_GLO_WIND_L4_NRT_OBSERVATIONS_012_004

# Open netcdf file 
wind <- nc_open("../CERSAT-GLO-BLENDED_WIND_L4-V6-OBS_FULL_TIME_SERIE_1626911972037.nc")
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

################ EXAMPLE 2: TEMPERATURE DATA (CRU) ################ 

# Example from https://pjbartlein.github.io/REarthSysSci/netCDF.html#introduction 
# The examples make use of a netCDF file of climate data from the Climate Research Unit http://www.cru.uea.ac.uk/data, 
# consisting of long-term mean values (1961-1990) of near-surface air temperature on a 0.5-degree grid (for land points). 
# The dimensions of the array are 720 (longitudes) x 360 (latitudes) x 12 (months), 
# thus forming a raster “stack” or “brick” consisting of 12 layers.
# File is located on IBEEM google drive > L0 > example_netcdf_files

# Load packages

library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)

# set path and filename
ncpath <- "~/example_netcdf_files/"
ncname <- "cru10min30_tmp"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "tmp"  #temperature 

# open a netCDF file
ncin <- nc_open(ncfname)
print(ncin)

# get longitude and latitude
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
head(lon)

# get longitude and latitude
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"lat")
nlat <- dim(lat)
head(lat)

print(c(nlon,nlat))

# get time
time <- ncvar_get(ncin,"time")
time
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt
tunits

# get temperature
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# get global attributes
title <- ncatt_get(ncin,0,"title")
institution <- ncatt_get(ncin,0,"institution")
datasource <- ncatt_get(ncin,0,"source")
references <- ncatt_get(ncin,0,"references")
history <- ncatt_get(ncin,0,"history")
Conventions <- ncatt_get(ncin,0,"Conventions")

# close netCDF file
nc_close(ncfname)

# convert time -- split the time units string into fields
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(unlist(tdstr)[3])
tyear <- as.integer(unlist(tdstr)[1])
chron(time,origin=c(tmonth, tday, tyear))

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
length(na.omit(as.vector(tmp_array[,,1])))

# get a single slice or layer (January)
m <- 1
tmp_slice <- tmp_array[,,m]

# quick map (to see if it's working)
image(lon,lat,tmp_slice, col=rev(brewer.pal(10,"RdBu")))

# create dataframe -- reshape data
# matrix (nlon*nlat rows by 2 cols) of lons and lats
lonlat <- as.matrix(expand.grid(lon,lat))
dim(lonlat)

# vector of `tmp` values
tmp_vec <- as.vector(tmp_slice)
length(tmp_vec)

# create dataframe and add names
tmp_df01 <- data.frame(cbind(lonlat,tmp_vec))
names(tmp_df01) <- c("lon","lat",paste(dname,as.character(m), sep="_"))
head(na.omit(tmp_df01), 10)

# reshape the whole array into vector
tmp_vec_long <- as.vector(tmp_array)
length(tmp_vec_long)

# reshape the vector into a matrix
tmp_mat <- matrix(tmp_vec_long, nrow=nlon*nlat, ncol=nt)
dim(tmp_mat)

head(na.omit(tmp_mat))

# create a dataframe
lonlat <- as.matrix(expand.grid(lon,lat))
tmp_df02 <- data.frame(cbind(lonlat,tmp_mat))
names(tmp_df02) <- c("lon","lat","tmpJan","tmpFeb","tmpMar","tmpApr","tmpMay","tmpJun",
                     "tmpJul","tmpAug","tmpSep","tmpOct","tmpNov","tmpDec")
# options(width=96)
head(na.omit(tmp_df02, 20))

# get the annual mean and mean temp of the warmest and coldest month (MTWA and MTCO)
tmp_df02$mtwa <- apply(tmp_df02[3:14],1,max) # mtwa
tmp_df02$mtco <- apply(tmp_df02[3:14],1,min) # mtco
tmp_df02$mat <- apply(tmp_df02[3:14],1,mean) # annual (i.e. row) means
head(na.omit(tmp_df02))

dim(na.omit(tmp_df02))

# There's more on converting data frames to arrays and creating netCDF file in the website (https://pjbartlein.github.io/REarthSysSci/netCDF.html#introduction)
# There's also a module on time-series plots: https://pjbartlein.github.io/REarthSysSci/rasterVis01.html#time-series-plots
# Also, a handy function to track processing time:

# time an R process
ptm <- proc.time() # start the timer
# ... some code ...
proc.time() - ptm # how long?