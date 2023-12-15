#---------------------------------#
# Figure 1. Conceptual figure ----
#---------------------------------#

## Environmental variability map ----

#adapted from: https://stackoverflow.com/questions/48572744/plot-a-bivariate-map-in-r
# and https://bluegreenlabs.org/post/map-building-3/

### Load packages ----

library(tidyverse)
library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
# devtools::install_github("wmurphyrd/colorplaner")
library(colorplaner)
library(ncdf4)

### Data ----

# Environmental data 

env_data <- read.csv(paste0(dir, 'data/L2/climate/era5/Env-main.csv')) %>%
  #only 'valid' cells (those over land and > -60S lat)
  dplyr::filter(valid == TRUE) %>%
  dplyr::select(-valid)

# Land outline

# download land outline (50 m) from natural earth
land_50m <- ne_download(scale = 50,
                        category = "physical",
                        type = "land",
                        returnclass = c("sf"))

# Map bounding box

bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), 
          crs = st_crs(4326)),
  n = 100))

### Make raster with temp data ----

temp_rast <- dplyr::select(env_data, lon, lat,
                           grep('temp', colnames(env_data), value = TRUE)) %>%
  # use "EPSG:4326" projection
  terra::rast(crs = "EPSG:4326")

# plot example:
plot(temp_rast$temp_sd_season)

### Make raster with precip data ----

precip_rast <- dplyr::select(env_data, lon, lat,
                           grep('precip', colnames(env_data), value = TRUE)) %>%
  # use "EPSG:4326" projection
  terra::rast(crs = "EPSG:4326")

# plot example:
plot(precip_rast$precip_cv_season)

### Transform data to Robinson projection ----

# Transform temp raster to Robinson ("ESRI:54030")
temp_rast_proj <- terra::project(temp_rast, "ESRI:54030")
plot(temp_rast_proj$temp_sd_season)

# Transform precip raster to Robinson ("ESRI:54030")
precip_rast_proj <- terra::project(precip_rast, "ESRI:54030")
plot(precip_rast_proj$precip_cv_season)

# Transform land vector to Robinson
land_50m_robinson <- st_transform(st_wrap_dateline(land_50m), "ESRI:54030")

# Transform bb to Robinson
bb_robinson <- st_transform(bb, crs= "ESRI:54030", type = "robin")

### Temp variability map ----

# Make data frame with projected env data 
temp_df <- terra::as.data.frame(temp_rast_proj, xy = T, cells = F, na.rm = T)
head(temp_df)

# Get temp quantiles to set vertical and horizontal color scales
ty_q <- quantile(temp_df$temp_sd_year, seq(0, 1, by = 0.05))
ts_q <- quantile(temp_df$temp_sd_season, seq(0, 1, by = 0.05))

# Map with ggplot 
ggplot() +
  # Map background
  geom_sf(data = land_50m_robinson,
          color = "gray80",
          fill = "gray95",
          linetype = "solid",
          size = 0.2) +
  geom_sf(data = bb_robinson,
          color = "white",
          fill = NA,
          linetype = "solid") +
  theme_minimal() +
  # Temp variation
  # geom_tile will through warning for fill2--it's ok, ignore.
  geom_tile(data = temp_df, 
            aes(x = x,
                y = y, 
                fill = temp_sd_year,
                fill2 = temp_sd_season)) +
  # remove X and Y labels from map (it also removes legend labels :S )  
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  # bivariate color scale 
  scale_fill_colourplane(name = "", 
                         na.color = "#FAFAFA",
                         color_projection = "interpolate", 
                         # vertical_color = "#FFAA00",
                         # horizontal_color = "#5555FF",
                         # zero_color = "#E8E8E8",
                         # vertical_color = "#FAE30C",
                         # horizontal_color = "#0E91BE",
                         # zero_color = "#E8E8E8",
                         # vertical_color = "#FF0000",
                         # horizontal_color = "#0000FF",
                         # zero_color = "#E8E8E8",
                         vertical_color = "#0071AA",
                         horizontal_color = "#00AA00",
                         zero_color = "#E8E8E8",
                         # use temp quantiles to define limits and breaks
                         limits = c(0, 2),
                         limits_y = c(0, 20),
                         breaks = c(0,0.5,1.0,1.5,2),
                         breaks_y = c(5,10,15,20)) 

### Precip variability map ----

# Make data frame with projected env data 
precip_df <- terra::as.data.frame(precip_rast_proj, xy = T, cells = F, na.rm = T)
head(precip_df)

precip_df %>%
  filter(precip_cv_season == Inf)
# 5613 Inf values

precip_df <- precip_df %>%
  mutate(precip_cv_season = ifelse(precip_cv_season == Inf, NA, precip_cv_season)) 

# Get precip quantiles to set vertical and horizontal color scales
py_q <- quantile(precip_df$precip_cv_year, seq(0, 1, by = 0.05))
ps_q <- quantile(precip_df$precip_cv_season, seq(0, 1, by = 0.05), na.rm = T)

# Map with ggplot 
ggplot() +
  # Map background
  geom_sf(data = land_50m_robinson,
          color = "gray80",
          fill = "gray95",
          linetype = "solid",
          size = 0.2) +
  geom_sf(data = bb_robinson,
          color = "white",
          fill = NA,
          linetype = "solid") +
  theme_minimal() +
  # Temp variation
  # geom_tile will through warning for fill2--it's ok, ignore.
  geom_tile(data = precip_df, 
            aes(x = x,
                y = y, 
                fill = precip_cv_year,
                fill2 = precip_cv_season)) +
  # remove X and Y labels from map (it also removes legend labels :S )  
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  # bivariate color scale 
  scale_fill_colourplane(name = "", 
                         na.color = "#FAFAFA",
                         color_projection = "interpolate",
                         vertical_color = "#FF0000",
                         horizontal_color = "#0000FF",
                         zero_color = "#E8E8E8",
                         limits = c(0, py_q[20]),
                         limits_y = c(0, ps_q[20]))

## Two species example ----

### Select species ----

# load bird data

bird_data <- read_csv("~/Documents/Documents/ibeem/main-bird-data-birdtree2.csv") #needs update

head(bird_data)
summary(bird_data)

# Gen length min = 1.42, max = 27.87, median = 3.041
# Temp SD season min = 0.19, max = 18.91, median = 1.49
# Temp SD year min = 0.124, max = 1.468, median = 0.391

# High inter annual, low intra annual, long gen length
bird_data %>%
  select(Birdtree_name, GenLength, Modeled_max_longevity, temp_sd_season, temp_sd_year, Min.Latitude, Max.Latitude) %>%
  arrange(-GenLength) %>%
  #arrange(-Modeled_max_longevity) %>%
  filter(temp_sd_year > 0.39) %>%
  filter(temp_sd_season < 1.49) %>%
  print(n= 30)

# amazona ochrocephala GenLength = 16.5, temp_sd_season = 0.827, temp_sd_year = 0.403

bird_data %>% 
  dplyr::select(ID, Birdtree_name) %>%
  dplyr::filter(Birdtree_name == "amazona ochrocephala")
# ID 9030

# High intra annual, low inter annual, short gen length

bird_data %>%
  select(Birdtree_name, Family, GenLength, Modeled_max_longevity, temp_sd_season, temp_sd_year, Min.Latitude, Max.Latitude) %>%
  arrange(GenLength) %>%
  #filter(Family == "Trochilidae") %>%
  #arrange(-Modeled_max_longevity) %>%
  #filter(temp_sd_year < 0.39) %>%
  filter(temp_sd_season > 1.49) %>%
  filter(Min.Latitude > 23) %>%
  print(n= 30)

# parus venustulus GenLength = 1.63, temp_sd_season = 8.38, temp_sd_year = 0.466
# stellula calliope (hummingbird)   GenLength = 1.95, temp_sd_season = 8.99, temp_sd_year 0.702

bird_data %>% 
  dplyr::select(ID, Birdtree_name) %>%
  dplyr::filter(Birdtree_name == "parus venustulus")
# ID 2712

bird_data %>% 
  dplyr::select(ID, Birdtree_name) %>%
  dplyr::filter(Birdtree_name == "stellula calliope")
# ID 366

### Extract species distribution ranges and environmental data ----

# distribution ranges in 'L1/range/'
# Read shapefiles
# files are in "/Volumes/home-219/uscanga1/Documents" for now--needs update!

dr_9030 <- sf::st_read("/Volumes/home-219/uscanga1/Documents/bird-breeding/birdtree-9030-breeding.shp")
dr_2712 <- sf::st_read("/Volumes/home-219/uscanga1/Documents/bird-breeding/birdtree-2712-breeding.shp")
dr_366 <- sf::st_read("/Volumes/home-219/uscanga1/Documents/bird-breeding/birdtree-366-breeding.shp")

plot(dr_9030)

# Read temporal env data

ncfiles <- list.files(path = "T2m", pattern = "moda", full.names = T)

for (i in 1:length(ncfiles))
{
  
  ncin <- ncfiles[i]
  print(paste0("Processing ", i, " out of ", length(ncfiles)))
  
  # Open netCDF file
  nc_temp <- nc_open(ncin)
  
  # get lon, lat, and time
  lon <- ncvar_get(nc_temp,"longitude")
  dlon <- dim(lon)
  lat <- ncvar_get(nc_temp,"latitude")
  dlat <- dim(lat)
  time <- ncvar_get(nc_temp,"time")
  
  #convert lon to -180 to 180
  lon[which(lon > 180)] <- (360 - lon[which(lon> 180)]) * -1
  
  #covert dates to ymd
  ymd_dates <- lubridate::ymd("1900-01-01") + lubridate::hours(time)
  
  # get temperature
  dname <- "VAR_2T"
  tmp_array <- ncvar_get(nc_temp, dname)
  
  # fill in NA values
  fillvalue_temp <- ncdf4::ncatt_get(nc_temp, "VAR_2T","_FillValue")
  tmp_array[tmp_array == fillvalue_temp$value] <- NA
  
  # convert from K to C
  tmp_array_c <- tmp_array - 273.15
  rm(tmp_array)
  
  # Close netcdf file
  ncdf4::nc_close(nc_temp)
  
  # convert to 2d vector
  tmp_vec <- round(as.vector(tmp_array_c), 2)
  
  # create a dataframe
  llt <- expand.grid(lon, 
                     lat, 
                     lubridate::year(ymd_dates)[1], 
                     lubridate::month(ymd_dates))
  tmp_df <- data.frame(llt, tmp_vec)
  
  colnames(tmp_df) <- c('lon', 'lat', 'year', 'month', 'temp')
  
  months <- lubridate::month(ymd_dates)
  months_ch <- as.character(months)
  
  # Make raster
  
  tmp_df_w <- tmp_df %>%
    tidyr::pivot_wider(names_from = month,
                       values_from = temp) %>%
    select(-year)
  
  temp_rast <- terra::rast(tmp_df_w, crs = "EPSG:4326")
}  
  # Extract mean temp across range per month
  
  

# open a netCDF file
ncin <- nc_open("T2m/e5p.moda.an.sfc.128_167_2t.ll025sc.1950010100_1950120100.nc")
#print(ncin)

# get longitude and latitude
lon <- ncvar_get(ncin,"longitude")
nlon <- dim(lon)
# head(lon)
# 
lat <- ncvar_get(ncin,"latitude")
nlat <- dim(lat)

head(lat)
tail(lat)

head(lon)
tail(lon)

print(c(nlon,nlat))

#convert lon to -180 to 180
lon[which(lon > 180)] <- (360 - lon[which(lon> 180)]) * -1

# # get time
time <- ncvar_get(ncin,"time")
# time
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
# nt
tunits
#
#covert dates to ymd
ymd_dates <- lubridate::ymd("1900-01-01") + lubridate::hours(time)
lubridate::month(ymd_dates)

# # get temperature
dname <- "VAR_2T"
tmp_array <- ncvar_get(ncin, dname)

fillvalue_temp <- ncdf4::ncatt_get(ncin, "VAR_2T","_FillValue")
tmp_array[tmp_array == fillvalue_temp$value] <- NA

#convert from K to C
head(tmp_array)
tmp_array_c <- tmp_array - 273.15
head(tmp_array_c)
rm(tmp_array)
# Close netcdf file
ncdf4::nc_close(ncin)

# get a single slice or layer (January)
m <- 1
tmp_slice <- tmp_array_c[,,m]

head(tmp_slice)

tmp_vec <- round(as.vector(tmp_array_c), 2)

# create a dataframe
llt <- expand.grid(lon, 
                   lat, 
                   lubridate::year(ymd_dates)[1], 
                   lubridate::month(ymd_dates))
tmp_df <- data.frame(llt, tmp_vec)

colnames(tmp_df) <- c('lon', 'lat', 'year', 'month', 'temp')

# ^Make loop to read in all files 
# In loop, transform annual env data to raster
# and extract info per spp

env_data_years <- env_data %>%
  group_by(year) %>%
  summarize(years = unique(year)) %>%
  select(years) %>%
  as.list()

years <- env_data_years[["years"]]

env_data_raster <- NULL

env_data_sp9030 <- data.frame(matrix(NA, nrow = length(years), ncol = 3))

rn <- 1

for (i in 1:length(years)) {
  
  current_year <- years[i]
  print(paste0("Year ", i, " out of ", length(years)))

  # Make raster per year
  
  env_data_raster <- env_data %>%
    dplyr::select(lon, lat, year, mean_temp) %>%
    dplyr::filter(year == current_year) %>%
    terra::rast(crs = "epsg:4326")
  
  # Extract mean temp and sd across dist range per year 
  
  avg_temp <- terra::extract(env_data_raster,
                             terra::vect(dr_9030),
                             touches = TRUE,
                             fun = function(x) median(x, na.rm = TRUE))
  
  # sd_temp <- terra::extract(env_data_raster,
  #                           terra::vect(dr_2712),
  #                           touches = TRUE,
  #                           fun = function(x) sd(x, na.rm = TRUE))
  
  env_data_sp9030[rn,] <- avg_temp
  
  rn <- rn + 1
  
}

# Transform env_data to raster

env_data_raster_y1 <- env_data %>%
  dplyr::select(lon, lat, year, mean_temp) %>%
  dplyr::filter(year == 1979) %>%
  terra::rast(crs = "epsg:4326")

# average temp over range

avg_temp <- terra::extract(env_data_raster_y1,
                           terra::vect(dr_2712),
                           touches = TRUE,
                           fun = function(x) median(x, na.rm = TRUE))
sd_temp <- terra::extract(env_data_raster_y1,
                          terra::vect(dr_2712),
                          touches = TRUE,
                          fun = function(x) sd(x, na.rm = TRUE))

### Plot time series ----

### Plot distribution range map ----

### Get silhouette 