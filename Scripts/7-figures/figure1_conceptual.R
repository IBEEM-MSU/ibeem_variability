# TITLE:            Figure 1: Conceptual figure  
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Env data (2b), model output (5), ranges shpfiles for example spp
# DATA OUTPUT:      Mean temp/precip for example spp (CSV), sub-plots for Figure 1
# DATE:             January 2024 
# OVERVIEW:         Individual sub-plots were aggregated in Adobe 
# NOTE:             Colorplaner runs on R version 4.1.2 (but not R/4.2.1). Ncdf4 package needs R version 4.2.1


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


rm(list = ls())

# load environmental variables ------------------------------------------------

source("./Scripts/0-config.R")

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

### Save output temp and precip raster for archiving --- 

archive_rast <- c(temp_rast$temp_sd_season, 
                  temp_rast$temp_sd_year, 
                  precip_rast$precip_cv_season, 
                  precip_rast$precip_cv_year)

terra::writeRaster(archive_rast,
                   filename = paste0(dir, 'data/L3/env_var_ras.tif'),
                   overwrite = TRUE)

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

temp_df_map <- temp_df %>%
  select(x, y, temp_sd_season, temp_sd_year) %>%
  mutate(temp_s = ifelse(temp_sd_season > ts_q[20], 16.7, temp_sd_season),
         temp_y = ifelse(temp_sd_year > ty_q[20], 1.26, temp_sd_year))
head(temp_df_map)


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
  geom_tile(data = temp_df_map, 
            aes(x = x,
                y = y, 
                fill = temp_y,
                fill2 = temp_s)) +
  # remove X and Y labels from map (it also removes legend labels :S )  
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  # bivariate color scale 
  scale_fill_colourplane(name = "", 
                         na.color = "#FAFAFA",
                         color_projection = "interpolate", 
                         # zero_color = "#E8E8E8",
                         # vertical_color = "#FAE30C",
                         # horizontal_color = "#0E91BE", 
                         zero_color = "#E8E8E8",
                         vertical_color = "#fff700",
                         horizontal_color = "#ff002f",
                         # use temp quantiles to define limits and breaks
                         limits = c(0, max(temp_df_map$temp_y)),
                         limits_y = c(0, max(temp_df_map$temp_s)))

                         # breaks = c(0,0.5,1.0,1.5,2),
                         # breaks_y = c(5,10,15,20)) 

### Precip variability map ----

# Make data frame with projected env data 
precip_df <- terra::as.data.frame(precip_rast_proj, xy = T, cells = F, na.rm = T)
head(precip_df)
summary(precip_df)

precip_df %>%
  filter(precip_cv_season == Inf)
# 5613 Inf values

precip_df <- precip_df %>%
  mutate(precip_cv_season = ifelse(precip_cv_season == Inf, 999, precip_cv_season)) 

summary(precip_df)

# Get precip quantiles to set vertical and horizontal color scales
py_q <- quantile(precip_df$precip_cv_year, seq(0, 1, by = 0.05))
ps_q <- quantile(precip_df$precip_cv_season, seq(0, 1, by = 0.05), na.rm = T)

precip_df_map <- precip_df %>%
  select(x, y, precip_cv_year, precip_cv_season) %>%
  mutate(precip_s = ifelse(precip_cv_season > ps_q[20], 1.7, precip_cv_season),
         precip_y = ifelse(precip_cv_year > py_q[20], 0.59, precip_cv_year))
head(precip_df_map)

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
  geom_tile(data = precip_df_map, 
            aes(x = x,
                y = y, 
                fill = precip_y,
                fill2 = precip_s)) +
  # remove X and Y labels from map (it also removes legend labels :S )  
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  # bivariate color scale 
  scale_fill_colourplane(name = "", 
                         na.color = "gray95",
                         color_projection = "interpolate",
                         vertical_color = "#0000FF",
                         horizontal_color = "#FF0000",
                         zero_color = "#E8E8E8",
                         limits = c(0, max(precip_df_map$precip_y)),
                         limits_y = c(0, max(precip_df_map$precip_s)))

## Two species example ----

### Select species ----

# load bird data

bird <- readRDS(paste0(dir, 'Results/bird-gl-phylo-vint-', gl_run_date, 
                            '/bird-gl-phylo-vint-data-', gl_run_date, '.rds'))$pro_data

head(bird)
summary(bird)

# lGl min= 0.3534, max = 3.3035, mean = 1.1671
# temp_sd_year min= 0.1247, max= 1.0807, mean= 0.3816
# temp_sd_season min= 0.191, max= 17.455, mean= 1.971

# High inter annual, low intra annual, long gen length
bird %>%
  select(ID, species, Order, Family, lGL, temp_sd_season, temp_sd_year) %>%
  arrange(-lGL) %>%
  filter(temp_sd_year > 0.39) %>%
  filter(temp_sd_season < 1.9) %>%
  print(n = 30)

# amazona ochrocephala GenLength = 16.5, temp_sd_season = 0.827, temp_sd_year = 0.403

bird %>% 
  dplyr::select(ID, species) %>%
  dplyr::filter(species == "Amazona_ochrocephala")
# ID 9030

# High intra annual, low inter annual, short gen length

# bird_data %>%
#   select(Birdtree_name, Family, GenLength, Modeled_max_longevity, temp_sd_season, temp_sd_year, Min.Latitude, Max.Latitude) %>%
#   arrange(GenLength) %>%
#   #filter(Family == "Trochilidae") %>%
#   #arrange(-Modeled_max_longevity) %>%
#   #filter(temp_sd_year < 0.39) %>%
#   filter(temp_sd_season > 1.49) %>%
#   filter(Min.Latitude > 23) %>%
#   print(n= 30)

bird %>%
  select(ID, species, Order, Family, lGL, temp_sd_season, temp_sd_year) %>%
  arrange(lGL) %>%
  arrange(-temp_sd_season) %>%
  filter(lGL < 0.4) %>%
  #filter(temp_sd_year < 0.39) %>%
  #filter(temp_sd_season > 1.49) %>%
  #filter(temp_sd_season > 7) %>%
  #filter(temp_sd_year < 0.39) %>%
  #filter(Family == "Elachuridae") %>%
  print(n= 50)

bird %>%
  select(ID, species, Order, Family, lGL, temp_sd_season, temp_sd_year) %>%
  arrange(lGL) %>%
  #filter(temp_sd_season > 1.49) %>%
  filter(temp_sd_season > 4) %>%
  #filter(temp_sd_year < 0.39) %>%
  #filter(Family == "Elachuridae") %>%
  filter(Family == "Trochilidae") %>%
  print(n= 50)

# Selasphorus platycercus (ID 9708), GenLength = 2.593, temp_sd_season = 8.40, temp_sd_year = 0.639

### Extract species distribution ranges and environmental data ----

# distribution ranges in 'L1/range/'
# Read shapefiles and convert to SpatVector

dr_9030 <- sf::st_read(paste0(dir, "/data/L1/range/bird-breeding/birdtree-9030-breeding.shp")) %>% 
  terra::vect() # parrot

dr_9708 <- sf::st_read(paste0(dir, "/data/L1/range/bird-breeding/birdtree-9708-breeding.shp")) %>% 
  terra::vect() 

# Read temporal env data

ncfiles <- list.files(path = paste0(dir, "/data/L0/climate/era5/T2m"), pattern = "moda", full.names = T)

# Make empty file
get_temp_precip <- function(bird_range){
  temp_mean_sd <- data.frame(matrix(NA, nrow = length(ncfiles), ncol = 25))
  colnames(temp_mean_sd) <- c("mean_temp_01",
                              "mean_temp_02",
                              "mean_temp_03",
                              "mean_temp_04",
                              "mean_temp_05",
                              "mean_temp_06",
                              "mean_temp_07",
                              "mean_temp_08",
                              "mean_temp_09",
                              "mean_temp_10",
                              "mean_temp_11",
                              "mean_temp_12",
                              "sd_temp_01",
                              "sd_temp_02",
                              "sd_temp_03",
                              "sd_temp_04",
                              "sd_temp_05",
                              "sd_temp_06",
                              "sd_temp_07",
                              "sd_temp_08",
                              "sd_temp_09",
                              "sd_temp_10",
                              "sd_temp_11",
                              "sd_temp_12",
                              "year")
  rn <- 1
  
  for (i in 1:length(ncfiles))
  {
    
    ncin <- ncfiles[i]
    print(paste0("Processing ", i, " out of ", length(ncfiles)))
    
    # Open netCDF file
    nc_temp <- nc_open(ncin)
    
    # get lon, lat, and time
    lon <- ncdf4::ncvar_get(nc_temp,"longitude")
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
    current_year <- tmp_df$year[1]
    
    # Make raster
    
    tmp_df_w <- tmp_df %>%
      tidyr::pivot_wider(names_from = month,
                         values_from = temp) %>%
      select(-year)
    
    temp_rast <- terra::rast(tmp_df_w, crs = "EPSG:4326")
    
    # Extract mean temp across range per month
    mean_temp <- terra::extract(temp_rast,
                               bird_range,
                               touches = TRUE,
                               fun = mean,
                               ID = F) 
    
    sd_temp <- terra::extract(temp_rast,
                              bird_range,
                              touches = TRUE,
                              fun = sd,
                              ID = F) %>%
      mutate(year = current_year)
    
    temp_mean_sd[rn,] <- c(mean_temp, sd_temp)
    
    rn <- rn + 1
  }  
  return(temp_mean_sd)
}

# Save files:

temp_mean_sd_sp9030 <- get_temp_precip(dr_9030)
write_csv(temp_mean_sd_sp9030, paste0(dir, "data/L3/temp_mean_sd_sp9030.csv"))  

temp_mean_sd_sp9708 <- get_temp_precip(dr_9708)
write_csv(temp_mean_sd_sp9708, paste0(dir, "data/L3/temp_mean_sd_sp9708.csv"))

### Plot time series ----

#Read in temp_mean_sd files for both spp
#parrot
monthly_temp_sp9030<- read_csv(paste0(dir, "data/L3/temp_mean_sd_sp9030.csv"))

# Broad-tailed hummingbird
monthly_temp_sp9708<- read_csv(paste0(dir, "data/L3/temp_mean_sd_sp9708.csv"))

#pivot longer
monthly_temp_sp9030_long_t <- pivot_longer(monthly_temp_sp9030,
                                         cols = starts_with("mean"),
                                         names_to = "month_temp",
                                         values_to = "temp") %>% 
  select(year, month_temp, temp) %>%
  mutate(month = as.character(substr(month_temp, 11, 12)),
         year = as.character(year)) %>%
  arrange(year, month) %>%
  unite("date_ym", c(year, month), sep = "-", remove = F) %>%
  mutate(date = lubridate::ym(date_ym)) %>%
  select(-month_temp, -date_ym) %>%
  relocate(temp, .after = date)

monthly_temp_sp9030_long_sd <- pivot_longer(monthly_temp_sp9030,
                                           cols = starts_with("sd"),
                                           names_to = "month_sd",
                                           values_to = "sd") %>% 
  select(year, month_sd, sd) %>%
  mutate(month = as.character(substr(month_sd, 9, 10)),
         year = as.character(year)) %>%
  arrange(year, month) %>%
  unite("date_ym", c(year, month), sep = "-", remove = F) %>%
  mutate(date = lubridate::ym(date_ym)) %>%
  select(-month_sd, -date_ym) %>%
  relocate(sd, .after = date)

monthly_temp_sp9030_long <- monthly_temp_sp9030_long_t %>%
  left_join(monthly_temp_sp9030_long_sd, by = c("year", "month", "date"))

rm(monthly_temp_sp9030_long_sd, monthly_temp_sp9030_long_t)

#pivot longer
monthly_temp_sp9708_long_t <- pivot_longer(monthly_temp_sp9708,
                                           cols = starts_with("mean"),
                                           names_to = "month_temp",
                                           values_to = "temp") %>% 
  select(year, month_temp, temp) %>%
  mutate(month = as.character(substr(month_temp, 11, 12)),
         year = as.character(year)) %>%
  arrange(year, month) %>%
  unite("date_ym", c(year, month), sep = "-", remove = F) %>%
  mutate(date = lubridate::ym(date_ym)) %>%
  select(-month_temp, -date_ym) %>%
  relocate(temp, .after = date)

monthly_temp_sp9708_long_sd <- pivot_longer(monthly_temp_sp9708,
                                            cols = starts_with("sd"),
                                            names_to = "month_sd",
                                            values_to = "sd") %>% 
  select(year, month_sd, sd) %>%
  mutate(month = as.character(substr(month_sd, 9, 10)),
         year = as.character(year)) %>%
  arrange(year, month) %>%
  unite("date_ym", c(year, month), sep = "-", remove = F) %>%
  mutate(date = lubridate::ym(date_ym)) %>%
  select(-month_sd, -date_ym) %>%
  relocate(sd, .after = date)

monthly_temp_sp9708_long <- monthly_temp_sp9708_long_t %>%
  left_join(monthly_temp_sp9708_long_sd, by = c("year", "month", "date"))

rm(monthly_temp_sp9708_long_sd, monthly_temp_sp9708_long_t)

# plots

monthly_temp_sp9030_long %>%
  ggplot(aes(x = date, y = temp)) +
  #geom_point() +
  geom_line() +
  geom_ribbon(aes(x = date,
                  ymin = temp - sd,
                  ymax = temp + sd),
              alpha = 0.2) +
  theme_classic() +
  scale_x_date(date_breaks = "5 years",
               date_labels = "%Y") +
  scale_y_continuous(limits = c(-30, 30))

monthly_temp_sp9708_long %>%
  ggplot(aes(x = date, y = temp)) +
  #geom_point() +
  geom_line() +
  geom_ribbon(aes(x = date,
                  ymin = temp - sd,
                  ymax = temp + sd),
              alpha = 0.2) +
  theme_classic() +
  scale_x_date(date_breaks = "5 years",
               date_labels = "%Y") +
  scale_y_continuous(limits = c(-30, 30))

ggplot() +
  geom_line(aes(x = monthly_temp_sp9708_long$date,
                y = monthly_temp_sp9708_long$temp),
            color = "#bf1da1") +
  geom_ribbon(aes(x = monthly_temp_sp9708_long$date,
                  ymin = monthly_temp_sp9708_long$temp - monthly_temp_sp9708_long$sd,
                  ymax = monthly_temp_sp9708_long$temp + monthly_temp_sp9708_long$sd),
              alpha = 0.2,
              fill = "#bf1da1") +
  geom_line(aes(x = monthly_temp_sp9030_long$date,
                y = monthly_temp_sp9030_long$temp),
            color = "#1dbf76") +
  geom_ribbon(aes(x = monthly_temp_sp9030_long$date,
                  ymin = monthly_temp_sp9030_long$temp - monthly_temp_sp9030_long$sd,
                  ymax = monthly_temp_sp9030_long$temp + monthly_temp_sp9030_long$sd),
              alpha = 0.2,
              fill = "#1dbf76") +
  geom_segment(aes(x = as_date(-7490),
                   xend = as_date(-6541),
                   y = -11,
                   yend = -11),
               color = "#bf1da1",
               #fill = "black",
               linewidth = 1.5) +
  geom_segment(aes(x = as_date(-7490),
                   xend = as_date(-1468),
                   y = 30,
                   yend = 30),
               color = "#1dbf76",
               linewidth = 1.5) +
  theme_classic() +
  scale_x_date(date_breaks = "5 years",
               date_labels = "%Y") +
  scale_y_continuous(limits = c(-12, 30)) +
  labs(y = "Temperature (°C)",
       x = NULL) +
  theme(legend.text = element_text(size = 16))

### Plot distribution range map ----

ggplot() +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "grey75",
               color = "gray10",
               size = 0.1) +
  geom_sf(data = dr_9030,
          color = "#1dbf76",
          fill = "#1dbf76",
          linetype = "solid",
          alpha = 0.8) +
  coord_sf(crs = 4326,
           xlim = c(-90,-30),
           ylim = c(-20,20)) +
  theme_bw() +
  labs(y = NULL,
       x = NULL)

ggplot() +
  geom_polygon(data = map_data("world"),
               aes(x = long, y = lat, group = group),
               fill = "grey75",
               color = "gray10",
               size = 0.1) +
  geom_sf(data = dr_9708,
          color = "#bf1da1",
          fill = "#bf1da1",
          linetype = "solid",
          alpha = 0.8) +
  coord_sf(crs = 4326,
           xlim = c(-120,-80),
           ylim = c(15,43)) +
  theme_bw() +
  labs(y = NULL,
       x = NULL)

### Get silhouette ----

# I downloaded silhouettes from phylopic
