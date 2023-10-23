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

### Transform data to Robinson projection ----

# Transform temp raster to Robinson ("ESRI:54030")
temp_rast_proj <- terra::project(temp_rast, "ESRI:54030")
plot(temp_rast_proj$temp_sd_season)

# Transform land vector to Robinson
land_50m_robinson <- st_transform(st_wrap_dateline(land_50m), "ESRI:54030")

# Transform bb to Robinson
bb_robinson <- st_transform(bb, crs= "ESRI:54030", type = "robin")

### Map ----

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



