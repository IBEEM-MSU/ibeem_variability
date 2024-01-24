#----------------------------#
# Pace of life map
#----------------------------#

### Load packages ----

library(tidyverse)
library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)

### Data ----

# Bird data

bird_ras <- terra::rast(paste0(dir, 'data/L3/raster-LH-nsp.tif'))
bird_ras <- terra::rast("~/Documents/Documents/ibeem/raster-gl-dh-nsp.tif") #needs update!
bird_ras <- tt
#bird_ras2 <- bird_ras[[c('median_gl', 'sd_gl', 
#                         'median_dh', 'sd_dh')]]

# Mask areas with fewer than 5 species
msk <- ifel(bird_ras[['n_sp']] < 10, NA, 1)
bird_ras2 <- terra::mask(bird_ras, msk, 
                         inverse = FALSE)

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

### Transform data to Robinson projection ----

# Transform raster to Robinson ("ESRI:54030")
bird_ras_proj <- terra::project(bird_ras2, "ESRI:54030")
plot(bird_ras_proj$median_gl)

# Transform land vector to Robinson
land_50m_robinson <- st_transform(st_wrap_dateline(land_50m), "ESRI:54030")

# Transform bb to Robinson
bb_robinson <- st_transform(bb, crs= "ESRI:54030", type = "robin")

### Median gl map ----

# Make data frame with projected env data 
bird_df <- terra::as.data.frame(bird_ras_proj, xy = T, cells = F, na.rm = T)
head(bird_df)

medgl_q <- quantile(bird_df$median_gl, seq(0, 1, by = 0.05))
n_sp <- quantile(bird_df$n_sp, seq(0, 1, by = 0.05))

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
  # Median gl
  geom_tile(data = bird_df, 
            aes(x = x,
                y = y, 
                fill = median_gl)) +
  scale_fill_gradient(low = "#fef0d9",
                      high = "#b30000",
                      na.value = "#FAFAFA",
                      limits = c(0,5)) +
  # remove X and Y labels from map   
theme(axis.title.x = element_blank(),
      axis.title.y = element_blank())
            