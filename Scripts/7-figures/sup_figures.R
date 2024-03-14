# TITLE:            Supplementary Figures   
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Peter Williams, Jeff Dozer, Adriana Uscanga, Lala Kounta, Kelly Kapsar, Phoebe Zarnetske, Pat Bills
# DATA INPUT:       
# DATA OUTPUT:      
# DATE:             January 2024 
# OVERVIEW:         
# Maps ----

# Load packages ----

library(tidyverse)
library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)

# Load raster ----

del_ras <- terra::rast(paste0(dir, 'data/L3/raster-gl-dT-dP-nsp.tif')) 

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

# Transform raster to Robinson ("ESRI:54030")
del_ras_proj <- terra::project(del_ras, "ESRI:54030")

# Transform land vector to Robinson
land_50m_robinson <- st_transform(st_wrap_dateline(land_50m), "ESRI:54030")

# Transform bb to Robinson
bb_robinson <- st_transform(bb, crs = "ESRI:54030", type = "robin")

# Make data frame with projected data 
del_df <- terra::as.data.frame(del_ras_proj, xy = T, cells = F, na.rm = T)
head(del_df)

# Calculate quantiles for color scale
gl_q <- quantile(log(del_df$mean_gl), seq(0, 1, by = 0.05))
n_sp_q <- quantile(del_df$n_sp, seq(0, 1, by = 0.05))

# Base map with ggplot 
base <- ggplot() +
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
  theme_minimal()


# Gen Length
base + geom_tile(data = del_df,
                 aes(x = x,
                     y = y,
                     fill = log(mean_gl))) +
  scale_fill_gradient(name = "Log Gen Length",
                      low = "#f7fbff",
                      high = "#08519c",
                      na.value = "#FAFAFA",
                      limits = range(gl_q)) +
  # remove X and Y labels from map   
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle('Mean generation length')

# Number of sp
base + geom_tile(data = del_df,
                 aes(x = x,
                     y = y,
                     fill = n_sp)) +
  scale_fill_gradient(name = "No. spp",
                      low = "#f7fcf5",
                      high = "#238443",
                      na.value = "#FAFAFA",
                      limits = range(n_sp_q)) +
  # remove X and Y labels from map   
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle('Number of species')
