# TITLE:            Figure 5: Histograms and maps of dT and dP   
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Model output (5), raster stack (6) 
# DATA OUTPUT:      Histograms and maps of dT and dP used in figure 5
# DATE:             January 2024 
# OVERVIEW:         Individual sub-plots were aggregated in Adobe 

rm(list = ls())

# load environmental variables ------------------------------------------------

source("./Scripts/0-config.R")

# Load packages ----

library(tidyverse)
library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)

# Read in data ------------------------------------------------------------

# df from results
bird_df <- readRDS(paste0(dir, 'Results/bird-gl-phylo-vint-berk-oe-', gl_run_date, 
                          '/bird-gl-phylo-vint-berk-oe-data-', gl_run_date, '.rds'))$pro_data

# stacked raster
del_ras <- terra::rast(paste0(dir, 'data/L3/raster-gl-dT-dP-nsp.tif')) 


# Delta histograms ------------------------------------------------------

## Temp delta:

# dt_hist <- hist(bird_df$temp_delta, breaks = 50, plot = F)
# dt_col <- ifelse(dt_hist$breaks >= 0.1 & dt_hist$breaks < 0.3, rgb(0.2, 0.8, 0.5, 0.5), 
#                  ifelse(dt_hist$breaks >= 0.3, 'purple', rgb(0.2,0.2,0.2,0.2)))
# range(bird_df$temp_delta)
# XLIM <- c(-0.2, 1)

# Alternative colors:
dt_hist <- hist(bird_df$temp_delta, breaks = 50, plot = F)
dt_col <- ifelse(dt_hist$breaks >= 0.1 & dt_hist$breaks < 0.3, "gray60", 
                 ifelse(dt_hist$breaks >= 0.3, "gray20", "gray85"))

plot(dt_hist, col = dt_col, border = FALSE, main = NULL,  
     xlab = 'Delta T') 

#      xlim = XLIM,
#      xaxt = 'n')
# axis(1, at = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))

## Precip delta:

# dp_hist <- hist(bird_df$precip_delta, breaks = 50, plot = F)
# dp_col <- ifelse((dp_hist$breaks >= 0.1 & dp_hist$breaks < 0.3) | dp_hist$breaks <= -0.1, rgb(0.2, 0.8, 0.5, 0.5), 
#                  ifelse(dp_hist$breaks >= 0.3, 'purple', rgb(0.2,0.2,0.2,0.2)))

# Alternative colors:
dp_hist <- hist(bird_df$precip_delta, breaks = 50, plot = F)
dp_col <- ifelse((dp_hist$breaks >= 0.1 & dp_hist$breaks < 0.3) | dp_hist$breaks <= -0.1, "gray60", 
                 ifelse(dp_hist$breaks >= 0.3, "gray20", "gray85"))

plot(dp_hist, col = dp_col, border = FALSE, main = NULL, 
     xlab = 'Delta P')

# xlim = XLIM,
# xaxt = 'n')
# range(bird_df$precip_delta)
# axis(1, at = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))


# Delta maps --------------------------------------------------------------

# Mask areas with fewer than NSP species
NSP <- 10
msk <- ifel(del_ras[['n_sp']] < NSP, NA, 1)
del_ras2 <- terra::mask(del_ras, msk, 
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

# Transform raster to Robinson ("ESRI:54030")
del_ras_proj <- terra::project(del_ras2, "ESRI:54030")

# Transform land vector to Robinson
land_50m_robinson <- st_transform(st_wrap_dateline(land_50m), "ESRI:54030")

# Transform bb to Robinson
bb_robinson <- st_transform(bb, crs = "ESRI:54030", type = "robin")

# Make data frame with projected data 
del_df <- terra::as.data.frame(del_ras_proj, xy = T, cells = F, na.rm = T)
head(del_df)

# Calculate quantiles for color scale
med_dT_q <- quantile(del_df$median_dT, seq(0, 1, by = 0.05))
med_dP_q <- quantile(del_df$median_dP, seq(0, 1, by = 0.05))
n_sp <- quantile(del_df$n_sp, seq(0, 1, by = 0.05))

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
  

# delta T
delta_T_map <- base + geom_tile(data = del_df, 
            aes(x = x,
                y = y, 
                fill = median_dT)) +
  scale_fill_gradient(name = "Median delta T",
                      low = "#ffffcc",
                      high = "#bd0026",
                      na.value = "#FAFAFA",
                      #limits = c(0.08, 0.20)) +
                      limits = range(med_dT_q)) +
  # remove X and Y labels from map   
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) #+
  #ggtitle('Delta T')

delta_T_map

# delta P
delta_P_map <- base + geom_tile(data = del_df, 
                 aes(x = x,
                     y = y, 
                     fill = median_dP)) +
  scale_fill_gradient(name = "Median delta P",
                      low = "#f7fcfd",
                      high = "#3f007d",
                      na.value = "#FAFAFA",
                      #limits = c(0,5)) +
                      limits = range(med_dP_q)) +
  # remove X and Y labels from map   
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) #+
  ggtitle('Delta P')

delta_P_map  

# correlation delta T and delta P -----------------------------------------

cor(bird_df$precip_delta, bird_df$temp_delta)

