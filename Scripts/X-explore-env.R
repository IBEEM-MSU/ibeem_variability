################
# Explore env variability metrics
#
################


# Specify dir --------------------------------------------------

#path CY machine
dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
# devtools::install_github("wmurphyrd/colorplaner")
library(colorplaner)


# read in data -------------------------------------------------

# Climate data (takes a couple of minutes to load)
# only valid cells (land and N of -60S lat)
env.dat <- read.csv(paste0(dir, 'data/L2/climate/era5/Env-main.csv')) %>%
  #only 'valid' cells (those over land and > -60S lat)
  dplyr::filter(valid == TRUE) %>%
  dplyr::select(-valid)

#rasterize env data for extraction (so ensure that ranges that don't intersect cell centers are captured)
env.dat.rast <- dplyr::select(env.dat, lon, lat,
                              grep('temp', colnames(env.dat), value = TRUE),
                              grep('precip', colnames(env.dat), value = TRUE),
                              grep('dhi', colnames(env.dat), value = TRUE)) %>%
  terra::rast(crs = "epsg:4326")


# bivariate map (inter- and intra-annual var) -----------------------------------

#adapted from: https://stackoverflow.com/questions/48572744/plot-a-bivariate-map-in-r

ty_q <- quantile(env.dat$temp_sd_year, seq(0, 1, by = 0.05))
ts_q <- quantile(env.dat$temp_sd_season, seq(0, 1, by = 0.05))
py_q <- quantile(env.dat$precip_cv_year, seq(0, 1, by = 0.05))
ps_q <- quantile(env.dat$precip_cv_season, seq(0, 1, by = 0.05))

ggplot(data = env.dat, 
       aes(lon, lat, 
           fill = temp_sd_year, 
           fill2 = temp_sd_season))+
  geom_tile() +
  scale_fill_colourplane(name = "", 
                         na.color = "white",
                         color_projection = "interpolate", 
                         vertical_color = "#FAE30C",
                         horizontal_color = "#0E91BE", 
                         zero_color = "#E8E6F2",
                         limits_y = c(0, ts_q[20]),
                         limits = c(0, ty_q[20])) +
  theme_minimal()


ggplot(data = env.dat, 
       aes(lon, lat, 
           fill = precip_cv_year, 
           fill2 = precip_cv_season))+
  geom_tile() +
  scale_fill_colourplane(name = "", 
                         na.color = "white",
                         color_projection = "interpolate", 
                         vertical_color = "red",
                         horizontal_color = "blue", 
                         zero_color = "#E8E6F2",
                         limits_y = c(0, ps_q[20]),
                         limits = c(0, py_q[20])) +
  theme_minimal()


# py_q <- quantile(env.dat$precip_cv_year, seq(0, 1, by = 0.05))
# ps_q <- quantile(env.dat$precip_cv_season, seq(0, 1, by = 0.05))
# ggplot(data = env.dat2, 
#        aes(lon, lat, 
#            fill = temp_sp_color_year, 
#            fill2 = temp_rel_slope))+
#   geom_tile() +
#   scale_fill_colourplane(name = "", 
#                          na.color = "white",
#                          color_projection = "interpolate", 
#                          vertical_color = "red",
#                          horizontal_color = "blue", 
#                          zero_color = "#E8E6F2"),
#                          limits_y = c(0, ps_q[20]),
#                          limits = c(0, py_q[20])) +
#   theme_minimal()


# plots -------------------------------------------------------------------

#temperature
plot(env.dat.rast[['temp_mean']], main = 'Temp mean')
plot(env.dat.rast[['temp_sd_year']], main = 'Temp inter-annual sd')
plot(env.dat.rast[['temp_sd_season']], main = 'Temp intra-annual sd')
plot(env.dat.rast[['temp_sp_color_month']], main = 'Temp color (month)')

plot(env.dat.rast[['temp_slope']], main = 'Temp slope')
plot(env.dat.rast[['temp_rel_slope']], main = 'Temp rel slope')

plot(env.dat.rast[['temp_skew']], main = 'Temp skew')
plot(env.dat.rast[['temp_kurt']], main = 'Temp kurtosis') #normal dist kurtosis is 3


#precip
plot(env.dat.rast[['precip_mean']], main = 'Precip mean')
plot(env.dat.rast[['precip_cv_year']], main = 'Precip inter-annual CV')
plot(env.dat.rast[['precip_cv_season']], main = 'Precip intra-annual CV')
plot(env.dat.rast[['precip_sp_color_month']], main = 'Precip color (month)')

plot(env.dat.rast[['precip_slope']], main = 'Precip slope')
plot(env.dat.rast[['precip_rel_slope']], main = 'Precip rel slope')

plot(env.dat.rast[['precip_skew']], main = 'Precip skew')
plot(env.dat.rast[['precip_kurt']], main = 'Precip kurtosis')


# stats -------------------------------------------------------------------

# #median over globe
terra::global(env.dat.rast[['temp_kurt']], fun = function(x) median(x, na.rm = TRUE))
terra::global(env.dat.rast[['temp_skew']], fun = function(x) median(x, na.rm = TRUE))
terra::global(env.dat.rast[['temp_sp_color_month']], fun = function(x) median(x, na.rm = TRUE))



#which areas are highly predictable (P; high temporal autocorrelation) on short time scales (S) and have low intrinsic variability (IV)?
#high P, short S, low IV = short LH
#low P, long S, high IV = long LH

#high P/short S = faster LH
#high IV = slower LH

#is P at lag = gen time the same across all species?

