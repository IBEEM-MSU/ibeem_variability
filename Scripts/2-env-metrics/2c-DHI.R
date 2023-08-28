################
# #process DHI - create tif with cv_year and cv_season
#
# args:
# - in dir
# - out dir
################

# get args fed to script --------------------------------------------------

# args <- commandArgs(trailingOnly = TRUE)
args <- rep(NA, 2)
# args[1] <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/data/L0/DHI/'
# args[2] <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/data/L1/DHI/'
args[1] <- '/mnt/research/ibeem/variability/data/L0/DHI'
args[2] <- '/mnt/research/ibeem/variability/data/L1/DHI'


# load packages -----------------------------------------------------------

library(tidyverse)
library(terra)


# read in data -------------------------------------------------

#list files
#source: https://silvis.forest.wisc.edu/data/dhis/
#Radeloff et al. 2019 Remote Sensing of Environment
#2003-2015
#1-km resolution
lf <- list.files(args[1], full.names = TRUE)

#read in all tifs
main_rast <- terra::rast(lf)


# process data ------------------------------------------------------------

#indices for each metric
idx <- 1:(length(lf) * 3)
cum_idx <- idx[seq(1, max(idx), 3)]
min_idx <- idx[seq(2, max(idx), 3)]
sea_idx <- idx[seq(3, max(idx), 3)]

#separate rasters
cum_rast <- main_rast[[cum_idx]]
min_rast <- main_rast[[min_idx]]
sea_rast <- main_rast[[sea_idx]]

#CV of cumulative NDVI
sd_cum <- terra::app(cum_rast, fun = sd)
mn_cum <- terra::app(cum_rast, fun = mean)
cv_cum <- sd_cum / mn_cum
#mean of seasonality (CV seasonal NDVI)
mn_sea <- terra::app(sea_rast, fun = mean)

#combine
dhi_met <- c(mn_cum, cv_cum, mn_sea)

#change names (dhi_cum_mean, dhi_cv_year, dhi_cv_season)
names(dhi_met) <- c('dhi_cum_mean', 'dhi_cv_year', 'dhi_cv_season')

#plot
#plot(dhi_met)

#visualize
# tt <- terra::values(dhi_met)
# pdf('~/Desktop/test.pdf')
# plot(tt[,2], tt[,3],
#      xlim = c(0,4),
#      ylim = c(0,4),
#      col = rgb(0,0,0,0.05),
#      pch = 19,
#      xlab = 'CV year',
#      ylab = 'CV season')
# abline(0, 1, col = 'red', lty = 2, lwd = 2)
# dev.off()


#write to file
terra::writeRaster(dhi_met, paste0(args[2], 'DHI_rast.tif'), overwrite = TRUE)
