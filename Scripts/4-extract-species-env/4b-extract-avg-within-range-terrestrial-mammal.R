# 4-example-species-extract.R: a script to start thinking about how we'll want
#                              to calculate metrics within each species range
#                              and output all the associated information for 
#                              running models of the form trait ~ climate

# Hello 
rm(list = ls())
library(tidyverse)
library(sf)

# Specify directory -------------------------------------------------------
MAM.dir <- './data/L1/range-mammal/'
env.dir <- './data/L2/climate/era5/'
out.dir <- './data/L2/range-env-pieces-mammal/'


# Get the current file to process -----------------------------------------
file.name <- commandArgs(trailingOnly = TRUE)
# Testing
# file.name <- 'BLIDsPiece-14.rda'
if(length(file.name) == 0) base::stop('Need to give the file name to process')

# Read in data ------------------------------------------------------------
# Loads in current set of ids (a vector called ids)
load(paste0(MAM.dir, file.name))
# Climate data (takes a couple of minutes to load
# NOTE: can change this to the GAM stuff if that's what we want.
env.dat <- read.csv(paste0(env.dir, 'Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv'))

##MIGHT NEED TO RASTERIZE AND THEN EXTRACT INSTEAD, TO MAKE SURE CELLS THAT JUST INTERSECT (BUT DON'T NECESSARILY HAVE CENTROIDS THAT OVERLAP RANGE) GET EXTRACTED
# VV I (Casey) think this should rasterize
#I think this needs to be done separately for temp and precip - maybe not though
# temp.dat.rast <- dplyr::filter(env.dat, 
#                               var = 'temp',
#                               valid == TRUE) %>% #valid cells are only those over land
#   terra::rast(crs = "epsg:4326")
# precip.dat.rast <- dplyr::filter(env.dat, 
#                                var = 'precip',
#                                valid == TRUE) %>%
#   terra::rast(crs = "epsg:4326")
#
# VV I think this should extract
# terra::extract(temp.dat.rast, 
#                terra::vect(curr.range),
#                touches = TRUE,
#                fun = mean)

# Convert to sf
env.dat.sf <- st_as_sf(env.dat, 
		       coords = c('lon', 'lat'),
		       crs = 4326)

# Loop through all the species in the current file ------------------------
# Make data frame to hold everything
# Number of climate variables
clim.var.names <- sort(unique(env.dat$var))
n.clim.var <- length(clim.var.names)
# Different metrics of the variables
metrics <- env.dat %>%
  dplyr::select(-cell_id, -var, -lon, -lat) %>%
  colnames()
n.metrics <- length(metrics)
avg.clim.df <- matrix(NA, nrow = length(ids), ncol = n.metrics * n.clim.var + 1)
# Definitely a better way to do this
tmp.names <- paste(rep(metrics, each = n.clim.var), clim.var.names, sep = "_")
colnames(avg.clim.df) <- c('id', tmp.names)
id.indx <- which(colnames(avg.clim.df) == 'id')
avg.clim.df <- as.data.frame(avg.clim.df)

start <- proc.time()			  
for (i in 1:length(ids)) {
  print(paste0("Currently on species ", i, " out of ", length(ids)))
  curr.sp <- ids[i]  
  curr.range <- st_read(paste0(MAM.dir, curr.sp, '.shp'))
  # Get pixels within current species range. Takes about 2 min per species
  env.curr.range <- st_crop(env.dat.sf, curr.range)
  
  ##CALC RANGE SIZE AND LAT/LON CENTROIDS HERE?
  ##RASTERIZE SHP FILE, RESAMPLE TO ENV GRID, REPLACE RASTER VALUES (1'S) WITH GENERATION LENGTH
  ##believe this should resample to env grid
  #terra::resample(range.rast, temp.dat.rast, method = 'near'))
  
  # Get average of all the environmental variables for the given species
  bad.cols <- which(names(env.curr.range) %in% c('cell_id', 'geometry'))
  avg.clim.vars <- env.curr.range %>%
    dplyr::select(-cell_id, -geometry) %>%
    st_drop_geometry() %>%
    group_by(var) %>%
    summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    arrange(var) %>%
    dplyr::select(-var) %>%
    as.matrix() %>%
    c()
  avg.clim.df[i, id.indx] <- curr.sp
  # Need this if, as there are some species whose ranges don't overlap the climate
  # data (e.g., Snares Penguin)
  if (length(avg.clim.vars) > 0) {
    avg.clim.df[i, -id.indx] <- avg.clim.vars
  }
  rm(curr.range, env.curr.range)
  gc()
}
proc.time() - start
save(avg.clim.df, file = paste0(out.dir, 'summarized-data-piece-', 
				str_extract(file.name, '\\d+')))
