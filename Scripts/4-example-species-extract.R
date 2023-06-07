# 4-example-species-extract.R: a script to start thinking about how we'll want
#                              to calculate metrics within each species range
#                              and output all the associated information for 
#                              running models of the form trait ~ climate
rm(list = ls())
library(tidyverse)
library(sf)

# Specify directory -------------------------------------------------------
BL.dir <- '/mnt/research/ibeem/data/L1/range/bird-breeding/'
trait.dir <- '/mnt/research/ibeem/data/L1/trait/'
env.dir <- '/mnt/research/ibeem/L2/climate/era5/'


# Read in data ------------------------------------------------------------
# Birdlife species ids (object called BL.unique.ids)
load(paste0(BL.dir, 'BL-ids.rda'))
# Generation time data from Bird et al
gen.dat <- read.csv(paste0(trait.dir, 'bird-et-al-data-with-id.csv'))
# Climate data
env.dat <- read.csv(paste0(env.dir, 'Env-var-1_2_3_4_5_6_7_8_9_10_11_12.csv'))

env.dat.sf <- st_as_sf(env.dat, 
		       coords = c('lon', 'lat'),
		       crs = 4326)

# Some preliminary analyses -----------------------------------------------
curr.sp <- BL.unique.ids[1]
curr.range <- st_read(paste0(BL.dir, curr.sp, '-breeding.shp'))
# Get pixels within current species range. 
# Takes a couple minutes for just 1 species. Do we want to save this out 
# somewhere?
env.curr.range <- st_crop(env.dat.sf, curr.range)
# Get average of all the environmental variables for the given species
bad.cols <- which(names(env.curr.range) %in% c('cell_id', 'geometry'))
# Mean values for all the calculated climate variables within a species range. 
# Do we want to output these and save somewhere along with species traits
# to avoid having to read in spatial data when fitting models? 
avg.clim.vars <- env.curr.range %>%
  dplyr::select(-cell_id, -geometry) %>%
  st_drop_geometry() %>%
  group_by(var) %>%
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
# Then can just filter trait data accordingly
curr.trait <- gen.dat %>% filter(id == curr.sp)
