################
# Resolve differences in species names and combine bird/mammal datasets:
# - range maps
# - generation time
# - other traits
# - phylogeny
# 
# NOTES:
# - should decide on taxonomic source to match to (itis? col?)
# - for species that don't have db match but common among data sources, just merge?
################


# specify dir -------------------------------------------------------------

#path for data on CY machine - remember trailing slash
# range_map_data_dir <- '~/Downloads/BOTW_2022_2/'
# life_history_dir <- '~/Downloads/'
# paths on HPCC
range_map_data_dir <- '/mnt/research/ibeem/variability/data/L0/ranges/'
life_history_dir <- '/mnt/research/ibeem/variability/data/L0/trait/'
avonet_dir <- '/mnt/research/ibeem/variability/data/L0/trait/'

#directory to save out intermediate file
# out_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/Sample_output/'
out_dir <- '/mnt/research/ibeem/variability/data/L1/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)


# read in data -------------------------------------------------

#Birdlife range maps - takes a few minutes (40 min on HPCC) to read in
#data location: https://drive.google.com/drive/u/1/folders/11eAFmFKUc7tU59IArNEpN5YP_ehGvAwZ
# Get all species
BL_data <- sf::st_read(dsn = paste0(range_map_data_dir, 'BOTW.gdb/a0000000a.gdbtable'))


#bird life history data (processed to L1)
#data location: https://drive.google.com/drive/u/1/folders/18p2Zn3dMA78hdwrCVQBtJ8dnnLEYa7_w
LH_data <- read.csv(paste0(life_history_dir, 'Bird_et_al_gen_length_birds.csv'))

# Avonet trait data
AV_data <- read.csv(paste0(avonet_dir, 'AVONET1_Birdlife.csv'))

#OR read in species names - SEE BELOW


# Process data to get consistent names ------------------------------------

BL_usp <- unique(BL_data$sci_name) # 11187
LH_usp <- unique(LH_data$Sci_name) # 11126
AV_usp <- unique(AV_data$Species1) # 11009
# Save for local manipulation
save(BL_usp, LH_usp, AV_usp, 
     file = paste0(out_dir, 'range/bird-breeding/species-names-all-data.rda'))

# Load in species names in each data set for local manipulation -----------
# load('data/species-names-all-data.rda')
# Convert all names to lower letters
AV.sp.names <- data.frame(name = tolower(AV_usp))
BL.sp.names <- data.frame(name = tolower(BL_usp))
LH.sp.names <- data.frame(name = tolower(LH_usp))


# Using BirdLife taxonomy as the overall master taxonomy ------------------ 

BL.name.df <- data.frame(name = BL.sp.names, 
                         BL.id = 1:dplyr::n_distinct(BL.sp.names))

# Complete DF
full.name.df <- dplyr::full_join(BL.name.df, AV.sp.names, by = 'name') %>%
  dplyr::full_join(LH.sp.names, by = 'name')
write.csv(full.name.df, file = paste0(out_dir, 'trait/bird-names-full.csv'),
          row.names = FALSE)


# Read in full CSV with all species names from all data sources -----------
# The mismatched names between BL and Avonet and Bird et al. were manually
# matched together in Google Sheets.
full.name.dat <- read.csv(paste0(out_dir, 'trait/bird-names-full-matched.csv'))


# Associate each record of the Avonet and Bird et al data with BirdLife ID
# and save to output file

# Bird et al data --------------------- 

LH_data <- LH_data %>%
  dplyr::mutate(Sci_name = tolower(Sci_name)) %>%
  dplyr::left_join(full.name.dat, by = c('Sci_name' = 'name'))

write.csv(LH_data, file = paste0(out_dir, 'trait/bird-et-al-data-with-id.csv'), 
	  row.names = FALSE)

AV_data <- AV_data %>%
  dplyr::mutate(Species1 = tolower(Species1)) %>%
  dplyr::left_join(full.name.dat, by = c('Species1' = 'name'))

write.csv(AV_data, file = paste0(out_dir, 'trait/avonet-with-id.csv'), 
	  row.names = FALSE)

BL_data <- BL_data %>%
  dplyr::mutate(sci_name = tolower(sci_name)) %>%
  dplyr::left_join(BL_data, full.name.dat, by = c('sci_name' = 'name'))


# Bird range data one by one ----------

# Grab breeding/resident ranges, only extant ranges, only native ranges
BL_data_breeding <- BL_data %>%
  dplyr::filter(seasonal %in% 1:2, presence == 1, origin == 1)

# Save out BL species names in case you want them without reading in the whole thing.
BL.unique.ids <- unique(BL_data_breeding$BirdLife_ID)
save(BL.unique.ids, file = paste0(out_dir, 'range/bird-breeding/BL-ids.rda'))

# Takes 40 min or so
for (i in 1:length(BL.unique.ids)) {
  print(paste0("Currently on ", i, " out of ", length(BL.unique.ids)))
  
  tmp <- BL_data_breeding %>%
    dplyr::filter(BirdLife_ID == BL.unique.ids[i])
  
  sf::st_write(tmp, paste0(out_dir, 'range/bird-breeding/', BL.unique.ids[i], 
                           '-breeding.shp'), 
               driver = 'ESRI Shapefile', quiet = TRUE, append = FALSE)
}
