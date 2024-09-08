# TITLE:            Harmonize bird taxonomy across data sets 
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Range data (BirdLife International), generation length (Bird et al. 2020), life history (AVONET), and phylogeny (BirdTree)
# DATA OUTPUT:      .shp files with bird ranges, and csv of corresponding IDs, bird trait data with IDs, manually matched bird names
# DATE:             May 2023 
# OVERVIEW:         Resolve differences in species names and combine bird dataset: range maps, generation time, other traits, phylogeny
# NOTE:             Requires manually uploaded file: bird-names-full-matched.csv 


# load environment variables ------------------------------------------------

source("./Scripts/0-config.R")


# specify dir -------------------------------------------------------------

#path for data on CY machine - remember trailing slash
# range_map_data_dir <- '~/Downloads/BOTW_2022_2/'
# life_history_dir <- '~/Downloads/'
# paths on HPCC
range_map_data_dir <- paste0(dir, 'data/L0/ranges/')
life_history_dir <- paste0(dir, 'data/L0/trait/')
avonet_dir <- paste0(dir, 'data/L0/trait/')

#directory to save out intermediate file
# out_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/Sample_output/'
range_dir <- paste0(dir, 'data/L1/range/bird-breeding/')
trait_dir <- paste0(dir, 'data/L1/trait/')

ifelse(!dir.exists(file.path(range_dir)), dir.create(file.path(range_dir), recursive=TRUE), FALSE)
ifelse(!dir.exists(file.path(trait_dir)), dir.create(file.path(trait_dir), recursive=TRUE), FALSE)


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)


# read in data -------------------------------------------------

#Birdlife range maps - takes a few minutes (40 min on HPCC) to read in
#data location: https://drive.google.com/drive/u/1/folders/11eAFmFKUc7tU59IArNEpN5YP_ehGvAwZ
# Get all species
BL_data <- sf::st_read(paste0(range_map_data_dir, 'BOTW.gdb/a0000000a.gdbtable'))


#bird life history data (processed to L1)
#data location: https://drive.google.com/drive/u/1/folders/18p2Zn3dMA78hdwrCVQBtJ8dnnLEYa7_w
#'Modeled' values are from Table 4 (some were 'modeled' and others 'hierarchically imputed' - see Bird et al. methods)
LH_data <- read.csv(paste0(life_history_dir, 'Bird_et_al_gen_length_birds.csv'))

# Avonet trait data
AV_data <- read.csv(paste0(avonet_dir, 'AVONET1_Birdlife.csv'))


# Process data to get consistent names ------------------------------------

BL_usp <- unique(BL_data$sci_name) # 11187
LH_usp <- unique(LH_data$Sci_name) # 11126
AV_usp <- unique(AV_data$Species1) # 11009
# Save for local manipulation
save(BL_usp, LH_usp, AV_usp, 
     file = paste0(range_dir, 'species-names-all-data.rda'))


# Load in species names in each data set for local manipulation -----------

# load('data/species-names-all-data.rda')
# Convert all names to lower letters
BL.sp.names <- data.frame(name = tolower(BL_usp))
LH.sp.names <- data.frame(name = tolower(LH_usp))
AV.sp.names <- data.frame(name = tolower(AV_usp))


# Using BirdLife taxonomy as the overall master taxonomy ------------------

BL.name.df <- data.frame(name = BL.sp.names,
                         BL.id = 1:dplyr::n_distinct(BL.sp.names))

# Complete DF
full.name.df <- dplyr::full_join(BL.name.df, AV.sp.names, by = 'name') %>%
  dplyr::full_join(LH.sp.names, by = 'name')
write.csv(full.name.df, file = paste0(trait_dir, 'bird-names-full.csv'),
          row.names = FALSE)


# Read in full CSV with all species names from all data sources -----------

# The mismatched names between BL and Avonet and Bird et al. were manually
# matched together in Google Sheets.
full.name.dat <- read.csv(paste0(trait_dir, 'bird-names-full-matched.csv'))
# Change names for more clarity with linking with different names
colnames(full.name.dat) <- c('birdlifeName', 'birdlifeID', 'notes')
full.name.dat <- full.name.dat %>%
  mutate(birdlifeName = tolower(birdlifeName))

# Read in the Cross-walk that links BirdLife names to BirdTree names. 
crosswalk <- read.csv(paste0(avonet_dir, "avonet-birdtree-crosswalk.csv"))
# Change names for linking with full.name.dat
colnames(crosswalk) <- c('birdlifeName', 'birdtreeName', 'matchType')
# Shoot everything to lowercase to minimize mismatching
crosswalk <- crosswalk %>%
  mutate(birdlifeName = tolower(birdlifeName),
	 birdtreeName = tolower(birdtreeName))

# Join the matched names with the BL names with the BirdTree crosswalk
full.name.dat <- left_join(full.name.dat, crosswalk, by = 'birdlifeName')
# How many BirdTree species are there in the data complete data set
n_distinct(full.name.dat$birdtreeName, na.rm = TRUE)
# 9988 (the full birdtree data set has 9993 unique species, so that's pretty good!)
# Generate a column with unique BirdTree ID
unique.birdtree.names <- sort(unique(full.name.dat$birdtreeName))
# Remove any NA values that that pulled
unique.birdtree.names <- unique.birdtree.names[which(!is.na(unique.birdtree.names))]
# Fill birdtreeIDs
full.name.dat$birdtreeID <- NA
for (i in 1:length(unique.birdtree.names))
{
  #i <- 3360
  tmp <- which(full.name.dat$birdtreeName == unique.birdtree.names[i])
  if (length(tmp) > 0)
  {
    full.name.dat$birdtreeID[tmp] <- i
  }
}


# Bird et al data --------------------- 

LH_data2 <- LH_data %>%
  dplyr::mutate(Sci_name = tolower(Sci_name)) %>%
  dplyr::left_join(full.name.dat, by = c('Sci_name' = 'birdlifeName'))

write.csv(LH_data2, file = paste0(trait_dir, 'bird-et-al-data-with-id.csv'), 
	  row.names = FALSE)

AV_data2 <- AV_data %>%
  dplyr::mutate(Species1 = tolower(Species1)) %>%
  dplyr::left_join(full.name.dat, by = c('Species1' = 'birdlifeName'))

write.csv(AV_data2, file = paste0(trait_dir, 'avonet-with-id.csv'), 
	  row.names = FALSE)

BL_data2 <- BL_data %>%
  dplyr::mutate(sci_name = tolower(sci_name)) %>%
  dplyr::left_join(full.name.dat, by = c('sci_name' = 'birdlifeName'))


# Bird range data one by one ----------

# Grab breeding/resident ranges, only extant ranges, only native ranges
BL_data_breeding <- BL_data2 %>%
  dplyr::filter(seasonal %in% 1:2, 
                presence == 1, 
                origin == 1,
                !is.na(birdtreeID))


# BirdTree ----------------------------

# Save out BirdTree species names in case you want them without reading in the whole thing.
birdTree.unique.ids <- unique(BL_data_breeding$birdtreeID)
save(birdTree.unique.ids, file = paste0(range_dir, 'birdTree-ids.rda'))

# Takes 40 min or so
for (i in 1:length(birdTree.unique.ids))
{
  print(paste0("Currently on ", i, " out of ", length(birdTree.unique.ids)))
  
  tmp <- BL_data_breeding %>%
    dplyr::filter(birdtreeID == birdTree.unique.ids[i])
  
  sf::st_write(tmp, paste0(range_dir, 'birdtree-', birdTree.unique.ids[i], 
                           '-breeding.shp'), 
               driver = 'ESRI Shapefile', quiet = TRUE, append = FALSE)
}
