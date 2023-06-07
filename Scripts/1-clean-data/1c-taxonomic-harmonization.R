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
range_map_data_dir <- '/mnt/research/ibeem/data/L0/ranges/'
life_history_dir <- '/mnt/research/ibeem/data/L0/trait/'
elton_dir <- '/mnt/research/ibeem/data/L0/trait/'
avonet_dir <- '/mnt/research/ibeem/data/L0/trait/'

#directory to save out intermediate file
# out_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/Sample_output/'
out_dir <- '/mnt/research/ibeem/data/L1/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(taxadb)


# read in data -------------------------------------------------

#Birdlife range maps - takes a few minutes (40 min on HPCC) to read in
#data location: https://drive.google.com/drive/u/1/folders/11eAFmFKUc7tU59IArNEpN5YP_ehGvAwZ
# Get all species
BL_data <- sf::st_read(dsn = paste0(range_map_data_dir, 'BOTW.gdb/a0000000a.gdbtable'))
# Only select first 100 rows from All_Species, just for testing.
# BL_data <- sf::st_read(dsn = paste0(range_map_data_dir, 'BOTW.gdb/a0000000a.gdbtable'), 
#                        query = "SELECT * FROM ALL_Species limit 100")

# usa <- st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE))
# wkt <- st_as_text(st_geometry(usa))
# # Can read in only ranges that overlap with the US.
# BL_usa <- sf::st_read(dsn = paste0(range_map_data_dir, 'BOTW.gdb/a0000000a.gdbtable'), 
# 		      wkt = wkt)



#bird life history data (processed to L1)
#data location: https://drive.google.com/drive/u/1/folders/18p2Zn3dMA78hdwrCVQBtJ8dnnLEYa7_w
LH_data <- read.csv(paste0(life_history_dir, 'Bird_et_al_gen_length_birds.csv'))

# Elton traits data
ET_data <- read.csv(paste0(elton_dir, 'BirdFuncDat.csv'))

# Avonet trait data
AV_data <- read.csv(paste0(avonet_dir, 'AVONET1_Birdlife.csv'))

#OR read in species names - SEE BELOW


# Process data to get consistent names ------------------------------------

BL_usp <- unique(BL_data$sci_name) # 11187
LH_usp <- unique(LH_data$Sci_name) # 11009
ET_usp <- unique(ET_data$Scientific) # 9994
AV_usp <- unique(AV_data$Species1) # 11126

#taxadb tutorial: https://docs.ropensci.org/taxadb/articles/intro.html
#paper: https://besjournals-onlinelibrary-wiley-com.libproxy.lib.unc.edu/doi/full/10.1111/2041-210X.13440

# For testing: grab the first 40 of each data set
# BL_test <- BL_data %>%
#   arrange(sci_name) %>%
#   slice(1:40)
# LH_test <- LH_data %>%
#   arrange(Sci_name) %>%
#   slice(1:40)
# ET_test <- ET_data %>%
#   arrange(Scientific) %>%
#   slice(1:40)
# AV_test <- AV_data %>%
#   arrange(Species1) %>%
#   slice(1:40)

# Using itis for now, but need to consider this more thoroughly
db <- 'itis' #see: https://docs.ropensci.org/taxadb/articles/data-sources.html
# Get associated db ids for bird-life species
BL_data <- BL_data %>%
  dplyr::mutate(id = taxadb::get_ids(sci_name, db)) %>%
  dplyr::mutate(accepted_name = taxadb::get_names(id, db))
sum(is.na(BL_data$id))
# Get associated db ids for Bird et al. data
LH_data <- LH_data %>%
  dplyr::mutate(id = taxadb::get_ids(Sci_name, db)) %>%
  dplyr::mutate(accepted_name = taxadb::get_names(id, db))
sum(is.na(LH_data$id))
# Elton Traits
ET_data <- ET_data %>%
  dplyr::mutate(id = taxadb::get_ids(Scientific, db)) %>%
  dplyr::mutate(accepted_name = taxadb::get_names(id, db))
sum(is.na(ET_data$id))
# Avonet
AV_data <- AV_data %>%
  dplyr::mutate(id = taxadb::get_ids(Species1, db)) %>%
  dplyr::mutate(accepted_name = taxadb::get_names(id, db))
sum(is.na(AV_data$id))

# Lots of missing values across the data sets

# Temporary processing ----------------------------------------------------
# Just save the data set with the ITIS ids to them, for easy linking across
# data sources. Leaving the missing values as missing values for now just 
# so we have something to work with. 
# Trait data
write.csv(ET_data, file = paste0(out_dir, 'trait/elton-traits-with-id.csv'),
	  row.names = FALSE)
write.csv(AV_data, file = paste0(out_dir, 'trait/avonet-with-id.csv'), 
	  row.names = FALSE)
write.csv(LH_data, file = paste0(out_dir, 'trait/bird-et-al-data-with-id.csv'), 
	  row.names = FALSE)
# Bird range data one by one ----------
# Grab breeding/resident ranges, only extant ranges, only native ranges
BL_data_breeding <- BL_data %>%
  filter(seasonal %in% 1:2, presence == 1, origin == 1)
# Save out BL species names in case you want them without reading in the whole thing.
BL.unique.ids <- unique(BL_data_breeding$id)
save(BL.unique.ids, file = paste0(out_dir, 'range/bird-breeding/BL-ids.rda'))
# Takes a 40 min or so
for (i in 1:length(BL.unique.ids)) {
  print(paste0("Currently on ", i, " out of ", length(BL.unique.ids)))
  tmp <- BL_data_breeding %>%
    filter(id == BL.unique.ids[i])
  st_write(tmp, paste0(out_dir, 'range/bird-breeding/', BL.unique.ids[i], '-breeding.shp'), 
           driver = 'ESRI Shapefile', quiet = TRUE, append = FALSE)
}
	
# Old stuff ---------------------------------------------------------------

#some names match betwen BL and LH but don't have a db id

# BL_NAs <- dplyr::filter(BL_ids, is.na(accepted_name))
# taxadb::filter_name(BL_NAs$sci_name[1])
#IF ITNERMEDIATE OBJECTS EXIST - read in
# BL_usp <- readRDS(paste0(out_dir, 'BL_names.rds'))
# LH_usp <- readRDS(paste0(out_dir, 'LH_names.rds'))
# 
#pretty close to same number of species...positive sign
# BL_usp <- unique(BL_data$sci_name)
# LH_usp <- unique(LH_data$Sci_name)
# length(BL_usp)
# length(LH_usp)
# 
# #clear up memory
# rm(BL_data)
# rm(LH_data)
# gc()
# 
# # #save out intermediate objects - to work with without loading in range maps
# # saveRDS(BL_usp, paste0(out_dir, 'BL_names.rds'))
# # saveRDS(LH_usp, paste0(out_dir, 'LH_names.rds'))
# 
# 
# #differences in the hundreds
# #which species in BL are not in LH
# '%ni%' <- Negate('%in%')
# BL_usp[which(BL_usp %ni% LH_usp)]
# 
# #which species in LH are not in BL
# LH_usp[which(LH_usp %ni% BL_usp)]



# taxize is currently offline, so not going to use that.

# match using taxize package -------------------------------------------------------------------

# BL_test <- sort(BL_usp)[1:40]
# LH_test <- sort(LH_usp)[1:40]
# 
# # https://docs.ropensci.org/taxize/articles/taxize.html
# #get classification and match
# BL_class <- taxize::classification(BL_test, db = 'ncbi')
# LH_class <- taxize::classification(LH_test, db = 'ncbi')
# 
# #find which don't match taxize database
# no_match_BL <- names(which(is.na(BL_class)))
# no_match_LH <- names(which(is.na(LH_class)))
# 
# #resolve names
# #if getting error, website might be down: https://github.com/ropensci/taxize/issues/895
# taxize::gnr_resolve(no_match_BL)
# taxize::gnr_resolve(no_match_LH)
# 
# #matches between BL and LH of first 50
# length(match(BL_class, LH_class))
