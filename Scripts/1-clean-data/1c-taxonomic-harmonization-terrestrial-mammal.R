################
# Resolve differences in species names:
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

#path for data on KK machine - remember trailing slash
range_map_data_dir <-"../terrestrial-mammal/"
life_history_dir <- "../terrestrial-mammal/"

# paths on HPCC
# range_map_data_dir <- '/mnt/research/ibeem/data/L1/range-mammal-clean/'
# life_history_dir <- '/mnt/research/ibeem/data/L0/trait-mammal/'

#directory to save out intermediate file
# out_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/Sample_output/'
# out_dir <- '/mnt/research/ibeem/data/L1/'
out_dir <- '../'


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(taxadb)


# read in data -------------------------------------------------

# Mammal range maps 
#data location: https://drive.google.com/drive/u/1/folders/11eAFmFKUc7tU59IArNEpN5YP_ehGvAwZ
# Get all species
MAM_data <- sf::st_read(dsn = paste0(range_map_data_dir, 'terrestrial-mammals-clean.shp'))


#mammal life history data
#data location: https://drive.google.com/drive/folders/11eAFmFKUc7tU59IArNEpN5YP_ehGvAwZ?usp=drive_link
LH_data <- read.csv(paste0(life_history_dir, 'Generation Length for Mammals.csv'))

#OR read in species names - SEE BELOW


# Process data to get consistent names ------------------------------------

#taxadb tutorial: https://docs.ropensci.org/taxadb/articles/intro.html
#paper: https://besjournals-onlinelibrary-wiley-com.libproxy.lib.unc.edu/doi/full/10.1111/2041-210X.13440

# Using itis for now, but need to consider this more thoroughly
db <- 'itis' #see: https://docs.ropensci.org/taxadb/articles/data-sources.html
# Get associated db ids for bird-life species
MAM_data <- MAM_data %>%
  dplyr::mutate(id = taxadb::get_ids(sci_name, db)) %>%
  dplyr::mutate(accepted_name = taxadb::get_names(id, db)) %>% 
  dplyr:: mutate(itis_id = substr(MAM_data$id, 6,20))
sum(is.na(MAM_data$id))
# Get associated db ids for Pacifici et al. data
LH_data <- LH_data %>%
  dplyr::mutate(id = taxadb::get_ids(Scientific_name, db)) %>%
  dplyr::mutate(accepted_name = taxadb::get_names(id, db)) %>% 
  dplyr:: mutate(itis_id = substr(LH_data$id, 6,20))
sum(is.na(LH_data$id))

MAM_usp <- unique(MAM_data$id) # 5161
LH_usp <- unique(LH_data$id) # 5205

sum(LH_usp %in% MAM_usp)/length(MAM_usp)*100 # 91.6% of mammal species have life history info


miss <- MAM_data[!(MAM_data$id %in% LH_data$id),]

# Temporary processing ----------------------------------------------------
# Just save the data set with the ITIS ids to them, for easy linking across
# data sources. Leaving the missing values as missing values for now just 
# so we have something to work with. 
# Trait data
write.csv(LH_data, file = paste0(out_dir, 'trait-mammal/pacifici-traits-with-id.csv'),
	  row.names = FALSE)

# Mammal range data one by one ----------

# Save out MAM species names in case you want them without reading in the whole thing.
MAM.unique <- MAM_data %>% filter(!is.na(itis_id)) 
MAM.unique.ids <- unique(MAM.unique$itis_id)
save(MAM.unique.ids, file = paste0(out_dir, 'range-mammal/MAM-ids.rda'))
# Takes a 40 min or so
for (i in 1:length(MAM.unique.ids)) {
  print(paste0("Currently on ", i, " out of ", length(MAM.unique.ids)))
  tmp <- MAM_data %>%
    filter(itis_id == MAM.unique.ids[i])
  st_write(tmp, paste0(out_dir, 'range-mammal/', MAM.unique.ids[i], '.shp'), 
           driver = 'ESRI Shapefile', quiet = TRUE, append = FALSE)
}
	
