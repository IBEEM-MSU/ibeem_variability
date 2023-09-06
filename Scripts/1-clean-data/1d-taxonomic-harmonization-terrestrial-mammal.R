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
range_map_data_dir <-"./data/L0/ranges/"
life_history_dir <- "./data/L0/trait/"

# paths on HPCC
# range_map_data_dir <- '/mnt/research/ibeem/data/L1/range-mammal-clean/'
# life_history_dir <- '/mnt/research/ibeem/data/L0/trait-mammal/'

#directory to save out intermediate file
# out_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/Sample_output/'
# out_dir <- '/mnt/research/ibeem/data/L1/'
out_dir <- './data/L1/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(rredlist)
library(progress)
#library(taxadb)


# read in data -------------------------------------------------

# Mammal range maps 
#data location: https://drive.google.com/drive/u/1/folders/11eAFmFKUc7tU59IArNEpN5YP_ehGvAwZ
# Get all species
MAM_data <- sf::st_read(dsn = paste0(range_map_data_dir, 'terrestrial-mammals-clean.shp'))


#mammal life history data
#data location: https://drive.google.com/drive/folders/11eAFmFKUc7tU59IArNEpN5YP_ehGvAwZ?usp=drive_link
LH_data <- read.csv(paste0(life_history_dir, 'Generation Length for Mammals.csv'))

NM_data <- read.csv(paste0(life_history_dir, 'phylacine_data/Synonymy_table_valid_species_only.csv'))

PH_data <- read.csv(paste0(life_history_dir, 'phylacine_data/Trait_data.csv'))

NM <- NM_data %>% 
  mutate(name_phylacine = paste(Genus.1.2, Species.1.2),
         name_phylacine_1.1 = paste(Genus.1.1, Species.1.1),
         name_phylacine_1.0 = paste(Genus.1.0, Species.1.0),
         name_EltonTraits = paste(EltonTraits.1.0.Genus, EltonTraits.1.0.Species),
         name_iucn_2016 = paste(IUCN.2016.3.Genus, IUCN.2016.3.Species)) %>% 
  mutate (name_phylacine = gsub("000 Species not accepted 000 Species not accepted", NA, name_phylacine),
          name_phylacine_1.1 = gsub("000 Species not accepted 000 Species not accepted", NA, name_phylacine_1.1),
          name_phylacine_1.0 = gsub("000 Species not accepted 000 Species not accepted", NA, name_phylacine_1.0),
          name_EltonTraits = gsub("000 Species not accepted 000 Species not accepted", NA, name_EltonTraits),
          name_iucn_2016 = gsub("000 Species not accepted 000 Species not accepted", NA, name_iucn_2016)) %>% 
  select(name_phylacine,
         name_phylacine_1.1,
         name_phylacine_1.0,
         name_EltonTraits,
         name_iucn_2016)

LH <- LH_data %>% 
  select(Scientific_name) %>% 
  rename(name_pacifici = Scientific_name)

MAM <- MAM_data %>% 
  st_drop_geometry() %>% 
  select(sci_name) %>% 
  mutate(id = 1:n_distinct(sci_name)) %>% 
  rename(name_iucn = sci_name) 

# Join data, use Phylacine synonyms to match names
t <- MAM %>% mutate("name_phylacine" = NA,
                    "name_pacifici" = NA)
for(i in 1:nrow(t)){
  j <- ifelse(t$name_iucn[i] %in% NM$name_phylacine,
              which(t$name_iucn[i] == NM$name_phylacine),
              ifelse(t$name_iucn[i] %in% NM$name_phylacine_1.1,
                     which(t$name_iucn[i] == NM$name_phylacine_1.1),
                     ifelse(t$name_iucn[i] %in% NM$name_phylacine_1.0,
                            which(t$name_iucn[i] == NM$name_phylacine_1.0),
                            ifelse(t$name_iucn[i] %in% NM$name_EltonTraits,
                                   which(t$name_iucn[i] == NM$name_EltonTraits),
                                   NA))))
  t$name_phylacine[i] <- ifelse(is.na(j), NA,
                                NM$name_phylacine[j])
}
t <- left_join(t, NM, by = "name_phylacine")
for(i in 1:nrow(t)){
  t$name_pacifici[i] <- ifelse(t$name_iucn[i] %in% LH$name_pacifici,
                               t$name_iucn[i],
                               ifelse(t$name_phylacine[i] %in% LH$name_pacifici,
                                      t$name_phylacine[i],
                                      ifelse(t$name_phylacine_1.1[i] %in% LH$name_pacifici,
                                             t$name_phylacine_1.1[i],
                                             ifelse(t$name_phylacine_1.0[i] %in% LH$name_pacifici,
                                                    t$name_phylacine_1.0[i],
                                                    ifelse(t$name_EltonTraits[i] %in% LH$name_pacifici,
                                                           t$name_EltonTraits[i], NA)))))
}
t <- t %>% select(name_iucn, name_phylacine, name_pacifici)

# Use IUCN Red List synonyms to match names
# Need IUCN Red List key to access API
t_complete <- t[rowSums(is.na(t)) == 0,]
t_no_match <- t[rowSums(is.na(t)) > 0,]
sum(!is.na(t_no_match$name_iucn))

pb <- progress_bar$new(
  format = "[:bar] :percent in :elapsed with :eta remaining",
  total = nrow(t_no_match), clear = FALSE, width= 60)
for(i in 1:nrow(t_no_match)){
  syn <- rredlist::rl_synonyms(name = t_no_match$name_iucn[i],
                               key = IUCN_REDLIST_KEY,
                               parse = TRUE)$result
  new_name_phylacine <- ifelse(length(syn)==0,
                               t_no_match$name_phylacine[i],
                               syn$synonym[which(syn$synonym %in% NM$name_phylacine)])
  new_name_pacifici <- ifelse(length(syn)==0,
                              t_no_match$name_pacifici[i],
                              syn$synonym[which(syn$synonym %in% LH$name_pacifici)])
  t_no_match$name_phylacine[i] <- ifelse(is.na(t_no_match$name_phylacine[i]),
                                         new_name_phylacine,
                                         t_no_match$name_phylacine[i])
  t_no_match$name_pacifici[i] <- ifelse(is.na(t_no_match$name_pacifici[i]),
                                        new_name_pacifici,
                                        t_no_match$name_pacifici[i])
  pb$tick()
}

# Add newly matched names to list of complete names, keep unmatched names
t_new_names <- t_no_match[rowSums(is.na(t_no_match)) == 0,]
t_complete <- bind_rows(t_complete, t_new_names)
t_no_match <- t_no_match[rowSums(is.na(t_no_match)) > 0,]
sum(!is.na(t_no_match$name_iucn))

# Save files
write.csv(t_complete, file = paste0(out_dir, 'mammal-names-matched.csv'), row.names = FALSE)
write.csv(t_no_match, file = paste0(out_dir, 'mammal-names-no-match.csv'), row.names = FALSE)
write.csv(NM, file = paste0(out_dir, 'mammal-names-phylacine.csv'), row.names = FALSE)
write.csv(LH, file = paste0(out_dir, 'mammal-names-pacifici.csv'), row.names = FALSE)

#
#
#
#
#
#





# test <- read.csv('./data/L1/trait/mammal-names-full.csv')


##################################
##################################
##################################



##################################
# for loop on iucn species and then get list of synonyms....
# see if those synonyms match 


###################################################################################

# Process data to get consistent names ------------------------------------

#taxadb tutorial: https://docs.ropensci.org/taxadb/articles/intro.html
#paper: https://besjournals-onlinelibrary-wiley-com.libproxy.lib.unc.edu/doi/full/10.1111/2041-210X.13440
# 
# # Using itis for now, but need to consider this more thoroughly
# db <- 'itis' #see: https://docs.ropensci.org/taxadb/articles/data-sources.html
# # Get associated db ids for bird-life species
# MAM <- MAM_data %>%
#   dplyr::mutate(id = taxadb::get_ids(sci_name, db)) %>%
#   dplyr::mutate(accepted_name = taxadb::get_names(id, db)) %>% 
#   dplyr:: mutate(itis_id = substr(id, 6,20)) 
# sum(is.na(MAM_data$id))
# # Get associated db ids for Pacifici et al. data
# LH <- LH_data %>%
#   dplyr::mutate(id = taxadb::get_ids(Scientific_name, db)) %>%
#   dplyr::mutate(accepted_name = taxadb::get_names(id, db)) %>% 
#   dplyr:: mutate(itis_id = substr(id, 6,20))
# 
# NM <- NM_data %>%
#   dplyr::mutate(id = taxadb::get_ids(Binomial.1.2, db)) %>%
#   dplyr::mutate(accepted_name = taxadb::get_names(id, db)) %>% 
#   dplyr:: mutate(itis_id = substr(id, 6,20))
# 
# MAM_usp <- unique(MAM_data$id) # 5161
# LH_usp <- unique(LH_data$id) # 5205
# NM_usp <- unique(NM_data$id) # 5304
# 
# sum(LH_usp %in% MAM_usp)/length(MAM_usp)*100 # 91.6% of mammal species have life history info
# sum(NM_usp %in% MAM_usp)/length(MAM_usp)*100 # 94.6% of mammal species have phylogenetic info
# 

# Temporary processing ----------------------------------------------------
# Just save the data set with the ITIS ids to them, for easy linking across
# data sources. Leaving the missing values as missing values for now just 
# so we have something to work with. 
# Trait data
# write.csv(LH_data, file = paste0(out_dir, 'trait-mammal/pacifici-traits-with-id.csv'),
# 	  row.names = FALSE)
# 
# # Mammal range data one by one ----------
# 
# # Save out MAM species names in case you want them without reading in the whole thing.
# MAM.unique <- MAM_data %>% filter(!is.na(itis_id)) 
# MAM.unique.ids <- unique(MAM.unique$itis_id)
# save(MAM.unique.ids, file = paste0(out_dir, 'range-mammal/MAM-ids.rda'))
# # Takes a 40 min or so
# for (i in 1:length(MAM.unique.ids)) {
#   print(paste0("Currently on ", i, " out of ", length(MAM.unique.ids)))
#   tmp <- MAM_data %>%
#     filter(itis_id == MAM.unique.ids[i])
#   st_write(tmp, paste0(out_dir, 'range-mammal/', MAM.unique.ids[i], '.shp'), 
#            driver = 'ESRI Shapefile', quiet = TRUE, append = FALSE)
# }
# 	
