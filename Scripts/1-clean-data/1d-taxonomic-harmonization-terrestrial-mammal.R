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
MAM_data <- sf::st_read(dsn = paste0(range_map_data_dir, 'terrestrial-mammals-clean.shp')) %>% 
  rename(name_iucn = sci_name) 


#mammal life history data
#data location: https://drive.google.com/drive/folders/11eAFmFKUc7tU59IArNEpN5YP_ehGvAwZ?usp=drive_link
LH_data <- read.csv(paste0(life_history_dir, 'Generation Length for Mammals.csv')) %>% 
  rename(name_pacifici = Scientific_name)

NM_data <- read.csv(paste0(life_history_dir, 'phylacine_data/Synonymy_table_valid_species_only.csv'))

PH_data <- read.csv(paste0(life_history_dir, 'phylacine_data/Trait_data.csv')) %>% 
  mutate(name_phylacine = gsub("_", " ", Binomial.1.2)) %>% select(-Binomial.1.2)


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
  select(name_pacifici) 

MAM <- MAM_data %>% 
  st_drop_geometry() %>% 
  select(name_iucn) %>% 
  mutate(id = 1:n_distinct(name_iucn)) 

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


# Create master names list and join to spatial data ------------------------------------

# Load in taxonomy lists 
full_names <- read.csv('./data/L1/trait/mammal-names-matched.csv') # Names with no problems (5146)
matched_names <- read.csv('./data/L1/trait/mammal-names-fixed.csv') %>% # Names Peter manually matched 
  select(name_iucn, name_pacifici, name_phylacine)

# Join and remove incomplete data (i.e., spp not present in >= 1 of the DBs)
master_names <- rbind(full_names, matched_names) 
master_names[master_names == ""] <- NA
master_names <- master_names[complete.cases(master_names),]

# Eval number of unique iucn species lumped together as one in pacifici/phylacine
lumps <- master_names %>% group_by(name_pacifici) %>% summarize(n = n()) %>% filter(n > 1) # 209
lumps2 <- master_names %>% group_by(name_phylacine) %>% summarize(n = n()) %>% filter(n > 1) # 173

# Join range polys to names list 
# (right join keeps spatial attributes while removing polys that we have range info, but no traits)
master_names_sf <- right_join(MAM_data, master_names, by = "name_iucn")

# # Mammal range data one by one ----------

# Save out MAM species names in case you want them without reading in the whole thing.
MAM.unique.ids <- unique(master_df$id_no)

save(MAM.unique.ids, file = paste0(out_dir, 'range-mammal/MAM-ids.rda'))

# Takes a 40 min or so
for (i in 1:length(MAM.unique.ids)) {
  print(paste0("Currently on ", i, " out of ", length(MAM.unique.ids)))
  tmp <- MAM_data %>%
    filter(id_no == MAM.unique.ids[i])
  st_write(tmp, paste0(out_dir, 'range-mammal/', MAM.unique.ids[i], '.shp'),
           driver = 'ESRI Shapefile', quiet = TRUE, append = FALSE)
}

