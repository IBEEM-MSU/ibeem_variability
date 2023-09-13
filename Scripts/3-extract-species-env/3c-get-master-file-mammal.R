# 3c-get-master-file.R: script to generate a single master file where each 
#                       row is a species and the columns contain the traits, 
#                       and climate variables

rm(list = ls())

# load packages ------------------------------------------------------------

library(tidyverse)


# Specify directories -----------------------------------------------------

dir <- '/mnt/research/ibeem/variability/'
# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# dir <- "./"

#path for data on KK machine - remember trailing slash
# life_history_dir <- "./data/L0/trait/"


# Get the climate data ----------------------------------------------------

# Loads the full set of current BL IDs 
load(paste0(dir, 'data/L1/range-mammal/MAM-ids.rda'))
# Load the first chunk of species and climate averages to get number of variables
# Loads object called "env.out"
load(paste0(dir, 'data/L2/range-env-pieces-mammal/summarized-data-piece-1'))

#build into list
climate.list <- list()
clim.files <- list.files(paste0(dir, 'data/L2/range-env-pieces-mammal/'), 
                         full.names = TRUE)
for (i in 1:length(clim.files)) {
  print(i)
  load(clim.files[i])
  climate.list[[i]] <- env.out
}
# Put it in a data frame
climate.df <- do.call("rbind", climate.list)
# How many species have NA for temp_sd_year? Just 2
sum(apply(dplyr::select(climate.df, temp_sd_year), 
          1, function(a) sum(is.na(a))) > 0)

# Add names to climate data 
master_names <- read.csv(paste0(dir, "data/L1/trait-mammal/mammal-names-master.csv"))
climate.df <- dplyr::left_join(climate.df, master_names, by = c("ID" = "id_no"))


# Load in the trait data and join ------------------------------------------

# Life history from Pacifici et al. 
LH_data <- read.csv(paste0(dir, 'data/L0/trait/Generation Length for Mammals.csv')) %>% 
  dplyr::rename(name_pacifici = Scientific_name)

# Trait data from Phylacine
PH_data <- read.csv(paste0(dir, 'data/L0/trait/phylacine_data/Trait_data.csv')) %>% 
  dplyr::mutate(name_phylacine = gsub("_", " ", Binomial.1.2)) %>% 
  dplyr::select(-Binomial.1.2)

# Specify which source the data came from 
# (LH = life history from Pacifici)
# (PH = Phylacine)
colnames(LH_data) <- paste0("LH_", colnames(LH_data))
colnames(PH_data) <- paste0("PH_", colnames(PH_data))

## Options for dealing with "splitters" (i.e., multiple IUCN spp, but only one in trait DBs)
# 1. Repeat the trait data across the split-up species. Chances are they're not that different in size/generation length anyway.

# Remove random duplicate rows in LH_data
LH_data <- LH_data %>% filter(LH_TaxID != 198920) # Duplicate row for Neophocaena phocaenoides
LH_data <- LH_data %>% filter(LH_TaxID != 8212) # Duplicate row for Neophocaena phocaenoides

climate.df <- dplyr::left_join(climate.df, LH_data, 
                        by = c("name_pacifici" = "LH_name_pacifici")) %>% 
  dplyr::select(-LH_TaxID, -LH_Order, -LH_Family, -LH_Genus)

climate.df <- dplyr::left_join(climate.df, PH_data, 
                        by = c("name_phylacine" = "PH_name_phylacine")) %>% 
  dplyr::select(-PH_Order.1.2, -PH_Family.1.2, -PH_Genus.1.2, -PH_Species.1.2)

# 2. Do a spatial join of all the ranges to effectively re-lump the split up species together again. 
# 3. Remove any species that's been lumped/split. (least favorite option)

# Clean up data 
# only species with values for mean temp
main.dat <- climate.df %>% 
  dplyr::select(-name_phylacine, -name_pacifici, -X) %>% 
  dplyr::rename(Accepted_name = name_iucn) %>% 
  dplyr::relocate(ID, Accepted_name) %>% 
  dplyr::mutate(precip_cv_space = precip_sd_space / precip_mean,
                dhi_cum_cv_space = dhi_cum_sd_space / dhi_cum_mean) %>%
  dplyr::filter(!is.na(temp_mean)) %>%
  dplyr::relocate(precip_cv_space, .after = precip_sd_space) %>%
  dplyr::relocate(dhi_cum_cv_space, .after = dhi_cum_sd_space) #%>%
# dplyr::filter(!is.na(range_size_km2))


# Save output -----------------------------------------------------

write.csv(main.dat, file = paste0(dir, 'data/L3/main-mammal-data.csv'), 
          row.names = FALSE)


