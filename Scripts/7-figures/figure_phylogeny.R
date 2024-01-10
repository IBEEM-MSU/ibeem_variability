#---------------------------------#
# Figure X. Phylogenetic tree ----
#---------------------------------#

# Code based on Monte Neate-Clegg (2023-10-30) - Neate-Clegg et al. 2023 Current Biology Fig. 2

### Load packages ----

library(ape)
library(phytools)
library(tidyverse)

### Data ----

# Load consensus tree: '/mnt/research/ibeem/variability/data/L3/bird-consensus-tree.rda'
load("~/bird-consensus-tree.rda")

# Load bird data
bird_data <- read_csv("~/Documents/Documents/ibeem/main-bird-data-birdtree2.csv") #needs update

head(bird_data)

bird_data %>%
  select(ID, Birdtree_name, Avonet_name, Family, GenLength, Modeled_max_longevity) 

# Extract the tree tip labels (in order), use to reorder dataframe
spp<-bird.phylo$tip.label

# Make spp data frame
spp <- sub("_", " ", spp)

spp_df <- as.data.frame(spp) %>%
  mutate(ID_spp = seq(1:9993)) %>%
  mutate(spp_match = tolower(spp)) %>%
  select(-spp)

# Join to bird_data

bird_data_short <- bird_data %>%
  select(ID, Birdtree_name, Avonet_name, Family, GenLength, Modeled_max_longevity) %>%
  mutate(spp_match = as.character(Birdtree_name)) %>%
  left_join(spp_df, by = "spp_match") %>%
  arrange(ID_spp)

head(bird_data_short)

# Replace tip labels with family names
bird.phylo$tip.label<-as.character(bird_data_short$Family) 

# Summarize data by family

bird_data_short %>%
  group_by(Family) %>%
  summarize(gen_length = mean(GenLength),
            spp_no = n())


