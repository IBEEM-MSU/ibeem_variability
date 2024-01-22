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
spp <- bird.phylo$tip.label

# Change to lower case and remove underscore so that names match "Birdtree_name" in bird dataset
spp <- sub("_", " ", spp)
spp <- tolower(spp)

# Make spp data frame with spp names and ordered ID
spp_df <- as.data.frame(spp) %>%
  mutate(ID_spp = seq(1:9993)) %>%
  mutate(spp_match = spp) %>%
  select(-spp)

# Remove spp that are not in bird dataset from tree
bird_data <- bird_data %>%
  mutate(spp_match = as.character(Birdtree_name))

spp_rm <- spp_df %>%
  anti_join(bird_data, by = "spp_match")

spp_rm_list <- spp_rm %>%
  select(spp_match) %>%
  as.list()
spp_rm_list <- spp_rm_list[["spp_match"]]

bird_tree$tip.label <- as.character(spp)

bird_tree <- ape::drop.tip(bird_tree, spp_rm_list)

spp_tree <- bird_tree$tip.label

spp_tree_df <- as.data.frame(spp_tree) %>%
  mutate(ID_spp = seq(1:9648)) %>%
  mutate(spp_match = spp_tree) %>%
  select(-spp_tree)

# Join to bird_data

bird_data_short <- inner_join(spp_tree_df, bird_data, by = "spp_match")
head(bird_data_short)

  

# Replace tip labels with family names
bird.phylo$tip.label<-as.character(bird_data_short$Family) 

# Summarize data by family

bird_data_short %>%
  group_by(Family) %>%
  summarize(gen_length = mean(GenLength),
            spp_no = n())


