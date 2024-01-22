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

# Extract spp list from tree and from bird_data
spp_tree <- bird.phylo$tip.label
spp_data <- bird_data$Birdtree_name

# Change to lower case and remove underscore so that names match "Birdtree_name" in bird dataset
spp_tree <- sub("_", " ", spp_tree)
spp_tree <- tolower(spp_tree)

# Make spp data frames with spp names and ordered ID
spp_tree_df <- as.data.frame(spp_tree) %>%
  mutate(ID_spp = seq(1:9993)) %>%
  mutate(spp_match = spp_tree) %>%
  select(-spp_tree)

spp_data_df <- as.data.frame(spp_data) %>%
  mutate(spp_match = spp_data) %>%
  select(-spp_data)

# Inner join
spp_list <- inner_join(spp_tree_df, spp_data_df)

# Anti join to get a list of spp that need to be removed from tree
spp_rm <- anti_join(spp_tree_df, spp_list)
spp_rm_list <- as.list(spp_rm)
spp_rm_list <- spp_rm_list[["spp_match"]]

# Remove spp with missing info from both tree and bird_data

bird.phylo$tip.label <- as.character(spp_tree)
bird_tree <- ape::drop.tip(bird.phylo, spp_rm_list)

# Extract tip labels from tree (in order)
spp_match <- bird_tree$tip.label

spp_match_df <- as.data.frame(spp_match) %>%
  mutate(ID_spp = seq(1:9648))

# Add bird data to spp list
bird_data_short <- bird_data %>%
  mutate(spp_match = as.character(Birdtree_name)) %>%
  select(ID, spp_match, Family, GenLength) %>%
  inner_join(spp_match_df)

bird_data_short <- bird_data_short %>%
  filter(!duplicated(ID_spp))

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


