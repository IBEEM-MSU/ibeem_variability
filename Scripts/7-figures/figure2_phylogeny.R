# TITLE:            Figure 2: Phylogeny  
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Bird consensus tree (4b),model output (5) 
# DATA OUTPUT:      Bird phylogeny figure (labels by order)
# DATE:             January 2024 
# OVERVIEW:         Individual sub-plots were aggregated in Adobe Code based on Monte Neate-Clegg (2023-10-30) - Neate-Clegg et al. 2023 Current Biology Fig. 2. 

### Load packages ----

library(ape)
library(phytools)
library(tidyverse)
library(viridis)
library(plotfunctions)

rm(list = ls())

# load environmental variables ------------------------------------------------

source("./Scripts/0-config.R")

### Data ----

# Load consensus tree: '/mnt/research/ibeem/variability/data/L3/bird-consensus-tree.rda'
load(paste0(dir, "data/L3/bird-consensus-tree.rda"))

# Load bird data

bird_data <- readRDS(paste0(dir, 'Results/bird-gl-phylo-vint-', gl_run_date, 
                          '/bird-gl-phylo-vint-data-', gl_run_date, '.rds'))$pro_data

head(bird_data)

bird_data <- bird_data %>%
  dplyr::mutate(GenLength = exp(lGL)) %>%
  select(ID, species, Order, Family, GenLength, Trophic_niche) 

# Extract spp list from tree and from bird_data
spp_tree <- bird.phylo$tip.label
spp_data <- bird_data$species

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

# Remove extra spp from tree

#bird.phylo$tip.label <- as.character(spp_tree)
bird_tree <- ape::drop.tip(bird.phylo, spp_rm_list)

# Extract tip labels from tree (in order)
spp_match <- bird_tree$tip.label

spp_match_df <- as.data.frame(spp_match) %>%
  mutate(ID_spp = seq(1:7476))

# Add bird data to spp list
bird_data_short <- bird_data %>%
  mutate(spp_match = as.character(species)) %>%
  #select(ID, spp_match, Family, GenLength) %>%
  inner_join(spp_match_df)

# Summarize info by family
bird_data_family <- bird_data_short %>%
  group_by(Family) %>%
  summarize(gen_length = mean(GenLength),
            spp_no = n())

# Change tree tip labels with family names and leave only one tip per family
family_labels <- bird_data_short %>%
  arrange(ID_spp) %>%
  select(Family) %>%
  group_by(Family) %>%
  mutate(fam_no = seq(1:n())) %>%
  mutate(fam_no = as.character(fam_no),
         Family = as.character(Family),
         family_id = paste0(Family, "_", fam_no, "."))

fam_labels_ls <- family_labels %>%
  ungroup() %>%
  select(family_id) %>%
  as.list()

fam_labels_ls <- fam_labels_ls[["family_id"]]

bird_tree$tip.label <- fam_labels_ls

# make list of tips to drop

fam_tips_tokeep <- grep("_1.", fam_labels_ls, fixed = T, value = T)

fam_tips_todrop <- setdiff(fam_labels_ls, fam_tips_tokeep)

# drop tips
tree <- drop.tip(bird_tree, fam_tips_todrop)

# remove id from label
fam_tips <- str_sub(fam_tips_tokeep, 1, -4)

tree$tip.label <- fam_tips
tree

# Extract gen length and species 

bird_data_fam_df <- as.data.frame(fam_tips) %>%
  mutate(Family = as.character(fam_tips)) %>%
  left_join(bird_data_family)

gen_length <- bird_data_fam_df$gen_length
names(gen_length)<-tree$tip.label

spp_no <- bird_data_fam_df$spp_no
names(spp_no)<-tree$tip.label

# Define color ramp with # spp
bird_data_fam_df <- bird_data_fam_df %>%
  mutate(logsp = log1p(spp_no))
logsp <- bird_data_fam_df$logsp
names(logsp) <- tree$tip.label
barcols <- round(logsp*100)-68

#barcols <- round(log1p(bird_data_fam_df$spp_no)*10)-6
# names(barcols) <- tree$tip.label
# ramp<-viridis(max(barcols), option = "A", direction = -1) # Define palette w species richness

# Plot tree with bars

plotTree.wBars(tree,
               x= gen_length,
               type= 'fan',
               tip.labels = F, fsize = .4,
               #col = ramp[barcols],
               col = viridis(max(barcols), option = "A")[barcols],
               border = F,
               width = 3,
               scale = 3)

ramp<-viridis(max(barcols), option = "A")[round(log1p(1:315)*100)-68] # Define palette w species richness
gradientLegend(1:315,
               color = ramp,
               pos = c(-210,70,-190,120),
               side = 2, 
               coords = T,
               length = 0.15,
               depth = 0.025,
               n.seg = c(100, 200),
               inside = T,
               cex = 1)
text(-200,130,'# Species',adj=0.5,pos=3,offset=.5,cex=1.5)

###################
# Tree by order

# Summarize info by order
bird_data_order <- bird_data_short %>%
  group_by(Order) %>%
  summarize(gen_length = mean(GenLength),
            spp_no = n())

# with trophic niche
bird_data_short %>%
  group_by(Order, Trophic_niche) %>%
  summarize(gen_length = mean(GenLength),
            spp_no = n())

# Change tree tip labels with order names and leave only one tip per order
order_labels <- bird_data_short %>%
  arrange(ID_spp) %>%
  select(Order) %>%
  group_by(Order) %>%
  mutate(order_no = seq(1:n())) %>%
  mutate(order_no = as.character(order_no),
         Order = as.character(Order),
         order_id = paste0(Order, "_", order_no, "."))

order_labels_ls <- order_labels %>%
  ungroup() %>%
  select(order_id) %>%
  as.list()

order_labels_ls <- order_labels_ls[["order_id"]]

bird_tree$tip.label <- order_labels_ls

# make list of tips to drop

order_tips_tokeep <- grep("_1.", order_labels_ls, fixed = T, value = T)

order_tips_todrop <- setdiff(order_labels_ls, order_tips_tokeep)

# drop tips
tree <- drop.tip(bird_tree, order_tips_todrop)

# remove id from label
order_tips <- str_sub(order_tips_tokeep, 1, -4)

tree$tip.label <- order_tips
tree

# Extract gen length and species 

bird_data_order_df <- as.data.frame(order_tips) %>%
  mutate(Order = as.character(order_tips)) %>%
  left_join(bird_data_order)

gen_length <- bird_data_order_df$gen_length
names(gen_length)<-tree$tip.label

spp_no <- bird_data_order_df$spp_no
names(spp_no)<-tree$tip.label

# Define color ramp with # spp
bird_data_order_df <- bird_data_order_df %>%
  mutate(logsp = log1p(spp_no))
logsp <- bird_data_order_df$logsp
names(logsp) <- tree$tip.label
barcols <- round(logsp*100)-68

#barcols <- round(log1p(bird_data_fam_df$spp_no)*10)-6
# names(barcols) <- tree$tip.label
# ramp<-viridis(max(barcols), option = "A", direction = -1) # Define palette w species richness

# Plot tree with bars

plotTree.wBars(tree,
               x= gen_length,
               type= 'fan',
               tip.labels = T, fsize = .6,
               #col = ramp[barcols],
               col = viridis(max(barcols), option = "A")[barcols],
               border = F,
               width = 1)

ramp<-viridis(max(barcols), option = "A")[round(log1p(1:4780)*100)-68] # Define palette w species richness
gradientLegend(1:4780,
               color = ramp,
               pos = c(0,2,5,7),
               side = 2, 
               coords = T,
               length = 0.15,
               depth = 0.025,
               n.seg = c(1000, 2000),
               inside = T,
               cex = 1)
text(5,8,'# Species',adj=0.5,pos=3,offset=.5,cex=1.5)



