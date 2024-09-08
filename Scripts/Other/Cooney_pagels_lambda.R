# TITLE:            Check Pagel's Lambda for gen length for Cooney et al. data
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Data from Cooney et al. 2020 Nat Comms
# DATA OUTPUT:      Pagel's Lambda
# DATE:             August 2024



# load packages -----------------------------------------------------------

library(tidyverse)
library(ape)
library(geiger)


# load data ---------------------------------------------------------------

coo <- read.csv('~/Work/Research/Projects/IBEEM_variability/data/L1/trait/Cooney_et_al_2020_Nat_comms.csv') %>%
  dplyr::filter(!is.na(log_generation_length))


# load tree ---------------------------------------------------------------

#load consensus tree - bird.phylo
load(paste0(dir, 'data/L3/bird-consensus-tree.rda'))


# prune tree/process ---------------------------------------------------------

#species not found in both datasets (species to drop from tree)
nm <- setdiff(bird.phylo$tip.label, coo$binomial)

#prune specified tips from tree
pr_tree <- ape::drop.tip(bird.phylo, nm)

#filter coo
coo2 <- dplyr::filter(coo, binomial %ni% nm) %>%
  rename(species = binomial)

#get idx
j_idx3 <- dplyr::left_join(data.frame(species = pr_tree$tip.label),
                           data.frame(idx = 1:NROW(coo2), coo2),
                           by = 'species')

#apply
coo3 <- coo2[j_idx3$idx,]

#make tree binary
pr_tree2 <- ape::multi2di(pr_tree)

#make response into matrix with species as rownames
dd <- dplyr::select(coo3,
                    log_generation_length) %>%
  as.matrix()
row.names(dd) <- coo3$species
NROW(coo3)


# calc Lambda -------------------------------------------------------------

#get estimate of Pagel's lambda
pl <- geiger::fitContinuous(pr_tree2, dd[,'log_generation_length'], 
                            model = "lambda")

