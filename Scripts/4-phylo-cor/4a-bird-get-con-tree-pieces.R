# TITLE:            Generate consensus tree for bird phylogeny 
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Phylogenetic data from BirdTree
# DATA OUTPUT:      Partial phylogenetic trees for bird taxonomy
# DATE:             October 2023 
# OVERVIEW:         This script calculates a consensus tree individually for 100 trees, then saves each of the resulting 10 consensus trees. The script does this one at a time for each piece,as each calculation can take a good chunk of RAM.
# NOTES:            Run in R version 4.2.1 


rm(list = ls())


# load environment variables ------------------------------------------------

source('./Scripts/0-config.R')


# load packages ------------------------------------------------

library(ape)
library(tidyverse)
library(phytools)

# Directory with BirdTree stuff
bt.dir <- paste0(dir, '/data/L0/phylogeny/')
out.dir <- paste0(dir, '/data/L1/phylogeny/')
bt.dat <- read.tree(paste0(bt.dir, 'BirdzillaHackett6.tre'))


# data --------------------------------------------------------------------

#tree data were compiled from `HackettStage2_XXXX.zip` downloads from Birdtree.org
#source data: 
#Jetz, W., G. H. Thomas, J. B. Joy, K. Hartmann, and A. O. Mooers. 2012. The global diversity of birds in space and time. Nature 491:444–448.(DOI: [https://doi.org/10.1038/nature11631](https://doi.org/10.1038/nature11631))


# Running from desktop (or interactive session) ----------

# Create consensus tree from 1000 phylogenies
# Calculate consensus trees of 10 chunks (less computing power)
chunks <- split(1:1000, ceiling(seq_along(1:1000)/100))

for(i in 1:length(chunks)){
  print(paste0("Processing chunk ", i))
  start <- Sys.time()
  curr.val <- i
  chunk <- chunks[[i]]
  bt.phylo <- consensus.edges(bt.dat[chunk],
                              method = "least.squares")
  write.tree(bt.phylo, paste0(out.dir, paste0('bird-phylo-',i,'.tre')))
  rm(bt.phylo)
  print(start-Sys.time())
}


# Running from command line ----------

# args <- commandArgs(trailingOnly = TRUE)
# # Current chain
# curr.val <- as.numeric(args[1])
# if(length(args) == 0) base::stop('Need to indicate the current value')
# print(curr.val)
# 
# # Create consensus tree from 1000 phylogenies
# # Calculate consensus trees of 10 chunks (less computing power)
# chunks <- split(1:1000, ceiling(seq_along(1:1000)/100))
# 
# start <- Sys.time()
# curr.val <- i
# chunk <- chunks[[i]]
# bt.phylo <- consensus.edges(bt.dat[chunk],
#                                method = "least.squares")
# write.tree(bt.phylo, paste0(out.dir, paste0('bird-phylo-',i,'.tre')))
# rm(bt.phylo)
# print(start-Sys.time())
