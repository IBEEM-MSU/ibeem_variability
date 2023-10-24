# 6a-mammal-get-con-tree-pieces.R: this script calculates a consensus tree
#                                  individually for 100 trees, then saves each 
#                                  of the resulting 10 consensus trees. The 
#                                  script does this one at a time for each piece,
#                                  as each calculation can take a good chunk of RAM.
rm(list = ls())
library(ape)
library(tidyverse)
library(phytools)

# Get chain number from command line for saving multiple chains ----------
# This code is to extract the current species name from the command line
# to easily run the script for different species
args <- commandArgs(trailingOnly = TRUE)
# Current chain
curr.val <- as.numeric(args[1])
# Alternatively, if not running the script from the command line:
# curr.val <- 1
if(length(args) == 0) base::stop('Need to indicate the current value')
print(curr.val)

# Directory with Phylacine stuff
phy.dir <- '/mnt/research/ibeem/variability/data/L0/trait/phylacine_data/'
out.dir <- '/mnt/research/ibeem/variability/data/L1/phylogeny/'
phy.dat <- read.nexus(paste0(phy.dir, 'Complete_phylogeny.nex'))

# Create consensus tree from 1000 phylogenies
# Calculate consensus trees of 10 chunks (less computing power)
start <- Sys.time()
phy.phylo.01 <- consensus.edges(phy.dat[1:100],
                               method = "least.squares")
write.tree(phy.phylo.01, paste0(out.dir, 'mammal-phylo-1.tre'))
rm(phy.phylo.01)
gc()
print("Currently on tree 2")
phy.phylo.02 <- consensus.edges(phy.dat[101:200],
                               method = "least.squares")
write.tree(phy.phylo.02, paste0(out.dir, 'mammal-phylo-2.tre'))
rm(phy.phylo.02)
gc()
print("Currently on tree 3")
phy.phylo.03 <- consensus.edges(phy.dat[201:300],
                             method = "least.squares")
write.tree(phy.phylo.03, paste0(out.dir, 'mammal-phylo-3.tre'))
rm(phy.phylo.03)
gc()
print("Currently on tree 4")
phy.phylo.04 <- consensus.edges(phy.dat[301:400],
                             method = "least.squares")
write.tree(phy.phylo.04, paste0(out.dir, 'mammal-phylo-4.tre'))
rm(phy.phylo.04)
gc()
print("Currently on tree 5")
phy.phylo.05 <- consensus.edges(phy.dat[401:500],
                             method = "least.squares")
write.tree(phy.phylo.05, paste0(out.dir, 'mammal-phylo-5.tre'))
rm(phy.phylo.05)
gc()
print("Currently on tree 6")
phy.phylo.06 <- consensus.edges(phy.dat[501:600],
                             method = "least.squares")
write.tree(phy.phylo.06, paste0(out.dir, 'mammal-phylo-6.tre'))
rm(phy.phylo.06)
gc()
print("Currently on tree 7")
phy.phylo.07 <- consensus.edges(phy.dat[601:700],
                             method = "least.squares")
write.tree(phy.phylo.07, paste0(out.dir, 'mammal-phylo-7.tre'))
rm(phy.phylo.07)
gc()
print("Currently on tree 8")
phy.phylo.08 <- consensus.edges(phy.dat[701:800],
                             method = "least.squares")
write.tree(phy.phylo.08, paste0(out.dir, 'mammal-phylo-8.tre'))
rm(phy.phylo.08)
gc()
print("Currently on tree 9")
phy.phylo.09 <- consensus.edges(phy.dat[801:900],
                             method = "least.squares")
write.tree(phy.phylo.09, paste0(out.dir, 'mammal-phylo-9.tre'))
rm(phy.phylo.09)
gc()
print("Currently on tree 10")
phy.phylo.10 <- consensus.edges(phy.dat[901:1000],
                             method = "least.squares")
write.tree(phy.phylo.10, paste0(out.dir, 'mammal-phylo-10.tre'))
rm(phy.phylo.10)
gc()
