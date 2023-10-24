# 6a-bird-get-con-tree-pieces.R: this script calculates a consensus tree
#                                individually for 100 trees, then saves each 
#                                of the resulting 10 consensus trees. The 
#                                script does this one at a time for each piece,
#                                as each calculation can take a good chunk of RAM.
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

# Directory with BirdTree stuff
bt.dir <- '/mnt/research/ibeem/variability/data/L0/phylogeny/'
out.dir <- '/mnt/research/ibeem/variability/data/L1/phylogeny/'
bt.dat <- read.tree(paste0(bt.dir, 'BirdzillaHackett6.tre'))

# Create consensus tree from 1000 phylogenies
# Calculate consensus trees of 10 chunks (less computing power)
start <- Sys.time()
if (curr.val == 1) {
  bt.phylo.01 <- consensus.edges(bt.dat[1:100],
                                 method = "least.squares")
  write.tree(bt.phylo.01, paste0(out.dir, 'bird-phylo-1.tre'))
}
if (curr.val == 2) {
  bt.phylo.02 <- consensus.edges(bt.dat[101:200],
                                 method = "least.squares")
  write.tree(bt.phylo.02, paste0(out.dir, 'bird-phylo-2.tre'))
}
if (curr.val == 3) {
  bt.phylo.03 <- consensus.edges(bt.dat[201:300],
                               method = "least.squares")
  write.tree(bt.phylo.03, paste0(out.dir, 'bird-phylo-3.tre'))
}
if (curr.val == 4) {
  bt.phylo.04 <- consensus.edges(bt.dat[301:400],
                               method = "least.squares")
  write.tree(bt.phylo.04, paste0(out.dir, 'bird-phylo-4.tre'))
}
if (curr.val == 5) {
  bt.phylo.05 <- consensus.edges(bt.dat[401:500],
                               method = "least.squares")
  write.tree(bt.phylo.05, paste0(out.dir, 'bird-phylo-5.tre'))
}
if (curr.val == 6) {
  bt.phylo.06 <- consensus.edges(bt.dat[501:600],
                               method = "least.squares")
  write.tree(bt.phylo.06, paste0(out.dir, 'bird-phylo-6.tre'))
}
if (curr.val == 7) {
  bt.phylo.07 <- consensus.edges(bt.dat[601:700],
                               method = "least.squares")
  write.tree(bt.phylo.07, paste0(out.dir, 'bird-phylo-7.tre'))
}
if (curr.val == 8) {
  bt.phylo.08 <- consensus.edges(bt.dat[701:800],
                               method = "least.squares")
  write.tree(bt.phylo.08, paste0(out.dir, 'bird-phylo-8.tre'))
}
if (curr.val == 9) {
  bt.phylo.09 <- consensus.edges(bt.dat[801:900],
                               method = "least.squares")
  write.tree(bt.phylo.09, paste0(out.dir, 'bird-phylo-9.tre'))
}
if (curr.val == 10) {
  bt.phylo.10 <- consensus.edges(bt.dat[901:1000],
                               method = "least.squares")
  write.tree(bt.phylo.10, paste0(out.dir, 'bird-phylo-10.tre'))
}
