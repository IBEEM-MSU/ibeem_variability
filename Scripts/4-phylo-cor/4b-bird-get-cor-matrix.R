# TITLE:            Generate correlation matrix for bird phylogeny 
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Peter Williams, Jeff Dozer, Adriana Uscanga, Lala Kounta, Kelly Kapsar, Phoebe Zarnetske, Pat Bills
# DATA INPUT:       Partial phylogenetic trees for bird taxonomy from 4a
# DATA OUTPUT:      Correlation matrix (rda) for full bird phylogeny
# DATE:             October 2023 
# OVERVIEW:         Script to calculate the final consensus tree and extract a correlation matrix from it. 

rm(list = ls())

library(ape)
library(tidyverse)
library(phytools)

bt.dir <- '/mnt/research/ibeem/variability/data/L1/phylogeny/'
out.dir <- '/mnt/research/ibeem/variability/data/L3/'

# Read in each of the 10 consensus tree chunks
bt.phylo.1 <- read.tree(paste0(bt.dir, "bird-phylo-1.tre"))
bt.phylo.2 <- read.tree(paste0(bt.dir, "bird-phylo-2.tre"))
bt.phylo.3 <- read.tree(paste0(bt.dir, "bird-phylo-3.tre"))
bt.phylo.4 <- read.tree(paste0(bt.dir, "bird-phylo-4.tre"))
bt.phylo.5 <- read.tree(paste0(bt.dir, "bird-phylo-5.tre"))
bt.phylo.6 <- read.tree(paste0(bt.dir, "bird-phylo-6.tre"))
bt.phylo.7 <- read.tree(paste0(bt.dir, "bird-phylo-7.tre"))
bt.phylo.8 <- read.tree(paste0(bt.dir, "bird-phylo-8.tre"))
bt.phylo.9 <- read.tree(paste0(bt.dir, "bird-phylo-9.tre"))
bt.phylo.10 <- read.tree(paste0(bt.dir, "bird-phylo-10.tre"))

# Then create consensus tree of the 10 consensus trees
bird.phylo <- consensus.edges(
  as.multiPhylo(c(bt.phylo.1, bt.phylo.2,
                  bt.phylo.3, bt.phylo.4,
                  bt.phylo.5, bt.phylo.6,
                  bt.phylo.7, bt.phylo.8,
                  bt.phylo.9, bt.phylo.10)),
  method = "least.squares")
 
# Check that consensus tree is ultrametric
is.ultrametric(bird.phylo)

# write out tree
save(bird.phylo, file = paste0(out.dir, 'bird-consensus-tree.rda'))

# Calculate covariance matrix
cov.mat <- ape::vcv.phylo(bird.phylo)

# Convert covariance matrix to correlation matrix
cor.mat <- cov2cor(cov.mat)

# Get corr matrix for species in the final main file, in the correct order 
bird.dat <- read.csv(paste0(out.dir, "main-bird-data-birdtree2.csv"))
dimnames(cor.mat)[[1]] <- tolower(str_replace_all(dimnames(cor.mat)[[1]], '_', ' '))
dimnames(cor.mat)[[2]] <- tolower(str_replace_all(dimnames(cor.mat)[[2]], '_', ' '))
final.cor.mat <- matrix(NA, nrow(bird.dat), nrow(bird.dat))
rownames(final.cor.mat) <- bird.dat$Birdtree_name
colnames(final.cor.mat) <- bird.dat$Birdtree_name
ordered.indx <- rep(NA, nrow(final.cor.mat))
for (i in 1:nrow(final.cor.mat)) {
  ordered.indx[i] <- which(dimnames(cor.mat)[[1]] == bird.dat$Birdtree_name[i])
}
for (i in 1:nrow(final.cor.mat)) {
  final.cor.mat[i, ] <- cor.mat[ordered.indx[i], ordered.indx]
}

save(final.cor.mat, file = paste0(out.dir, 'bird-phylo-cor-matrix.rda'))

