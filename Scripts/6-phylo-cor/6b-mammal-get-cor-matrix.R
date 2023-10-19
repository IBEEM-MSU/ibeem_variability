# 6b-mammal-get-cor-matrix.R: script to calculate the final consensus tree
#                             and extract a correlation matrix from it for
#                             the mammal data
rm(list = ls())
library(ape)
library(tidyverse)
library(phytools)

mammal.dir <- '/mnt/research/ibeem/variability/data/L1/phylogeny/'
out.dir <- '/mnt/research/ibeem/variability/data/L3/'

# Read in each of the 10 consensus tree chunks
mammal.phylo.1 <- read.tree(paste0(mammal.dir, "mammal-phylo-1.tre"))
mammal.phylo.2 <- read.tree(paste0(mammal.dir, "mammal-phylo-2.tre"))
mammal.phylo.3 <- read.tree(paste0(mammal.dir, "mammal-phylo-3.tre"))
mammal.phylo.4 <- read.tree(paste0(mammal.dir, "mammal-phylo-4.tre"))
mammal.phylo.5 <- read.tree(paste0(mammal.dir, "mammal-phylo-5.tre"))
mammal.phylo.6 <- read.tree(paste0(mammal.dir, "mammal-phylo-6.tre"))
mammal.phylo.7 <- read.tree(paste0(mammal.dir, "mammal-phylo-7.tre"))
mammal.phylo.8 <- read.tree(paste0(mammal.dir, "mammal-phylo-8.tre"))
mammal.phylo.9 <- read.tree(paste0(mammal.dir, "mammal-phylo-9.tre"))
mammal.phylo.10 <- read.tree(paste0(mammal.dir, "mammal-phylo-10.tre"))

# Then create consensus tree of the 10 consensus trees
mammal.phylo <- consensus.edges(
  as.multiPhylo(c(mammal.phylo.1, mammal.phylo.2,
                  mammal.phylo.3, mammal.phylo.4,
                  mammal.phylo.5, mammal.phylo.6,
                  mammal.phylo.7, mammal.phylo.8,
                  mammal.phylo.9, mammal.phylo.10)),
  method = "least.squares")
 
# Check that consensus tree is ultrametric
is.ultrametric(mammal.phylo)

# write out tree
save(mammal.phylo, file = paste0(out.dir, 'mammal-consensus-tree.rda'))

# Calculate covariance matrix
cov.mat <- ape::vcv.phylo(mammal.phylo)

# Convert covariance matrix to correlation matrix
cor.mat <- cov2cor(cov.mat)

# Get corr matrix for species in the final main file, in the correct order 
mammal.dat <- read.csv(paste0(out.dir, "main-mammal-data.csv"))
dimnames(cor.mat)[[1]] <- tolower(str_replace_all(dimnames(cor.mat)[[1]], '_', ' '))
dimnames(cor.mat)[[2]] <- tolower(str_replace_all(dimnames(cor.mat)[[2]], '_', ' '))
final.cor.mat <- matrix(NA, nrow(mammal.dat), nrow(mammal.dat))
rownames(final.cor.mat) <- mammal.dat$Accepted_name
colnames(final.cor.mat) <- mammal.dat$Accepted_name
ordered.indx <- rep(NA, nrow(final.cor.mat))
for (i in 1:nrow(final.cor.mat)) {
  ordered.indx[i] <- which(dimnames(cor.mat)[[1]] == tolower(mammal.dat$Accepted_name[i]))
}
for (i in 1:nrow(final.cor.mat)) {
  final.cor.mat[i, ] <- cor.mat[ordered.indx[i], ordered.indx]
}

save(final.cor.mat, file = paste0(out.dir, 'mammal-phylo-cor-matrix.rda'))

