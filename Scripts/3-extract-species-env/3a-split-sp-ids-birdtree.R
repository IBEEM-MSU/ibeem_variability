# split-sp-ids-birdtree.R: split the saved rda file of BirdTree IDs into multiple 
#                          files for easily running multiple jobs.
rm(list = ls())

# Specify directory -------------------------------------------------------
BT.dir <- '/mnt/research/ibeem/variability/data/L1/range/bird-breeding/'

# Load in full set of current IDs -----------------------------------------
load(paste0(BT.dir, 'birdTree-ids.rda'))
# Remove NA id from BL.unique.ids
BT.unique.ids <- birdTree.unique.ids[!is.na(birdTree.unique.ids)]
n.ids <- length(BT.unique.ids)
n.pieces <- 100
vals <- split(BT.unique.ids, ceiling(seq_along(1:n.ids)/n.pieces))

for (i in 1:length(vals)) {
  ids <- vals[[i]]
  save(ids, file = paste0(BT.dir, 'BTIDsPiece-', i, '.rda'))
}
