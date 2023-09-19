# split-sp-ids.R: split the saved rda file of BirdLife IDs into multiple 
#                 files for easily running multiple jobs.
rm(list = ls())

# Specify directory -------------------------------------------------------
BL.dir <- '/mnt/research/ibeem/variability/data/L1/range/bird-breeding/'

# Load in full set of current IDs -----------------------------------------
load(paste0(BL.dir, 'BL-ids.rda'))
# Remove NA id from BL.unique.ids
BL.unique.ids <- BL.unique.ids[!is.na(BL.unique.ids)]
n.ids <- length(BL.unique.ids)
n.pieces <- 100
vals <- split(BL.unique.ids, ceiling(seq_along(1:n.ids)/n.pieces))

for (i in 1:length(vals)) {
  ids <- vals[[i]]
  save(ids, file = paste0(BL.dir, 'BLIDsPiece-', i, '.rda'))
}
