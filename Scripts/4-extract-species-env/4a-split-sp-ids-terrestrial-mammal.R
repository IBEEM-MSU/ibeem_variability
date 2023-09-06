# split-sp-ids.R: split the saved rda file of BirdLife IDs into multiple 
#                 files for easily running multiple jobs.
rm(list = ls())

# Specify directory -------------------------------------------------------
# MAM.dir <- "./data/L1/range-mammal/"
MAM.dir <- "/mnt/research/ibeem/data/L1/range-mammal/"

# Load in full set of current IDs -----------------------------------------
load(paste0(MAM.dir, 'MAM-ids.rda'))
# Remove NA id from MAM.unique.ids
MAM.unique.ids <- MAM.unique.ids[!is.na(MAM.unique.ids)]
n.ids <- length(MAM.unique.ids)
n.pieces <- 100
vals <- split(MAM.unique.ids, ceiling(seq_along(1:n.ids)/n.pieces))

for (i in 1:length(vals)) {
  ids <- vals[[i]]
  save(ids, file = paste0(MAM.dir, 'MAMIDsPiece-', i, '.rda'))
}
