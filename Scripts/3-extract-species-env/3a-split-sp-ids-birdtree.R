# TITLE:            Generate bird data processing batches 
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Rda file with IDs of all birds in analysis
# DATA OUTPUT:      Rda files with lists of IDs for batch processing
# DATE:             May 2023 
# OVERVIEW:         Split the saved rda file of BirdTree IDs into multiple files for easily running multiple jobs.


rm(list = ls())


# load environment variables ------------------------------------------------

source("./Scripts/0-config.R")


# Specify directory -------------------------------------------------------

BT.dir <- paste0(dir, '/data/L1/range/bird-breeding/')


# Load in full set of current IDs -----------------------------------------

load(paste0(BT.dir, 'birdTree-ids.rda'))
# Remove NA id from birdTree.unique.ids
BT.unique.ids <- birdTree.unique.ids[!is.na(birdTree.unique.ids)]
n.ids <- length(BT.unique.ids)
n.pieces <- 10000
vals <- split(BT.unique.ids, ceiling(seq_along(1:n.ids)/n.pieces))

for (i in 1:length(vals)) {
  ids <- vals[[i]]
  save(ids, file = paste0(BT.dir, 'BTIDsPiece-', i, '.rda'))
}
