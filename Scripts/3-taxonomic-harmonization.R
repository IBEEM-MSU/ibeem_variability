################
# Resolve differences in species names and combine bird/mammal datasets:
# - range maps
# - generation time
# - other traits
# - phylogeny
# 
# NOTES:
# - should decide on taxonomic source to match to
# - create taxonomic key to match to from each source?
################


# specify dir -------------------------------------------------------------

#path for data on CY machine - remember trailing slash
range_map_data_dir <- '~/Downloads/BOTW_2022_2/'
life_history_dir <- '~/Downloads/'

#directory to save out intermediate file
out_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/Sample_output/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)
library(taxize)


# read in data -------------------------------------------------

#Birdlife range maps - takes a few minutes (10 min?) to read in
#data location: https://drive.google.com/drive/u/1/folders/11eAFmFKUc7tU59IArNEpN5YP_ehGvAwZ
BL_data <- sf::st_read(dsn = paste0(range_map_data_dir, 'BOTW.gdb/a0000000a.gdbtable'))

#bird life history data (processed to L1)
#data location: https://drive.google.com/drive/u/1/folders/18p2Zn3dMA78hdwrCVQBtJ8dnnLEYa7_w
LH_data <- read.csv(paste0(life_history_dir, 'Bird_et_al_gen_length_birds.csv'))

#OR read in species names - SEE BELOW


# process data ------------------------------------------------------------

#IF ITNERMEDIATE OBJECTS EXIST - read in
# BL_usp <- readRDS(paste0(out_dir, 'BL_names.rds'))
# LH_usp <- readRDS(paste0(out_dir, 'LH_names.rds'))

#pretty close to same number of species...positive sign
BL_usp <- unique(BL_data$sci_name)
LH_usp <- unique(LH_data$Sci_name)
length(BL_usp)
length(LH_usp)

#clear up memory
rm(BL_data)
rm(LH_data)
gc()

# #save out intermediate objects - to work with without loading in range maps
# saveRDS(BL_usp, paste0(out_dir, 'BL_names.rds'))
# saveRDS(LH_usp, paste0(out_dir, 'LH_names.rds'))


#differences in the hundreds
#which species in BL are not in LH
'%ni%' <- Negate('%in%')
BL_usp[which(BL_usp %ni% LH_usp)]

#which species in LH are not in BL
LH_usp[which(LH_usp %ni% BL_usp)]


# match -------------------------------------------------------------------

BL_test <- sort(BL_usp)[1:50]
LH_test <- sort(LH_usp)[1:50]

# https://docs.ropensci.org/taxize/articles/taxize.html
#get classification and match
BL_class <- taxize::classification(BL_test, db = 'ncbi')
LH_class <- taxize::classification(LH_test, db = 'ncbi')

#find which don't match taxize database
no_match_BL <- names(which(is.na(BL_class)))
no_match_LH <- names(which(is.na(LH_class)))

#resolve names
#if getting error, website might be down: https://github.com/ropensci/taxize/issues/895
taxize::gnr_resolve(no_match_BL)
taxize::gnr_resolve(no_match_LH)

#matches between BL and LH of first 50
length(match(BL_class, LH_class))

