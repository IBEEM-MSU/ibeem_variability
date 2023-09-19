################
# Initial exploration mammals
#
#
################


# Specify dir --------------------------------------------------

#path CY machine
dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(viridis)


# read in data -------------------------------------------------

mam_df <- read.csv(paste0(dir, 'Data/L2/main-mammal-data.csv')) %>%
  dplyr::mutate(fac_Family = factor(Family),
                fac_Order = factor(Order),
                lMass = log(Mass),
                lGL = log(GenLength))

# or_excl <- c('Sphenisciformes', #penguins 
#              'Procellariiformes', #tubenoses
#              'Pelecaniformes', #pelicans
#              'Suliformes', #gannets/boobies
#              'Phaethontiformes', #tropicbirds
#              'Charadriiformes', #skuas, gulls, terns, skimmers, auks
#              'Anseriformes', #waterfowl
#              'Ciconiiformes', #storks
#              'Gaviiformes', #aquatic birds (loons and divers)
#              'Gruiformes', #cranes, rails - Family Psophiidae are not waterbirds, but there are few members (~6 species)
#              'Phoenicopteriformes', #flamingos and relatives
#              'Podicipediformes') #grebes
# 
# bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data.csv')) %>%
#   dplyr::arrange(Accepted_name) %>%
#   dplyr::filter(Order %ni% or_excl,
#                 Migration == 1)
