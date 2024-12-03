# TITLE:            Generate master data file 
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Combine environmental 
# DATA OUTPUT:      Master bird ID list (rda), species-climate data from 3b, trait data (Bird et al. and AVONET) from 1d
# DATE:             September 2023 
# OVERVIEW:         Script to generate a single master file where each row is a species and the columns contain the traits, and climate variables. This file generates the file using the birdtree taxonomy


rm(list = ls())


# load environment variables ------------------------------------------------

source('./Scripts/0-config.R')


# load packages ------------------------------------------------------------

library(tidyverse)


# Get the climate data ----------------------------------------------------

# Loads the full set of current BL IDs
load(paste0(dir, 'data/L1/range/bird-breeding/birdTree-ids.rda'))
# Load the first chunk of species and climate averages to get number of variables
# Loads object called "env.out"
load(paste0(dir, 'data/L2/range-env-pieces-birdtree/summarized-data-piece-BT-1'))

#build into list
climate.list <- list()
clim.files <- list.files(paste0(dir, 'data/L2/range-env-pieces-birdtree/'),
                         full.names = TRUE)
for (i in 1:length(clim.files))
{
  print(i)
  load(clim.files[i])
  climate.list[[i]] <- env.out
}
# Put it in a data frame
climate.df <- do.call("rbind", climate.list)


# Load in trait data ------------------------------------------

# AVONET data - keep only necessary fields
avonet.dat1 <- read.csv(paste0(dir, 'data/L1/trait/avonet-with-id.csv')) %>%
  dplyr::filter(!is.na(birdtreeID)) %>%
  dplyr::select(-Sequence, -Avibase.ID1,
                -Total.individuals, -Female, -Male, -Unknown,
                -Complete.measures, -Inference, -Traits.inferred,
                -Reference.species, -Centroid.Longitude,
                -Range.Size, -birdlifeID, -notes, -Mass.Source, -Mass.Refs.Other,
                -matchType) %>%
  dplyr::rename(ID = birdtreeID, Avonet_name = Species1,
                Birdtree_name = birdtreeName, Family = Family1, Order = Order1)

avonet.dat2 <- avonet.dat1 %>% 
  dplyr::group_by(ID) %>%
  # Get means of continuous variables and the unique categorical variable for those
  # BirdTree species that are split into multiple BirdLife species.
  # Migration and Trophic Niche same across species that have same ID
  dplyr::summarize(Migration = min(Migration),
                   Trophic.Niche = min(Trophic.Niche),
                   Mass = mean(Mass),
                   Lat = mean(Centroid.Latitude, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  #first Family/Order of each ID group (should be the same across all sp ID)
  dplyr::left_join(dplyr::select(avonet.dat1, ID, Avonet_name, Birdtree_name,
                                 Family, Order) %>%
                     dplyr::group_by(ID) %>%
                     dplyr::slice_head(), 
                   by = 'ID') %>%
  dplyr::relocate(Birdtree_name, Avonet_name, Family, Order, .after = ID)

# Bird et al data
gen.time.dat <- read.csv(paste0(dir, 'data/L1/trait/bird-et-al-data-with-id.csv')) %>%
  dplyr::filter(!is.na(birdtreeID)) %>%
  dplyr::select(birdtreeID,
                GenLength) %>%
  dplyr::rename(ID = birdtreeID) %>%
  #average genlength across dups
  dplyr::group_by(ID) %>%
  dplyr::summarize(GenLength = mean(GenLength)) %>%
  dplyr::ungroup()


# combine trait and climate data-----------------------------------------------

#only species with values for mean temp
main.dat <- dplyr::full_join(avonet.dat2, gen.time.dat,
                             by = 'ID') %>%
  dplyr::full_join(climate.df, by = 'ID') %>%
  dplyr::filter(!is.na(GenLength)) %>%
  # add CV for precip over species range (spatial CV)
  dplyr::mutate(precip_mean_cv_space = precip_mean_sd_space / precip_mean) %>%
  # must have at least mean temp, gen length, and Migration
  # dplyr::filter(!is.na(temp_mean), !is.na(GenLength), !is.na(Migration)) %>%
  dplyr::relocate(precip_mean_cv_space, .after = precip_mean_sd_space) %>%
  dplyr::select(-temp_rng_space,-precip_rng_space, -cen_lon,-cen_lat)

# write to file
write.csv(main.dat, file = paste0(dir, 'data/L3/main-bird-data-birdtree2.csv'),
          row.names = FALSE)
