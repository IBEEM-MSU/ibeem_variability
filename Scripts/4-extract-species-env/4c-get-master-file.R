# 4c-get-master-file.R: script to generate a single master file where each 
#                       row is a species and the columns contain the traits, 
#                       and climate variables

rm(list = ls())

# load packages ------------------------------------------------------------

library(tidyverse)


# Specify directories -----------------------------------------------------

dir <- '/mnt/research/ibeem/variability/'
# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


# Get the climate data ----------------------------------------------------

# Loads the full set of current BL IDs 
load(paste0(dir, 'data/L1/range/bird-breeding/BL-ids.rda'))
# Load the first chunk of species and climate averages to get number of variables
# Loads object called "env.out"
load(paste0(dir, 'data/L2/range-env-pieces/summarized-data-piece-1'))

#build into list
climate.list <- list()
clim.files <- list.files(paste0(dir, 'data/L2/range-env-pieces/'), 
                         full.names = TRUE)
for (i in 1:length(clim.files)) {
  print(i)
  load(clim.files[i])
  climate.list[[i]] <- env.out
}
# Put it in a data frame
climate.df <- do.call("rbind", climate.list)
# How many species have NA for temp_sd_year? Just 14
sum(apply(dplyr::select(climate.df, temp_sd_year), 
          1, function(a) sum(is.na(a))) > 0)


# Load in trait data ------------------------------------------

# AVONET data - keep only necessary fields
avonet.dat <- read.csv(paste0(dir, 'data/L1/trait/avonet-with-id.csv')) %>%
  dplyr::filter(!is.na(BirdLife_ID)) %>%
  dplyr::select(-Sequence, -Avibase.ID1,
                -Total.individuals, -Female, -Male, -Unknown, 
                -Complete.measures, -Inference, -Traits.inferred, 
                -Reference.species, -Centroid.Latitude, -Centroid.Longitude,
                -Range.Size, -Notes, -Mass.Source, -Mass.Refs.Other) %>%
  dplyr::rename(ID = BirdLife_ID, Accepted_name = Species1,
                Family = Family1, Order = Order1)

# Bird et al data
gen.time.dat <- read.csv(paste0(dir, 'data/L1/trait/bird-et-al-data-with-id.csv')) %>%
  dplyr::filter(!is.na(BirdLife_ID)) %>%
  dplyr::select(BirdLife_ID, Sci_name, GenLength) %>%
  dplyr::rename(ID = BirdLife_ID, Accepted_name = Sci_name)


# combine trait and climate data-----------------------------------------------

#only species with values for mean temp
main.dat <- dplyr::full_join(avonet.dat, gen.time.dat,
                             by = c('ID', 'Accepted_name')) %>%
  dplyr::full_join(climate.df, trait.dat, by = c('ID')) %>%
  # add CV for precip and dhi over species range (spatial CV)
  dplyr::mutate(precip_cv_space = precip_sd_space / precip_mean,
                dhi_cum_cv_space = dhi_cum_sd_space / dhi_cum_mean) %>%
  dplyr::filter(!is.na(temp_mean)) %>%
  dplyr::relocate(precip_cv_space, .after = precip_sd_space) %>%
  dplyr::relocate(dhi_cum_cv_space, .after = dhi_cum_sd_space) #%>%
  # dplyr::filter(!is.na(range_size_km2))
 
# which species are missing?
# Antarctic-only seabirds and remote island birds
#dplyr::filter(main.dat, is.na(temp_mean))

#write to file
write.csv(main.dat, file = paste0(dir, 'data/L2/main-bird-data.csv'), 
          row.names = FALSE)

# write.csv(main.dat, file = paste0('~/main-bird-data.csv'),
#           row.names = FALSE)
