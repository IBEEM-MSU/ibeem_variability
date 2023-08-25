# 4c-get-master-file.R: script to generate a single master file where each 
#                       row is a species and the columns contain the traits, 
#                       and average climate variables

rm(list = ls())
library(tidyverse)


# Specify directories -----------------------------------------------------

climate.dir <- '/mnt/research/ibeem/variability/data/L2/range-env-pieces/'
BL.dir <- '/mnt/research/ibeem/variability/data/L1/range/bird-breeding/'
trait.dir <- '/mnt/research/ibeem/variability/data/L1/trait/'
out.dir <- '/mnt/research/ibeem/variability/data/L2/'


# Get the climate data ----------------------------------------------------

# Loads the full set of current BL IDs 
load(paste0(BL.dir, 'BL-ids.rda'))
# Load the first chunk of species and climate averages to get number of variables
# Loads object called "env.out"
load(paste0(climate.dir, 'summarized-data-piece-1'))

#build into list
climate.list <- list()
clim.files <- list.files(climate.dir)
for (i in 1:length(clim.files)) {
  print(i)
  load(paste0(climate.dir, clim.files[i]))
  climate.list[[i]] <- env.out
}
# Put it in a data frame
climate.df <- do.call("rbind", climate.list)
# Some species don't overlap climate data (e.g., Antarctic-only seabirds)
# How many species are NA? Just 16
sum(apply(climate.df, 1, function(a) sum(is.na(a))) > 0)


# Load in trait data ------------------------------------------

# AVONET data - keep only necessary fields
avonet.dat <- read.csv(paste0(trait.dir, 'avonet-with-id.csv')) %>%
  dplyr::filter(!is.na(BirdLife_ID)) %>%
  dplyr::select(-Sequence, -Avibase.ID1,
                -Total.individuals, -Female, -Male, -Unknown, 
                -Complete.measures, -Inference, -Traits.inferred, 
                -Reference.species, -Centroid.Latitude, -Centroid.Longitude,
                -Range.Size, -Notes, -Mass.Source, -Mass.Refs.Other) %>%
  dplyr::rename(ID = BirdLife_ID, Accepted_name = Species1,
                Family = Family1, Order = Order1)

# Bird et al data
gen.time.dat <- read.csv(paste0(trait.dir, 'bird-et-al-data-with-id.csv')) %>%
  dplyr::filter(!is.na(BirdLife_ID)) %>%
  dplyr::select(BirdLife_ID, Sci_name, GenLength) %>%
  dplyr::rename(ID = BirdLife_ID, Accepted_name = Sci_name)


# combine trait and climate data-----------------------------------------------

#only species with values for mean temp
main.dat <- dplyr::full_join(avonet.dat, gen.time.dat,
                             by = c('ID', 'Accepted_name')) %>%
  dplyr::full_join(climate.df, trait.dat, by = c('ID')) %>%
  dplyr::filter(!is.na(temp_mean))

#XXXX
#mutate precip_cv_space


#write to file
write.csv(main.dat, file = paste0(out.dir, 'main-bird-data.csv'), 
          row.names = FALSE)

