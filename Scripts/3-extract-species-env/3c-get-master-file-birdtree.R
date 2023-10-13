# 3c-get-master-file-birdtree.R: script to generate a single master file where each 
#                                row is a species and the columns contain the traits, 
#                                and climate variables. This file generates the 
#                                file using the birdtree taxonomy. 

rm(list = ls())

# load packages ------------------------------------------------------------

library(tidyverse)


# Specify directories -----------------------------------------------------

dir <- '/mnt/research/ibeem/variability/'
# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'


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
for (i in 1:length(clim.files)) {
  print(i)
  load(clim.files[i])
  climate.list[[i]] <- env.out
}
# Put it in a data frame
climate.df <- do.call("rbind", climate.list)
# How many species have NA for temp_sd_year? Just 11
sum(apply(dplyr::select(climate.df, temp_sd_year), 
          1, function(a) sum(is.na(a))) > 0)


# Load in trait data ------------------------------------------

# AVONET data - keep only necessary fields
avonet.dat1 <- read.csv(paste0(dir, 'data/L1/trait/avonet-with-id.csv')) %>%
  dplyr::filter(!is.na(birdtreeID)) %>%
  dplyr::select(-Sequence, -Avibase.ID1,
                -Total.individuals, -Female, -Male, -Unknown, 
                -Complete.measures, -Inference, -Traits.inferred, 
                -Reference.species, -Centroid.Latitude, -Centroid.Longitude,
                -Range.Size, -birdlifeID, -notes, -Mass.Source, -Mass.Refs.Other, 
		-matchType) %>%
  dplyr::rename(ID = birdtreeID, Avonet_name = Species1, 
                Birdtree_name = birdtreeName, Family = Family1, Order = Order1)

avonet.dat2 <- dplyr::group_by(avonet.dat1, ID) %>%
  # Get means of continuous variables and the unique categorical variable for those
  # BirdTree species that are split into multiple BirdLife species.
  dplyr::summarize(Beak.Length_Culmen = mean(Beak.Length_Culmen),
		   Beak.Length_Nares = mean(Beak.Length_Nares),
		   Beak.Width = mean(Beak.Width),
		   Beak.Depth = mean(Beak.Depth),
		   Tarsus.Length = mean(Tarsus.Length),
		   Wing.Length = mean(Wing.Length),
		   Kipps.Distance = mean(Kipps.Distance),
		   Secondary1 = mean(Secondary1),
		   Hand.Wing.Index = mean(Hand.Wing.Index),
		   Tail.Length = mean(Tail.Length),
		   Mass = mean(Mass),
		   Habitat = unique(Habitat),
		   Habitat.Density = unique(Habitat.Density),
		   Migration = unique(Migration),
		   Trophic.Level = unique(Trophic.Level),
		   Trophic.Niche = unique(Trophic.Niche),
		   Primary.Lifestyle = unique(Primary.Lifestyle),
		   Min.Latitude = min(Min.Latitude),
		   Max.Latitude = max(Max.Latitude)) %>%
  dplyr::ungroup() %>%
  #first Family/Order of each ID group (should be the same across all sp ID)
  dplyr::left_join(dplyr::select(avonet.dat1, ID, Avonet_name, Birdtree_name,
                                 Family, Order) %>%
                     dplyr::group_by(ID) %>%
                     dplyr::slice_head()) %>%
  dplyr::relocate(Birdtree_name, Avonet_name, Family, Order, .after = ID)
  

# Bird et al data
gen.time.dat <- read.csv(paste0(dir, 'data/L1/trait/bird-et-al-data-with-id.csv')) %>%
  dplyr::filter(!is.na(birdtreeID)) %>%
  dplyr::select(birdtreeID, 
                GenLength,
                Measured_survival,
                Measured_age_first_breeding,
                Measured_max_longevity,
                Modeled_survival,
                Modeled_age_first_breeding,
                Modeled_max_longevity) %>%
  #correct for incorrect scaling in data
  dplyr::mutate(Measured_survival = Measured_survival * 0.01) %>%
  dplyr::rename(ID = birdtreeID) %>%
  #average genlength across dups
  dplyr::group_by(ID) %>%
  dplyr::summarize(GenLength = mean(GenLength),
                   Measured_survival = mean(Measured_survival, 
                                            na.rm = TRUE),
                   Measured_age_first_breeding = mean(Measured_age_first_breeding, 
                                                      na.rm = TRUE),
                   Measured_max_longevity = mean(Measured_max_longevity, 
                                                 na.rm = TRUE),
                   Modeled_survival = mean(Modeled_survival, 
                                           na.rm = TRUE),
                   Modeled_age_first_breeding = mean(Modeled_age_first_breeding, 
                                                     na.rm = TRUE),
                   Modeled_max_longevity = mean(Modeled_max_longevity, 
                                                na.rm = TRUE)) %>%
  dplyr::ungroup()


# combine trait and climate data-----------------------------------------------

#only species with values for mean temp
main.dat <- dplyr::full_join(avonet.dat2, gen.time.dat,
                             by = 'ID') %>%
  dplyr::full_join(climate.df, by = 'ID') %>%
  # add CV for precip and dhi over species range (spatial CV)
  dplyr::mutate(precip_cv_space = precip_sd_space / precip_mean,
                dhi_cum_cv_space = dhi_cum_sd_space / dhi_cum_mean) %>%
  # must have at least mean temp, gen length, and Migration
  dplyr::filter(!is.na(temp_mean), !is.na(GenLength), !is.na(Migration)) %>%
  dplyr::relocate(precip_cv_space, .after = precip_sd_space) %>%
  dplyr::relocate(dhi_cum_cv_space, .after = dhi_cum_sd_space) #%>%
  # dplyr::filter(!is.na(range_size_km2))


# dplyr::filter(main.dat, is.na(GenLength) | 
#                 is.na(temp_mean) | 
#                 is.na(Migration)) %>%
#   dplyr::select(Accepted_name, ID, 
#                 GenLength, temp_mean, Migration)


#write to file
write.csv(main.dat, file = paste0(dir, 'data/L3/main-bird-data-birdtree2.csv'), 
          row.names = FALSE)

# write.csv(main.dat, file = paste0('~/main-bird-data.csv'),
#           row.names = FALSE)
