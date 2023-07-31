# 4c-get-master-file.R: script to generate a single master file where each 
#                       row is a species and the columns contain the traits, 
#                       and average climate variables
rm(list = ls())
library(tidyverse)

# Specify directories -----------------------------------------------------
climate.dir <- '/mnt/research/ibeem/data/L2/range-env-pieces/'
BL.dir <- '/mnt/research/ibeem/data/L1/range/bird-breeding/'
trait.dir <- '/mnt/research/ibeem/data/L1/trait/'
out.dir <- '/mnt/research/ibeem/data/L2/'

# Get the climate data ----------------------------------------------------
# Loads the full set of current BL IDs 
load(paste0(BL.dir, 'BL-ids.rda'))
# Load the first chunk of species and climate averages to get number of variables
# Loads object called "avg.clim.df"
load(paste0(climate.dir, 'summarized-data-piece-1'))
# Yeah, I know this is super inefficient, but it gets the job done.
climate.list <- list()
clim.files <- list.files(climate.dir)
for (i in 1:length(clim.files)) {
  print(i)
  load(paste0(climate.dir, clim.files[i]))
  climate.list[[i]] <- avg.clim.df
}
# Put it in a data frame
climate.df <- do.call("rbind", climate.list)
# Note that some of the species are waterbirds, and their ranges don't overlap
# with the extracted climate data. 
# How many species are NA? 
sum(apply(climate.df, 1, function(a) sum(is.na(a))) > 0)
# So, still a good amount of species!

# Load in some of the trait data ------------------------------------------
# AVONET data
avonet.dat <- read.csv(paste0(trait.dir, 'avonet-with-id.csv'))
# Bird et al data
gen.time.dat <- read.csv(paste0(trait.dir, 'bird-et-al-data-with-id.csv'))
# Get rid of species without an id and only grab columns of relevance
gen.time.dat <- gen.time.dat %>%
  dplyr::filter(!is.na(id)) %>%
  dplyr::select(starts_with('Measured'), starts_with('Modeled'), GenLength, accepted_name, id) %>%
  setNames(paste0('bird.et.al.', names(.))) %>%
  rename(id = bird.et.al.id, accepted_name = bird.et.al.accepted_name)
avonet.dat <- avonet.dat %>%
  dplyr::filter(!is.na(id)) %>%
  dplyr::select(-Sequence, -Species1, -Family1, -Order1, -Avibase.ID1) %>%
  setNames(paste0('avonet.', names(.))) %>%
  rename(id = avonet.id, accepted_name = avonet.accepted_name)
# Join the two trait data sources
trait.dat <- full_join(gen.time.dat, avonet.dat, by = c('id', 'accepted_name'))

# Put it all together -----------------------------------------------------
main.dat <- full_join(climate.df, trait.dat, by = c('id'))
write.csv(main.dat, file = paste0(out.dir, 'main-bird-data.csv'), row.names = FALSE)

