# 4c-get-master-file.R: script to generate a single master file where each 
#                       row is a species and the columns contain the traits, 
#                       and average climate variables
rm(list = ls())
library(dplyr)

# Specify directories -----------------------------------------------------
climate.dir <- '/mnt/research/ibeem/data/L2/range-env-pieces-mammal/'
MAM.dir <- '/mnt/research/ibeem/data/L1/range-mammal/'
trait.dir <- '/mnt/research/ibeem/data/L1/trait-mammal/'
# trait.dir <- '../trait/'
out.dir <- '/mnt/research/ibeem/data/L2/'

# Get the climate data ----------------------------------------------------
# Loads the full set of current MAM IDs 
load(paste0(MAM.dir, 'MAM-ids.rda'))
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
climate.df <- do.call("rbind", climate.list) %>% 
  mutate(id = paste0("ITIS:", id))

# How many species are NA? 
sum(apply(climate.df, 1, function(a) sum(is.na(a))) > 0)
# So, still a good amount of species!

# Load in some of the trait data ------------------------------------------
# Pacifici et al data
gen.time.dat <- read.csv(paste0(trait.dir, 'pacifici-traits-with-id.csv'))
# Get rid of species without an id and only grab columns of relevance
trait.dat <- gen.time.dat %>%
  dplyr::filter(!is.na(id)) %>%
  dplyr::select(-TaxID, -Order, -Family, -Genus, -Scientific_name) %>%
  setNames(paste0('pacifici.et.al.', names(.))) %>%
  rename(id = pacifici.et.al.id, accepted_name = pacifici.et.al.accepted_name)

# Put it all together -----------------------------------------------------
main.dat <- full_join(climate.df, trait.dat, by = c("id"))
write.csv(main.dat, file = paste0(out.dir, 'main-mammal-data.csv'), row.names = FALSE)

