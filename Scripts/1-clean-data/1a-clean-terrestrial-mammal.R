################
# Combine/clean mammal range data
#
################


# specify dir -------------------------------------------------------------

#path KK machine
dir <- "/mnt/research/ibeem/data/L0/range/MAMMALS_TERRESTRIAL_ONLY.shp"

# load packages -----------------------------------------------------------

library(tidyverse)
library(sf)


# read in data ------------------------------------------------------------

mam <- st_read(dir)


# combine -----------------------------------------------------------------


# Select known areas of occurrence 
# Select only resident spp
# Select native spp 
# See IUCN_Standard_attributes_for_spatial_data-v1.19_2021.xlsx for 
# "IUCN Codes" tab for other code meanings

mam1 <- mam %>% filter(presence == 1, seasonal == 1, origin == 1)

#### Testing #### 
# sci_name and id_no both unique identifiers for each species 
# test <- unique(mam[,c("sci_name", "id_no")])
# length(unique(test$sci_name))
# length(unique(test$id_no))
#################

# Merge species with multiple range polys into one multipoly 
# i.e. one row per species 
# A lot of the duplicate rows seem to be spp that live on multiple islands
mam2 <- mam1 %>% 
  rename(order = order_) %>% 
  group_by(sci_name, id_no, kingdom, phylum, class, order, family, genus) %>% 
  summarize(geometry = st_combine(geometry))


# write out ---------------------------------------------------------------

st_write(mam2, "../terrestrial-mammal/terrestrial-mammals-clean.shp")

          