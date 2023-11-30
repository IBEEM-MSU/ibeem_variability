####################
# IUCN status ~ gen length + env change + range size
####################


# specify dir -------------------------------------------------------------

# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# sc_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
dir <- '/mnt/research/ibeem/variability/'
sc_dir <- '/mnt/home/ccy/variability/'
run_date <- '2023-10-27'


# load packages -----------------------------------------------------------

library(tidyverse)
library(MASS)


# load bird data ---------------------------------------------------------------

or_excl <- c('Sphenisciformes', #penguins 
             'Procellariiformes', #tubenoses
             'Pelecaniformes', #pelicans
             'Suliformes', #gannets/boobies
             'Phaethontiformes', #tropicbirds
             'Charadriiformes')#, #skuas, gulls, terns, skimmers, auks
#'Anseriformes', #waterfowl
#'Ciconiiformes', #storks
#'Gaviiformes', #aquatic birds (loons and divers)
#'Gruiformes', #cranes, rails - Family Psophiidae are not waterbirds, but there are few members (~6 species)
#'Phoenicopteriformes', #flamingos and relatives
#'Podicipediformes') #grebes

'%ni%' <- Negate('%in%')
#subset out just traits of interest
#calculate temp change by generation and relative temp change by generation
bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data-birdtree2.csv')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::filter(Order %ni% or_excl,
                Migration == 1,
                range_size_km2 > 0,
                iucn != 'DD',
                !is.na(iucn)) %>%
  dplyr::mutate(species = stringr::str_to_title(gsub(' ', '_', Birdtree_name)),
                iucn_level = ordered(iucn, levels = c('LC', 'NT', 'VU', 'EN', 'CR', 'EW', 'EX'))) %>%
  dplyr::select(species,
                Birdtree_name,
                Mass,
                GenLength,
                range_size_km2,
                temp_slope,
                temp_sd_year,
                temp_rel_slope,
                iucn_level) %>%
  #keep duplicated species (for now)
  # dplyr::group_by(Birdtree_name) %>%
  # dplyr::slice_head() %>%
  # dplyr::ungroup() %>%
  dplyr::mutate(temp_slope_gen = temp_slope*GenLength,
                temp_rel_slope_gen = temp_rel_slope*GenLength)

# run ordinal regressions ---------------------------------------------------------------

plot(bird_df$iucn_level, log(bird_df$range_size_km2))
f1 <- MASS::polr(iucn_level ~ scale(log(Mass)) + scale(log(GenLength)) + scale(log(range_size_km2)) + 
                   scale(temp_slope), data=bird_df)
f2 <- MASS::polr(iucn_level ~ scale(log(Mass)) + scale(log(GenLength)) + scale(log(range_size_km2)) + 
                   scale(temp_rel_slope), data=bird_df)
f3 <- MASS::polr(iucn_level ~ scale(log(Mass)) + scale(log(GenLength)) + scale(log(range_size_km2)) +
                   scale(temp_slope_gen), data=bird_df)
f4 <- MASS::polr(iucn_level ~ scale(log(Mass)) + scale(log(GenLength)) + scale(log(range_size_km2)) +
                   scale(temp_rel_slope_gen), data=bird_df)
summary(f1)
summary(f2)
summary(f3)
summary(f4)

# Check temp_sd_year
f0 <- MASS::polr(iucn_level ~ scale(log(Mass)) + scale(log(GenLength)) + scale(log(range_size_km2)) +
                   scale(temp_sd_year), data=bird_df)
summary(f0)

# run logistic regressions ---------------------------------------------------------------

bird_df2 <- bird_df %>%
  dplyr::filter(iucn_level %in% c("LC","NT","VU","EN","CR")) %>%
  dplyr::mutate(threatened = ifelse(iucn_level %in% c("LC","NT"),0,1))
f1 <- glm(threatened ~ scale(log(Mass)) + scale(log(GenLength)) + scale(log(range_size_km2)) +
                   scale(temp_slope), data=bird_df2, family = "binomial")
f2 <- glm(threatened ~ scale(log(Mass)) + scale(log(GenLength)) + scale(log(range_size_km2)) +
            scale(temp_rel_slope), data=bird_df2, family = "binomial")
f3 <- glm(threatened ~ scale(log(Mass)) + scale(log(GenLength)) + scale(log(range_size_km2)) +
            scale(temp_slope_gen), data=bird_df2, family = "binomial")
f4 <- glm(threatened ~ scale(log(Mass)) + scale(log(GenLength)) + scale(log(range_size_km2)) +
            scale(temp_rel_slope_gen), data=bird_df2, family = "binomial")
summary(f1)
summary(f2)
summary(f3)
summary(f4)

# Check temp_sd_year
f0 <- glm(threatened ~ scale(log(Mass)) + scale(log(GenLength)) + scale(log(range_size_km2)) +
            scale(temp_sd_year), data=bird_df2, family = "binomial")
summary(f0)

################# earlier code ###############

#more relative temp change, more likely to be endangered
plot(bd3$is, bd3$temp_rel_slope)
plot(bd3$is, log(bd3$range_size_km2))
f1 <- MASS::polr(bd3$is ~ bd3$temp_rel_slope + log(bd3$range_size_km2))
summary(f1)
#same as above, though accounting for range size (smaller ranges more likely to be endangered)
f2 <- MASS::polr(bd3$is ~ bd3$temp_slope + log(bd3$range_size_km2))
summary(f2)
#more delta haldane, more likely to be endangered
f3 <- MASS::polr(bd3$is ~ scale(bd3$dh, scale = TRUE) + log(bd3$range_size_km2))
summary(f3)
#better AIC, longer gens and more rapid temp change, more likely to be endangered
#range size used in RL assessment
#GL is used in the sense that it's threshold for pop decline in 10 years or 3 gens, whichever is longer
# http://datazone.birdlife.org/species/spcredcrit
f4 <- MASS::polr(bd3$is ~ scale(bd3$lGL, scale = TRUE) + 
             scale(bd3$temp_rel_slope, scale = TRUE) + 
             log(bd3$range_size_km2))
summary(f4)
exp(coef(f4))

f5 <- MASS::polr(bd3$is ~ scale(bd3$lGL, scale = TRUE) + 
                   # scale(bd3$S, scale = TRUE) +
                   # scale(bd3$lAb, scale = TRUE) +
                   # scale(bd3$lMl, scale = TRUE) +
                   scale(bd3$temp_sd_year, scale = TRUE) +
                   # scale(bd3$temp_sd_season, scale = TRUE) +
                   # scale(bd3$precip_cv_year, scale = TRUE) +
                   # scale(bd3$precip_cv_season, scale = TRUE) +
                   scale(bd3$temp_rel_slope, scale = TRUE) +
                   # scale(bd3$temp_slope, scale = TRUE) +
                   log(bd3$range_size_km2))
summary(f5)
car::vif(f2)
