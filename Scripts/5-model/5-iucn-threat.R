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
bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data-birdtree2.csv')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::filter(Order %ni% or_excl,
                Migration == 1) %>%
  #filter out precip outlier
  dplyr::filter(precip_cv_season < 2.5) %>%
  dplyr::mutate(lMass = log(Mass),
                lGL = log(GenLength),
                lAb = log(Modeled_age_first_breeding),
                lMl = log(Modeled_max_longevity),
                lCs = log(Measured_clutch_size),
                S = Modeled_survival,
                # lMl = log(Measured_max_longevity),
                species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  #drop duplicated species (for now)
  dplyr::group_by(Birdtree_name) %>%
  dplyr::slice_head() %>%
  dplyr::ungroup()

#subset out just traits of interest
bird_df2 <- dplyr::select(bird_df, species,
                          lMass,
                          lGL,
                          S,
                          Measured_survival,
                          Measured_max_longevity,
                          Measured_age_first_breeding,
                          lAb,
                          lMl,
                          lCs,
                          Trophic_niche = Trophic.Niche,
                          Order,
                          Family,
                          cen_lon,
                          cen_lat,
                          temp_mean,
                          temp_sd_year,
                          temp_sd_season,
                          precip_cv_year,
                          precip_cv_season)

bd2 <- dplyr::filter(bird_df, iucn != 'DD', !is.na(iucn))
bd3 <- dplyr::filter(bd2, range_size_km2 > 0)
bd3$dh <- bd3$temp_rel_slope * exp(bd3$lGL)
bd3$temp_rel_slope_season <- bd3$temp_slope / bd3$temp_sd_season
bd3$is <- ordered(bd3$iucn, levels = c('LC', 'NT', 'VU', 'EN', 'CR', 'EW', 'EX'))


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
car::vif(f5)


# SEM ---------------------------------------------------------------------

library(lavaan)

bd4 <- dplyr::mutate(bd3,
                     iv = as.numeric(bd3$is),
                     sc_temp_sd_year = scale(bd3$temp_sd_year, 
                                            scale = TRUE),
                     sc_temp_rel_slope = scale(bd3$temp_rel_slope, 
                                               scale = TRUE),
                     sc_log_range_size_km2 = scale(log(bd3$range_size_km2), 
                                                   scale = TRUE))

model1 <- '
  iv ~ lGL + sc_temp_sd_year + sc_temp_rel_slope + sc_log_range_size_km2
  sc_temp_rel_slope ~ sc_temp_sd_year
  lGL ~ sc_temp_sd_year'

model1.fit <- lavaan::sem(model1, data = bd4) 
summary(model1.fit)

