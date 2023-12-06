####################
# Model residuals - no clear axes of unexplained variation
####################


# specify dir -------------------------------------------------------------

dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
sc_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# dir <- '/mnt/research/ibeem/variability/'
# sc_dir <- '/mnt/home/ccy/variability/'


# load packages -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(MCMCvis)


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
                mlAb = log(Measured_age_first_breeding),
                lMl = log(Modeled_max_longevity),
                mlMl = log(Measured_max_longevity),
                species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  #drop duplicated species (for now)
  dplyr::group_by(Birdtree_name) %>%
  dplyr::slice_head() %>%
  dplyr::ungroup()

#subset out just traits of interest
bird_df2 <- dplyr::select(bird_df, 
                          species,
                          Order,
                          Family,
                          lGL,
                          # S = Modeled_survival,
                          # lAb,
                          # lMl,
                          # mS = Measured_survival,
                          # mlAb,
                          # mlMl,
                          Trophic_niche = Trophic.Niche,
                          lMass,
                          temp_sd_year,
                          temp_sd_season,
                          precip_cv_year,
                          precip_cv_season) %>% 
  dplyr::left_join(read.csv(paste0(dir, 'data/L0/trait/pnas.2121467120.sd01(2).csv')), 
                 by = c('species' = 'tip_label'))


# load fit model ----------------------------------------------------------

#SURVIVAL
# run_date <- '2023-10-28'
# fit <- readRDS(paste0(dir, '/Results/bird-s-phylo-vint-oe-', run_date,
#                       '/bird-s-phylo-vint-oe-fit-', run_date, '.rds'))
# DATA <- readRDS(paste0(dir, '/Results/bird-s-phylo-vint-oe-', run_date,
#                       '/bird-s-phylo-vint-oe-data-', run_date, '.rds'))

#AGE FIRST BREEDING
run_date <- '2023-10-26'
fit <- readRDS(paste0(dir, '/Results/bird-ab-phylo-vint-oe-', run_date,
                      '/bird-ab-phylo-vint-oe-fit-', run_date, '.rds'))
DATA <- readRDS(paste0(dir, '/Results/bird-ab-phylo-vint-oe-', run_date,
                      '/bird-ab-phylo-vint-oe-data-', run_date, '.rds'))

#MAX LONGEVITY
# run_date <- '2023-10-28'
# fit <- readRDS(paste0(dir, '/Results/bird-ml-phylo-vint-oe-', run_date,
#                       '/bird-ml-phylo-vint-oe-fit-', run_date, '.rds'))
# DATA <- readRDS(paste0(dir, '/Results/bird-ml-phylo-vint-oe-', run_date,
#                       '/bird-ml-phylo-vint-oe-data-', run_date, '.rds'))


# get resids --------------------------------------------------------------

alpha_mn <- MCMCvis::MCMCpstr(fit, params = 'alpha')[[1]]
kappa_mn <- MCMCvis::MCMCpstr(fit, params = 'kappa')[[1]]
gamma_mn <- MCMCvis::MCMCpstr(fit, params = 'gamma')[[1]]
beta_mn <- MCMCvis::MCMCpstr(fit, params = 'beta')[[1]]

#combine linear predictor
mu_mn <- rep(NA, DATA$N)
mu_mn[DATA$obs_idx] <- kappa_mn + gamma_mn[DATA$niche_obs_idx] + 
  alpha_mn[DATA$obs_idx] + (as.matrix(DATA$X_obs) %*% beta_mn)[,1]
mu_mn[DATA$mod_idx] <- kappa_mn + gamma_mn[DATA$niche_mod_idx] +
  alpha_mn[DATA$mod_idx] + (as.matrix(DATA$X_mod) %*% beta_mn)[,1]

#combine reponse
Y_comb <- rep(NA, length(DATA$obs_idx) + length(DATA$mod_idx))
Y_comb[DATA$obs_idx] <- DATA$Y_obs
Y_comb[DATA$mod_idx] <- DATA$Y_mod

#residuals
resid <- Y_comb - mu_mn
bird_df2$resid <- resid


# look at resids ----------------------------------------------------------

plot(factor(bird_df$Habitat), resid)
plot(factor(bird_df$Habitat.Density), resid)
plot(factor(bird_df$Trophic.Level), resid)
plot(bird_df$Hand.Wing.Index, resid)
plot(factor(bird_df$Primary.Lifestyle), resid)
plot(log(bird_df$range_size_km2), resid)
plot(bird_df$cen_lat, resid)
plot(bird_df$cen_lon, resid)
plot(bird_df$temp_skew, resid)
plot(bird_df$Measured_clutch_size, resid)


plot(factor(bird_df2$devo_mode), resid)
plot(log(bird_df2$egg_mass), resid)
#brain size does not explain
bd3 <- dplyr::filter(bird_df2, !is.na(brain))
plot(residuals(lm(log(brain) ~ log(weight), data = bd3)), bd3$resid)
summary(lm(bd3$resid ~ residuals(lm(log(brain) ~ log(weight), data = bd3))))
plot(bird_df2$time_fed, resid)
summary(lm(resid ~ bird_df2$time_fed))
plot(bird_df2$caretakers, resid)
summary(lm(resid ~ bird_df2$caretakers))
plot(factor(bird_df2$social_bonds), resid)
plot(bird_df2$food_energy, resid)
summary(lm(resid ~ bird_df2$food_energy))
plot(bird_df2$food_h_level, resid)
plot(bird_df2$fibres, resid)
summary(lm(resid ~ bird_df2$fibres))
plot(factor(bird_df2$foraging), resid)
plot(factor(bird_df2$sedentariness), resid)
plot(factor(bird_df2$insularity), resid)
plot(factor(bird_df2$colonial), resid)
plot(factor(bird_df2$grouping), resid)

par(mar = c(7, 4, 2, 2) + 1)
plot(factor(bird_df$Order), resid, las = 2)



plot(factor(bd3$Trophic_niche), 
     residuals(lm(log(brain) ~ log(weight), data = bd3)), 
     las = 2)

plot(factor(bd3$Order), 
     residuals(lm(log(brain) ~ log(weight), data = bd3)), 
     las = 2)

ggplot(bd3, aes(log(weight), log(brain), color = Trophic_niche)) + 
  geom_point()

ggplot(bd3, aes(log(weight), log(brain), color = Order)) + 
  geom_point()
