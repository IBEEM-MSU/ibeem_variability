####################
# Fit Bayes model - clutch size ~ env + phylo (max like)
####################


# specify dir -------------------------------------------------------------

# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# sc_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
dir <- '/mnt/research/ibeem/variability/'
sc_dir <- '/mnt/home/ccy/variability/'
run_date <- '2023-10-12'


# load packages -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(MCMCvis)
library(ape)
library(Rphylopars)


# load bird data ---------------------------------------------------------------

#laod bird data, merge with surv and clutch size data
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
                lGL = log(GenLength)) %>%
  #drop duplicated species (for now)
  dplyr::group_by(Birdtree_name) %>%
  dplyr::slice_head() %>%
  dplyr::ungroup()

#clutch size from Bird et al.
bcs <- read.csv(paste0(dir, 'data/L1/trait/bird_et_al_clutch_size.csv'))

#join
bird_df2 <- dplyr::mutate(bird_df, Scientific.name = stringr::str_to_sentence(Birdtree_name)) %>%
  dplyr::left_join(dplyr::select(bcs, -Order, -Family, -Genus), 
                   by = 'Scientific.name') %>%
  dplyr::select(ID, Birdtree_name, Avonet_name, Family, Order, Trophic.Niche,
                temp_mean, temp_sd_year, temp_sd_season, temp_sp_color_month,
                precip_mean, precip_cv_year, precip_cv_season, precip_sp_color_month,
                Mass, GenLength, lMass, lGL, Mean.clutch.size, Measured_survival, 
                Measured_age_first_breeding, Measured_max_longevity, 
                Modeled_survival, Modeled_age_first_breeding, Modeled_max_longevity) %>%
  dplyr::mutate(species = stringr::str_to_title(gsub(' ', '_', Birdtree_name)))

#subset out just traits of interest for phylo imputation
bird_df3 <- dplyr::mutate(bird_df2, Measured_log_age_first_breeding = log(Measured_age_first_breeding),
                     Measured_log_max_longevity = log(Measured_max_longevity),
                     Measured_log_clutch_size = log(Mean.clutch.size)) %>%
        dplyr::filter(!is.na(Measured_log_clutch_size))
                

#load consensus tree - bird.phylo
load(paste0(dir, 'data/L3/bird-consensus-tree.rda'))

#prune tree
#df with names and idx
idx_df2 <- data.frame(idx = 1:NROW(bird_df3),
                      name = stringr::str_to_title(gsub(' ', '_', bird_df3$species)))

#species not found in both datasets (species to drop from tree)
nm <- setdiff(bird.phylo$tip.label, bird_df3$species)

#prune specified tips from tree
pr_tree <- ape::drop.tip(bird.phylo, nm)

#get idx
j_idx3 <- dplyr::left_join(data.frame(species = pr_tree$tip.label), 
                           data.frame(idx = 1:NROW(bird_df3), bird_df3), 
                           by = 'species')

#apply
bird_df4 <- bird_df3[j_idx3$idx,]

#make tree binary
pr_tree2 <- ape::multi2di(pr_tree)

#make response into matrix with species as rownames
dd <- dplyr::select(bird_df4, 
                    Measured_log_clutch_size) %>%
  as.matrix()
row.names(dd) <- bird_df4$species

#get estimate of Pagel's kappa to scale phylogeny
fit_ka <- geiger::fitContinuous(pr_tree2, dd[,'Measured_log_clutch_size'], model = "kappa")

#rescale tree
pr_tree_k <- phytools::rescale(pr_tree, 'kappa', 
                     kappa = fit_ka$opt$kappa, sigsq = fit_ka$opt$sigsq)

#get corr matrix of rescaled tree
Rho <- ape::vcv.phylo(pr_tree_k, corr = TRUE)

CM <- nlme::corSymm(Rho[lower.tri(Rho)], fixed = T)


pgls_fit <- nlme::gls(Measured_log_clutch_size ~ lMass + temp_sd_season + temp_sd_year + precip_cv_season + precip_cv_year, correlation = CM,
data = bird_df4,
method = "REML")

summary(pgls_fit)
