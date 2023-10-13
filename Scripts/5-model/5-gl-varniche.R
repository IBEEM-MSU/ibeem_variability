####################
# Fit Bayes model - gen length ~ env (varying intercepts and slopes by Trophic niche)
# 
####################


# specify dir -------------------------------------------------------------

# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# sc_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
dir <- '/mnt/research/ibeem/variability/'
sc_dir <- '/mnt/home/ccy/variability/'
run_date <- '2023-10-13'


# load packages -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(MCMCvis)
library(ape)


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

#subset out just traits of interest
bird_df2 <- dplyr::mutate(bird_df,
                          species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  dplyr::select(species,
                lMass,
                lGL,
                Modeled_survival,
                Modeled_age_first_breeding,
                Modeled_max_longevity,
                Trophic_niche = Trophic.Niche,
                temp_mean,
                temp_sd_year,
                temp_sd_season,
                temp_sp_color_month,
                precip_cv_year,
                precip_cv_season,
                precip_sp_color_month)


# phylo -------------------------------------------------------------------

#load consensus tree - bird.phylo
load(paste0(dir, 'data/L3/bird-consensus-tree.rda'))

#subset of imp
# set.seed(1)
# Nsel <- 1000
# stidx <- sample(x = 1:NROW(bird_df2), size = Nsel)
# bird_df3 <- bird_df2[stidx,]
bird_df3 <- bird_df2

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
#niche levels
#assign missing species (owls) to Vertivore
bird_df4$Trophic_niche[which(is.na(bird_df4$Trophic_niche))] <- "Vertivore"
bird_df4$niche_idx <- as.numeric(factor(bird_df4$Trophic_niche))
niche_names <- levels(factor(bird_df4$Trophic_niche))


# scale/prep data ---------------------------------------------------------

#scalars for data
lMass_scalar <- 1
temp_sd_season_scalar <- 1
temp_sd_year_scalar <- 1
precip_cv_season_scalar <- 1
precip_cv_year_scalar <- 1
y_scalar <- 1

dd <- dplyr::group_by(bird_df4, niche_idx) %>%
  dplyr::mutate(sc_lMass = scale(lMass, scale = FALSE)[,1],
                sc_temp_sd_season = scale(temp_sd_season, scale = FALSE)[,1],
                sc_temp_sd_year = scale(temp_sd_year, scale = FALSE)[,1],
                sc_precip_cv_season = scale(precip_cv_season, scale = FALSE)[,1],
                sc_precip_cv_year = scale(precip_cv_year, scale = FALSE)[,1])


# Run Stan model --------------------------------------------------------------

#data for model
DATA <- list(N = NROW(bird_df4),
             J = length(unique(bird_df4$niche_idx)),
             y = bird_df4$lGL,
             niche_idx = bird_df4$niche_idx, #niche id for each data point
             lMass = dd$sc_lMass,
             temp_sd_season = dd$sc_temp_sd_season,
             temp_sd_year = dd$sc_temp_sd_year,
             precip_cv_season = dd$sc_precip_cv_season,
             precip_cv_year = dd$sc_precip_cv_year,
             pro_data = bird_df4)

options(mc.cores = parallel::detectCores())

# DELTA <- 0.92
# TREE_DEPTH <- 12
# STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 2000

#compile model
mod <- cmdstanr::cmdstan_model(paste0(sc_dir, 'Scripts/Model_files/5-vint-vsl.stan'))

#sample
fit <- mod$sample(
  data = DATA,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 500)


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fit, 
                  round = 4,
                  file_name = paste0('bird-gl-vint-vsl-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('bird-gl-vint-vsl-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('bird-gl-vint-vsl-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('bird-gl-vint-vsl-data-', run_date),
                  cp_file = c(paste0(sc_dir, 'Scripts/Model_files/5-vint-vsl.stan'), 
                              paste0(sc_dir, 'Scripts/5-model/5-gl-varniche.R')),
                  cp_file_names = c(paste0('5-vint-vsl', run_date, '.stan'),
                                    paste0('5-gl-varniche-', run_date, '.R')))

fig_dir <- paste0(dir, 'Results/bird-gl-vint-vsl-', run_date, '/')

# library(shinystan)
# shinystan::launch_shinystan(fit)
# fit <- readRDS(paste0(dir, '/Results/se-bird-novar-oe-sep-', run_date,
#                       '/se-bird-novar-oe-sep-fit-', run_date, '.rds'))


# phylo signal in resids --------------------------------------------------

#extract resids
mu_mn <- MCMCvis::MCMCpstr(fit, params = 'mu')[[1]]
resids <- DATA$y - mu_mn

#K ~ 0.15
library(phytools)
phytools::phylosig(pr_tree, resids, method = 'K', test = TRUE) #quite slow


# Summary -----------------------------------------------------------------

#model summary
MCMCvis::MCMCsummary(fit, round = 3, 
                     params = c('mu_beta', #effect log mass
                                'mu_gamma1', #effect temp seasonality
                                'mu_gamma2', #effect temp interannual var
                                'mu_theta1', #effect precip seasonality
                                'mu_theta2'), #effect precip interranvual var
                     pg0 = TRUE)


# covariate effect on gen length -----------------------------------------------

#INTERPRETATION
#((e^param) - 1) * 100 = percent change in trait for every one unit change in covariate
#((e^(param * L)) - 1) * 100 = percent change in trait for every L unit change in covariate
mu_beta_ch <- MCMCvis::MCMCchains(fit, params = 'mu_beta')
mu_gamma1_ch <- MCMCvis::MCMCchains(fit, params = 'mu_gamma1')
mu_gamma2_ch <- MCMCvis::MCMCchains(fit, params = 'mu_gamma2')
mu_theta1_ch <- MCMCvis::MCMCchains(fit, params = 'mu_theta1')
mu_theta2_ch <- MCMCvis::MCMCchains(fit, params = 'mu_theta2')

# median((exp(mu_beta_ch * diff(range(DATA$lMass))) - 1) * 100)
# median((exp(mu_gamma1_ch * diff(range(DATA$temp_sd_season))) - 1) * 100)
# median((exp(mu_gamma2_ch * diff(range(DATA$temp_sd_year))) - 1) * 100)
# median((exp(mu_theta1_ch * diff(range(DATA$precip_cv_season))) - 1) * 100)
# median((exp(mu_theta2_ch * diff(range(DATA$precip_cv_year))) - 1) * 100)

#% change in gen length for 1 sd change in covariate
mu_beta_rs_ch <- (exp(mu_beta_ch * sd(bird_df4$lMass)) - 1) * 100
mu_gamma1_rs_ch <- (exp(mu_gamma1_ch * sd(bird_df4$temp_sd_season)) - 1) * 100
mu_gamma2_rs_ch <- (exp(mu_gamma2_ch * sd(bird_df4$temp_sd_year)) - 1) * 100
mu_theta1_rs_ch <- (exp(mu_theta1_ch * sd(bird_df4$precip_cv_season)) - 1) * 100
mu_theta2_rs_ch <- (exp(mu_theta2_ch * sd(bird_df4$precip_cv_year)) - 1) * 100


# cat plots ---------------------------------------------------------------

nn <- unique(bird_df4[,c('niche_idx', 'Trophic_niche')]) %>%
  dplyr::arrange(niche_idx)

#raw params
pdf(paste0(fig_dir, 'param-cat-raw-', run_date, '.pdf'),
    height = 14, width = 5)
MCMCvis::MCMCplot(fit,
                  params = c('mu_beta', 
                             'mu_gamma1', 'mu_gamma2',
                             'mu_theta1', 'mu_theta2'),
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  ref_ovl = TRUE,
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'Param estimates',
                  guide_lines = TRUE)
dev.off()

pdf(paste0(fig_dir, 'param-cat-rs-', run_date, '.pdf'),
    height = 14, width = 5)
MCMCvis::MCMCplot(cbind(#mu_beta_rs_ch, 
                        mu_gamma1_rs_ch,
                        mu_gamma2_rs_ch,
                        mu_theta1_rs_ch,
                        mu_theta2_rs_ch),
                  labels = c(#'log(Mass)',
                             'T season',
                             'T interannual var',
                             'P seasonality',
                             'P interannual var'),
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  ref_ovl = TRUE,
                  main = '% change in Gen Length for 1 sd change in cov',
                  guide_lines = TRUE)
dev.off()


pdf(paste0(fig_dir, 'gamma1-cat-', run_date, '.pdf'),
    height = 14, width = 5)
MCMCvis::MCMCplot(fit,
                  params = 'gamma1',
                  labels = nn$Trophic_niche,
                  main = 'Temp Seasonality',
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  ref_ovl = TRUE,
                  sz_thick = 3,
                  sz_thin = 3,
                  guide_lines = TRUE)
dev.off()

pdf(paste0(fig_dir, 'gamma2-cat-', run_date, '.pdf'),
    height = 14, width = 5)
MCMCvis::MCMCplot(fit,
                  params = 'gamma2',
                  labels = nn$Trophic_niche,
                  main = 'Temp Interannual',
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  ref_ovl = TRUE,
                  sz_thick = 3,
                  sz_thin = 3,
                  guide_lines = TRUE)
dev.off()

pdf(paste0(fig_dir, 'theta1-cat-', run_date, '.pdf'),
    height = 14, width = 5)
MCMCvis::MCMCplot(fit,
                  params = 'theta1',
                  labels = nn$Trophic_niche,
                  main = 'Precip Seasonality',
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  ref_ovl = TRUE,
                  sz_thick = 3,
                  sz_thin = 3,
                  guide_lines = TRUE)
dev.off()

pdf(paste0(fig_dir, 'theta2-cat-', run_date, '.pdf'),
    height = 14, width = 5)
MCMCvis::MCMCplot(fit,
                  params = 'theta2',
                  labels = nn$Trophic_niche,
                  main = 'Precip Interannual',
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  ref_ovl = TRUE,
                  sz_thick = 3,
                  sz_thin = 3,
                  guide_lines = TRUE)
dev.off()

pdf(paste0(fig_dir, 'alpha-cat-', run_date, '.pdf'),
    height = 14, width = 5)
MCMCvis::MCMCplot(fit,
                  params = 'alpha',
                  labels = nn$Trophic_niche,
                  main = 'Intercepts (mean Gen)',
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  ref_ovl = TRUE,
                  sz_thick = 3,
                  sz_thin = 3,
                  guide_lines = TRUE)
dev.off()

#by group
for (i in 1:NROW(niche_names))
{
  #i <- 10
  pp <- paste0(c('gamma1', 'gamma2', 'theta1', 'theta2'), 
             '[', i, ']')
  pdf(paste0(fig_dir, niche_names$Trophic.Niche[i], '-cat-', run_date, '.pdf'),
      height = 14, width = 5)
  MCMCvis::MCMCplot(fit,
                    params = pp,
                    exact = TRUE,
                    ISB = FALSE,
                    labels = c(#'log(Mass)',
                      'T season',
                      'T interannual var',
                      'P seasonality',
                      'P interannual var'),
                    main = niche_names$Trophic.Niche[i],
                    sz_labels = 1.5,
                    ci = c(89, 89),
                    ref_ovl = TRUE,
                    sz_thick = 3,
                    sz_thin = 3,
                    guide_lines = TRUE)
  dev.off()
}


pp <- c('gamma1', 'gamma2', 'theta1', 'theta2')
MCMCvis::MCMCplot(fit,
                  params = pp,
                  excl = 'mu',
                  exact = FALSE,
                  ISB = FALSE,
                  # labels = c(#'log(Mass)',
                  #   'T season',
                  #   'T interannual var',
                  #   'T spectra',
                  #   'P seasonality',
                  #   'P interannual var',
                  #   'P spectra'),
                  # main = niche_names$Trophic.Niche[i],
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  ref_ovl = TRUE,
                  sz_thick = 3,
                  sz_thin = 3,
                  guide_lines = TRUE)


# among group -------------------------------------------------------------

#seems to be fine within groups 
tdf <- dplyr::filter(bird_df4, Trophic_niche == nn$Trophic_niche[1])
tf1 <- lm(lGL ~ lMass + temp_sd_season + temp_sd_year +
            precip_cv_season + precip_cv_year, data = tdf)
car::vif(tf1)


# PPC ---------------------------------------------------------------------

#posterior predictive check
mu_ch <- MCMCvis::MCMCchains(fit, params = 'mu')
sigma_ch <- MCMCvis::MCMCchains(fit, params = 'sigma')
nu_ch <- MCMCvis::MCMCchains(fit, params = 'nu')

#500 iterations
sidx <- sample(1:NROW(mu_ch), size = 500)
y_rep <- matrix(NA, nrow = length(sidx), ncol = NCOL(mu_ch))
for (i in 1:length(sidx))
{
  #i <- 1
  print(paste0('iter: ', i, ' of ', length(sidx)))
  for (j in 1:NCOL(mu_ch))
  {
    #t-dis
    eps <- rt(n = 1, df = nu_ch[sidx[i],1]) * sigma_ch[sidx[i],1]
    y_rep[i,j] <- mu_ch[sidx[i],j] + eps
    # y_rep[i,j] <- rnorm(1, mu_ch[sidx[i],j], sigma_ch[sidx[i],1])
  }
}

pdf(paste0(fig_dir, 'varniche_PPC-', run_date, '.pdf'), height = 5, width = 5)
plot(density(DATA$y), col = 'black', lwd = 3, xlim = c(0, 3.5), ylim = c(0, 1.5))
for (i in 1:500)
{
  lines(density(y_rep[i,]), col = rgb(1,0,0,0.2))
}
dev.off()


# R^2 ---------------------------------------------------------------------

#Gelman et al. 2019 - American Statistician
#var predicted / (var predicted + var residuals)
#in other words:
#explained var / (explained var + resid var)

#with mass - 0.74
var_pred <- apply(mu_ch, 1, var)
var_resid <- apply(sweep(mu_ch, 2, DATA$y), 1, var)
r2_ch <- var_pred / (var_pred + var_resid)
hist(r2_ch)

#no mass - by group
alpha_mn <- MCMCvis::MCMCpstr(fit, params = 'alpha')[[1]]
beta_mn <- MCMCvis::MCMCpstr(fit, params = 'beta')[[1]]
gamma1_mn <- MCMCvis::MCMCpstr(fit, params = 'gamma1')[[1]]
gamma2_mn <- MCMCvis::MCMCpstr(fit, params = 'gamma2')[[1]]
theta1_mn <- MCMCvis::MCMCpstr(fit, params = 'theta1')[[1]]
theta2_mn <- MCMCvis::MCMCpstr(fit, params = 'theta2')[[1]]

#500 iterations
pmu <- rep(NA, length(DATA$y))
for (i in 1:length(DATA$y))
{
  #i <- 1
  pmu[i] <- alpha_mn[DATA$f_id[i]] + 
    #beta_mn[DATA$f_id[i]] * DATA$lMass +
    gamma1_mn[DATA$f_id[i]] * DATA$temp_sd_season[i] +
    gamma2_mn[DATA$f_id[i]] * DATA$temp_sd_year[i] +
    theta1_mn[DATA$f_id[i]] * DATA$precip_cv_season[i] +
    theta2_mn[DATA$f_id[i]] * DATA$precip_cv_year[i]
}

ug <- unique(DATA$f_id)
r2_vec <- rep(NA, length(ug))
for (i in 1:length(ug))
{
  #i <- 1
  tidx <- which(DATA$f_id == ug[i])
  tres <- DATA$y[tidx] - pmu[tidx]
  tpmu <- pmu[tidx]
  r2_vec[i] <- var(tpmu) / (var(tpmu) + var(tres))
}

dplyr::left_join(niche_names, data.frame(f_id = ug,
                                         r2 = round(r2_vec, 3)),
                 by = 'f_id')


# VIF ---------------------------------------------------------------------

#covariates as a function of other covariates
tf1 <- lm(lMass ~ temp_sd_year + temp_sd_season +
            precip_cv_year + precip_cv_season,data = bird_df4)
stf1 <- summary(tf1)

tf2 <- lm(temp_sd_year ~ lMass + temp_sd_season + 
            precip_cv_year + precip_cv_season, data = bird_df4)
stf2 <- summary(tf2)

tf3 <- lm(temp_sd_season ~ lMass + temp_sd_year + 
            precip_cv_year + precip_cv_season, data = bird_df4)
stf3 <- summary(tf3)

tf4 <- lm(precip_cv_year ~ lMass + temp_sd_year + temp_sd_season + 
            precip_cv_season, data = bird_df)
stf4 <- summary(tf4)

tf5 <- lm(precip_cv_season ~ lMass + temp_sd_year + temp_sd_season + 
            precip_cv_year, data = bird_df)
stf5 <- summary(tf5)

#calc VIF per covariate
1 / (1 - stf1$r.squared) #lMass
1 / (1 - stf2$r.squared) #temp_sd_year
1 / (1 - stf3$r.squared) #temp_sd_season
1 / (1 - stf3$r.squared) #precip_cv_year
1 / (1 - stf5$r.squared) #precip_cv_season

#correlation
cor(as.matrix(dplyr::select(bird_df,
                            lMass, temp_sd_year, temp_sd_season,
                            precip_cv_year, precip_cv_season)))
