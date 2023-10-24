####################
# Fit Bayes model - gl ~ env + phylo + vint (gamma)
####################


# specify dir -------------------------------------------------------------

# dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# sc_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
dir <- '/mnt/research/ibeem/variability/'
sc_dir <- '/mnt/home/ccy/variability/'
run_date <- '2023-10-17'


# load packages -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(MCMCvis)
library(ape)
library(geiger)
library(phytools)


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
                # lGL = log(GenLength),
                lGL = GenLength,
                # lGL = 1/GenLength,
                lAb = log(Modeled_age_first_breeding),
                lMl = log(Modeled_max_longevity),
                species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  #drop duplicated species (for now)
  dplyr::group_by(species) %>%
  dplyr::slice_head() %>%
  dplyr::ungroup()

#subset out just traits of interest
bird_df2 <- dplyr::select(bird_df, 
                          species,
                          lMass,
                          lGL,
                          Modeled_survival,
                          lAb,
                          lMl,
                          Trophic_niche = Trophic.Niche,
                          temp_mean,
                          temp_sd_year,
                          temp_sd_season,
                          precip_cv_year,
                          precip_cv_season)


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
                    lGL) %>%
  as.matrix()
row.names(dd) <- bird_df4$species

#get estimate of Pagel's kappa to scale phylogeny
fit_ka <- geiger::fitContinuous(pr_tree2, dd[,'lGL'], model = "kappa")

#rescale tree
pr_tree_k <- phytools::rescale(pr_tree, 'kappa', 
                               kappa = fit_ka$opt$kappa, sigsq = fit_ka$opt$sigsq)

#get corr matrix of rescaled tree
Rho <- ape::vcv.phylo(pr_tree_k, corr = TRUE)


# niche levels ------------------------------------------------------------

#assign missing species (owls) to Vertivore
bird_df4$Trophic_niche[which(is.na(bird_df4$Trophic_niche))] <- "Vertivore"

bird_df4$niche_idx <- as.numeric(factor(bird_df4$Trophic_niche))
niche_names <- levels(factor(bird_df4$Trophic_niche))


# scale/prep data ---------------------------------------------------------

#scalars for data - smaller number for larger param value (opposite for y)
lMass_scalar <- 1
temp_sd_season_scalar <- 0.05
temp_sd_year_scalar <- 0.1
precip_cv_season_scalar <- 0.2
precip_cv_year_scalar <- 0.5
y_scalar <- 1

#center predictors
tt <- data.frame(lMass = bird_df4$lMass * 
                   lMass_scalar,
                 temp_sd_season = bird_df4$temp_sd_season * 
                   temp_sd_season_scalar,
                 temp_sd_year = bird_df4$temp_sd_year * 
                   temp_sd_year_scalar,
                 precip_cv_season = bird_df4$precip_cv_season * 
                   precip_cv_season_scalar,
                 precip_cv_year = bird_df4$precip_cv_year * 
                   precip_cv_year_scalar) %>%
  apply(2, function(x) scale(x, scale = FALSE)[,1])


# fit model ---------------------------------------------------------------

DATA <- list(N = NROW(bird_df4),
             Y = bird_df4$lGL * y_scalar,
             K = NCOL(tt),
             J = length(unique(bird_df4$niche_idx)),
             X = tt,
             niche_idx = bird_df4$niche_idx,
             mu_kappa = 0,
             sigma_kappa = 1,
             Rho = Rho) #corr matrix

# summary(lm(DATA$Y ~ DATA$X[,1] +
#              DATA$X[,2] +
#              DATA$X[,3] +
#              DATA$X[,4] +
#              DATA$X[,5]))

options(mc.cores = parallel::detectCores())

# DELTA <- 0.92
# TREE_DEPTH <- 10
# STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 2000

#compile model
mod <- cmdstanr::cmdstan_model(paste0(sc_dir, 'Scripts/Model_files/5-phylo-vint-gamma.stan'))

#sample
fit <- mod$sample(
  data = DATA,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 500)
# max_treedepth = TREE_DEPTH
# adapt_delta = DELTA
# step_size = STEP_SIZE


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fit, 
                  round = 4,
                  file_name = paste0('bird-gl-phylo-vint-gamma-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('bird-gl-phylo-vint-gamma-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('bird-gl-phylo-vint-gamma-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('bird-gl-phylo-vint-gamma-data-', run_date),
                  cp_file = c(paste0(sc_dir, 'Scripts/Model_files/5-phylo-vint-gamma.stan'), 
                              paste0(sc_dir, 'Scripts/5-model/5-gl-phylo-vint-gamma.R')),
                  cp_file_names = c(paste0('5-phylo-vint-gamma-', run_date, '.stan'),
                                    paste0('5-gl-phylo-vint-gamma-', run_date, '.R')))

fig_dir <- paste0(dir, 'Results/bird-gl-phylo-vint-gamma-', run_date, '/')

# fit <- readRDS(paste0(dir, '/Results/bird-gl-phylo-vint-gamma-', run_date,
#                       '/bird-gl-phylo-vint-gamma-fit-', run_date, '.rds'))
# library(shinystan)
# shinystan::launch_shinystan(fit)


# residuals ---------------------------------------------------------------
# 
# alpha_mn <- MCMCvis::MCMCpstr(fit, params = 'alpha')[[1]]
# kappa_mn <- MCMCvis::MCMCpstr(fit, params = 'kappa')[[1]]
# gamma_mn <- MCMCvis::MCMCpstr(fit, params = 'gamma')[[1]]
# beta_mn <- MCMCvis::MCMCpstr(fit, params = 'beta')[[1]]
# 
# mu_mn <- kappa_mn + gamma_mn[DATA$niche_idx] + alpha_mn + (DATA$X %*% beta_mn)[,1]
# resid <- DATA$Y - mu_mn


# Summary -----------------------------------------------------------------

#model summary
MCMCvis::MCMCsummary(fit, round = 3, 
                     params = c('beta', 'kappa', 'gamma', 'a', 'sigma_phylo'),
                     pg0 = TRUE)


# covariate effect on LH trait -----------------------------------------------

#INTERPRETATION
#((e^param) - 1) * 100 = percent change in trait for every one unit change in covariate
#((e^(param * L)) - 1) * 100 = percent change in trait for every L unit change in covariate
beta1_ch <- MCMCvis::MCMCchains(fit, params = 'beta[1]', 
                                exact = TRUE, ISB = FALSE) * 
  lMass_scalar * y_scalar
beta2_ch <- MCMCvis::MCMCchains(fit, params = 'beta[2]', 
                                exact = TRUE, ISB = FALSE) * 
  temp_sd_season_scalar * y_scalar
beta3_ch <- MCMCvis::MCMCchains(fit, params = 'beta[3]', 
                                exact = TRUE, ISB = FALSE) * 
  temp_sd_year_scalar * y_scalar
beta4_ch <- MCMCvis::MCMCchains(fit, params = 'beta[4]', 
                                exact = TRUE, ISB = FALSE) * 
  precip_cv_season_scalar * y_scalar
beta5_ch <- MCMCvis::MCMCchains(fit, params = 'beta[5]', 
                                exact = TRUE, ISB = FALSE) * 
  precip_cv_year_scalar * y_scalar

# median((exp(beta_ch * diff(range(DATA$lMass))) - 1) * 100)
# median((exp(gamma1_ch * diff(range(DATA$temp_sd_season))) - 1) * 100)
# median((exp(gamma2_ch * diff(range(DATA$temp_sd_year))) - 1) * 100)
# median((exp(theta1_ch * diff(range(DATA$precip_cv_season))) - 1) * 100)
# median((exp(theta2_ch * diff(range(DATA$precip_cv_year))) - 1) * 100)

#scaling cov to measured scale bc transformed param est
#% change in LH trait for 1 sd change in covariate
beta1_rs_ch <- (exp(beta1_ch * sd(tt[,1] / lMass_scalar)) - 1) * 100
beta2_rs_ch <- (exp(beta2_ch * sd(tt[,2] / temp_sd_season_scalar)) - 1) * 100
beta3_rs_ch <- (exp(beta3_ch * sd(tt[,3] / temp_sd_year_scalar)) - 1) * 100
beta4_rs_ch <- (exp(beta4_ch * sd(tt[,4] / precip_cv_season_scalar)) - 1) * 100
beta5_rs_ch <- (exp(beta5_ch * sd(tt[,5] / precip_cv_year_scalar)) - 1) * 100


# added variable and partial resid plots ------------------------------------------------

# https://www.wikiwand.com/en/Partial_residual_plot
pr_fun <- function(num, nm)
{
  names <- c('Mass', 'Temp seasonality', 'Temp interannual',
             'Precip seasonality', 'Precip interannual')
  
  #partial residuals
  pr <- resid + (beta_mn[num] * DATA$X[,num])
  
  pdf(paste0(fig_dir, nm, '-pr-', run_date, '.pdf'),
      height = 5, width = 5)
  plot(DATA$X[,num], pr, col = rgb(0,0,0,0.2), pch = 19,
       xlab = 'Predictor',
       ylab = 'Partial residual',
       main = names[num])
  abline(h = 0, col = 'grey', lwd = 4, lty = 2)
  abline(a = 0, b = beta_mn[num], col = rgb(1,0,0,0.5), lwd = 4)
  dev.off()
}

pr_fun(num = 1, nm = 'mass') #Mass
pr_fun(num = 2, nm = 'temp-season') #temp season
pr_fun(num = 3, nm = 'temp-year') #temp year
pr_fun(num = 4, nm = 'precip-season') #precip season
pr_fun(num = 5, nm = 'precip-year') #precip year


# cat plots ---------------------------------------------------------------

pdf(paste0(fig_dir, 'param-cat-raw-', run_date, '.pdf'),
    height = 5, width = 5)
MCMCvis::MCMCplot(fit,
                  params = 'beta',
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  ref_ovl = TRUE,
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'Param estimates',
                  guide_lines = TRUE)
dev.off()

pdf(paste0(fig_dir, 'param-cat-rs-', run_date, '.pdf'),
    height = 5, width = 5)
MCMCvis::MCMCplot(cbind(beta2_rs_ch,
                        beta3_rs_ch,
                        beta4_rs_ch,
                        beta5_rs_ch),
                  labels = c('T seasonality',
                             'T interannual var',
                             'P seasonality',
                             'P interannual var'),
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = '% change GL for 1 sd change in cov',
                  guide_lines = TRUE)
dev.off()


pdf(paste0(fig_dir, 'gamma-cat-', run_date, '.pdf'),
    height = 5, width = 5)
MCMCvis::MCMCplot(fit, 
                  params = 'gamma',
                  labels = niche_names,
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = 'Niche group intercept',
                  guide_lines = TRUE)
dev.off()


# PPC ---------------------------------------------------------------------

# # PPC - normal model
# sigma_ch <- MCMCvis::MCMCchains(fit, params = 'sigma')
# alpha_ch <- MCMCvis::MCMCchains(fit, params = 'alpha')
# kappa_ch <- MCMCvis::MCMCchains(fit, params = 'kappa')
# gamma_ch <- MCMCvis::MCMCchains(fit, params = 'gamma')
# beta_ch <- MCMCvis::MCMCchains(fit, params = 'beta')
# 
# #500 iterations
# sidx <- sample(1:NROW(sigma_ch), size = 500)
# mu_rep <- matrix(NA, nrow = length(sidx), ncol = length(DATA$Y))
# y_rep <- matrix(NA, nrow = length(sidx), ncol = length(DATA$Y))
# for (i in 1:length(sidx))
# {
#   #i <- 1
#   print(paste0('iter: ', i, ' of ', length(sidx)))
#   for (j in 1:length(DATA$Y))
#   {
#     #j <- 1
#     #t-dis
#     # eps <- rt(n = 1, df = nu_ch[sidx[i],1]) * sigma_ch[sidx[i],1]
#     # y_rep[i,j] <- mu_ch[sidx[i],j] + eps
#     
#     mu_rep[i, j] <- kappa_ch[sidx[i], 1] + gamma_ch[sidx[i], DATA$niche_idx[j]] + 
#       alpha_ch[sidx[i], 1] + (DATA$X[j,] %*% beta_ch[sidx[i],])[,1]
#     
#     y_rep[i,j] <- rnorm(1, mu_rep[i, j], sigma_ch[sidx[i], 1])
#   }
# }
# 
# pdf(paste0(fig_dir, 'PPC-', run_date, '.pdf'), height = 5, width = 5)
# plot(density(DATA$Y), col = 'black', lwd = 3, ylim = c(0, 1))#, xlim = c(0, 3.5))
# for (i in 1:500)
# {
#   lines(density(y_rep[i,]), col = rgb(1,0,0,0.05))
# }
# dev.off()


# R^2 ---------------------------------------------------------------------

#Gelman et al. 2019 - American Statistician
#var predicted / (var predicted + var residuals)
#in other words:
#explained var / (explained var + resid var)

# var_pred <- apply(mu_rep, 1, var)
# var_resid <- apply(sweep(mu_rep, 2, DATA$Y), 1, var)
# r2_ch <- var_pred / (var_pred + var_resid)
# hist(r2_ch)


# VIF ---------------------------------------------------------------------

#covariates as a function of other covariates
# tf1 <- lm(lMass ~ temp_sd_year + temp_sd_season +
#             precip_cv_year + precip_cv_season, 
#           data = bird_df)
# stf1 <- summary(tf1)
# 
# tf2 <- lm(temp_sd_year ~ lMass + temp_sd_season +
#             precip_cv_year + precip_cv_season, data = 
#             bird_df)
# stf2 <- summary(tf2)
# 
# tf3 <- lm(temp_sd_season ~ lMass + temp_sd_year +
#             precip_cv_year + precip_cv_season, 
#           data = bird_df)
# stf3 <- summary(tf3)
# 
# tf4 <- lm(precip_cv_year ~ lMass + temp_sd_year + temp_sd_season + 
#             precip_cv_season, 
#           data = bird_df)
# stf4 <- summary(tf4)
# 
# tf5 <- lm(precip_cv_season ~ lMass + temp_sd_year + temp_sd_season + 
#             precip_cv_year, 
#           data = bird_df)
# stf5 <- summary(tf5)
# 
# 
# #calc VIF per covariate
# 1 / (1 - stf1$r.squared) #lMass
# 1 / (1 - stf2$r.squared) #temp_sd_year
# 1 / (1 - stf3$r.squared) #temp_sd_season
# 1 / (1 - stf4$r.squared) #precip_cv_year
# 1 / (1 - stf5$r.squared) #precip_cv_season
# 
# #correlation
# cor(as.matrix(dplyr::select(bird_df3,
#                             lMass, temp_sd_year, temp_sd_season,
#                             precip_cv_year, precip_cv_season)))


