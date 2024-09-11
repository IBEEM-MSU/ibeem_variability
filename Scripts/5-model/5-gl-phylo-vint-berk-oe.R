# TITLE:            Fit Bayesian model 
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Kelly Kapsar, Adriana Uscanga, Peter J. Williams, Jeffrey W. Doser, Lala Kounta, Phoebe L. Zarnetske
# DATA INPUT:       Bird consensus tree (4b), 
# DATA OUTPUT:      Model output summary and figures (caterpillar plots, partial residual plots, etc.)
# DATE:             July 2024
# OVERVIEW:         Fits Bayesian model to data (gl ~ env + phylo + vint + berk + oe)


rm(list = ls())


# load environment variables ------------------------------------------------

source("./Scripts/0-config.R")


# set directories and date ------------------------------------------------

sc_dir <- "./" # Assumes you're running the code from dir sourced in 0-config.Rsc_dir <- dir
run_date <- Sys.Date()


# load packages -----------------------------------------------------------

library(tidyverse)
library(cmdstanr)
# set_cmdstan_path("/mnt/home/kapsarke/.cmdstan/cmdstan-2.34.1")
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

'%ni%' <- Negate('%in%')
bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data-birdtree2.csv')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::filter(Order %ni% or_excl,
                Migration == 1,
                !is.na(temp_sd_year)) %>%
  dplyr::mutate(lMass = log(Mass),
                lGL = log(GenLength),
                species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  dplyr::rename(Trophic_niche = Trophic.Niche) %>%
  #how much var will change in 1 generation (in sds)
  #((degrees / year) * (year / generation) * (1 / degrees) = degrees (sd) / generation
  dplyr::mutate(temp_delta = (temp_slope * GenLength) / temp_sd_year,
                precip_delta = (precip_slope * GenLength) / precip_sd_year)

#subset out just traits of interest
bird_df2 <- dplyr::select(bird_df,
                          ID,
                          Birdtree_name, 
                          Avonet_name, 
                          species, #formatting in line with Birdtree phylo
                          Family,
                          Order,
                          lGL,
                          lMass,
                          Trophic_niche,
                          temp_sd_year,
                          temp_sd_year_sd_space,
                          temp_sd_season,
                          temp_sd_season_sd_space,
                          precip_cv_year,
                          precip_cv_year_sd_space,
                          precip_cv_season,
                          precip_cv_season_sd_space,
                          temp_delta,
                          precip_delta)

# write.csv(bird_df2, paste0(dir, "/data/L3/final-bird-data-for-archival.csv"),
#           row.names = FALSE)


# phylo -------------------------------------------------------------------

#load consensus tree - bird.phylo
load(paste0(dir, 'data/L3/bird-consensus-tree.rda'))

#prune tree
#species not found in both datasets (species to drop from tree)
nm <- setdiff(bird.phylo$tip.label, bird_df2$species)

#prune specified tips from tree
pr_tree <- ape::drop.tip(bird.phylo, nm)

#get idx
j_idx3 <- dplyr::left_join(data.frame(species = pr_tree$tip.label),
                           data.frame(idx = 1:NROW(bird_df2), bird_df2),
                           by = 'species')

#apply
bird_df3 <- bird_df2[j_idx3$idx,]

#make tree binary
pr_tree2 <- ape::multi2di(pr_tree)

#make response into matrix with species as rownames
dd <- dplyr::select(bird_df3,
                    lGL) %>%
  as.matrix()
row.names(dd) <- bird_df3$species

#get estimate of Pagel's lambda to scale phylogeny
#https://lukejharmon.github.io/pcm/chapter6_beyondbm/
fit_la <- geiger::fitContinuous(pr_tree2, dd[,'lGL'], model = "lambda")

#rescale tree
pr_tree_k <- phytools::rescale(pr_tree, 'lambda',
                               lambda = fit_la$opt$lambda, sigsq = fit_la$opt$sigsq)

#get corr matrix of rescaled tree
Rho <- ape::vcv.phylo(pr_tree_k, corr = TRUE)


# niche levels ------------------------------------------------------------

#assign missing species (owls) to Vertivore
bird_df3$Trophic_niche[which(is.na(bird_df3$Trophic_niche))] <- "Vertivore"

bird_df3$niche_idx <- as.numeric(factor(bird_df3$Trophic_niche))
niche_names <- levels(factor(bird_df3$Trophic_niche))


# summary metrics ---------------------------------------------------------

#7477 species
NROW(bird_df3)
#29 orders
length(unique(bird_df3$Order))
#198 families
length(unique(bird_df3$Family))

#min and max gen lengths
min_idx <- which.min(bird_df3$lGL)
max_idx <- which.max(bird_df3$lGL)
bird_df3[c(min_idx, max_idx),] %>%
  dplyr::mutate(GL = exp(lGL)) %>%
  dplyr::select(species, GL)

#IQR
qs <- quantile(exp(bird_df3$lGL), probs = c(0.25, 0.5, 0.75))
iqr <- qs[3] - qs[1]


# scale/prep data ---------------------------------------------------------

#scalars for data - smaller number for larger param value (opposite for y)
lMass_scalar <- 1
temp_sd_season_scalar <- 0.2
temp_sd_year_scalar <- 0.2
precip_cv_season_scalar <- 0.1
precip_cv_year_scalar <- 0.5
y_scalar <- 2

# #center predictors
# tt <- data.frame(lMass = bird_df3$lMass *
#                    lMass_scalar,
#                  temp_sd_season = bird_df3$temp_sd_season *
#                    temp_sd_season_scalar,
#                  temp_sd_year = bird_df3$temp_sd_year *
#                    temp_sd_year_scalar,
#                  precip_cv_season = bird_df3$precip_cv_season *
#                    precip_cv_season_scalar,
#                  precip_cv_year = bird_df3$precip_cv_year *
#                    precip_cv_year_scalar) %>%
#   apply(2, function(x) scale(x, scale = FALSE)[,1])


# delta -------------------------------------------------------------------

median(bird_df3$temp_delta)
t_delta_qs <- quantile(bird_df3$temp_delta, probs = c(0.25, 0.5, 0.75))
t_delta_iqr <- t_delta_qs[3] - t_delta_qs[1]
sum(abs(bird_df3$temp_delta) > 0.1) / NROW(bird_df3)
sum(abs(bird_df3$temp_delta) > 0.3) / NROW(bird_df3)

median(bird_df3$precip_delta)
p_delta_qs <- quantile(bird_df3$precip_delta, probs = c(0.25, 0.5, 0.75))
p_delta_iqr <- p_delta_qs[3] - p_delta_qs[1]
sum(abs(bird_df3$precip_delta) > 0.1) / NROW(bird_df3)
sum(abs(bird_df3$precip_delta) > 0.3) / NROW(bird_df3)

#correlation between delta temp and delta precip
cor(bird_df3$temp_delta, bird_df3$precip_delta)

#correlation between delta and gen length
cor(bird_df$temp_delta, bird_df$GenLength)
cor(bird_df$precip_delta, bird_df$GenLength)

#correlation between delta and relative change
rc_temp <- bird_df$temp_slope / bird_df$temp_sd_year
cor(bird_df$temp_delta, rc_temp)
rc_precip <- bird_df$precip_slope / bird_df$precip_sd_year
cor(bird_df$precip_delta, rc_precip)


# uncertainty GL ----------------------------------------------------------

#r2 = 1 - (resid var / total var)
#resid var / total var = (1 - r2)
#(1 - r2) * total var = resid var
#r2 from 1d-uncertainty-gl.R - derived from Bird et al. uncertainty
sd_Y <- sqrt((1 - 0.87) * var(bird_df$lGL))


# ensure no zeros for cov sd ----------------------------------------------

zidx_t_seas <- which(bird_df3$temp_sd_season_sd_space == 0)
bird_df3$temp_sd_season_sd_space[zidx_t_seas] <- 0.001
zidx_t_year <- which(bird_df3$temp_sd_year_sd_space == 0)
bird_df3$temp_sd_year_sd_space[zidx_t_year] <- 0.001
zidx_p_seas <- which(bird_df3$precip_cv_season_sd_space == 0)
bird_df3$precip_cv_season_sd_space[zidx_p_seas] <- 0.001
zidx_p_year <- which(bird_df3$precip_cv_year_sd_space == 0)
bird_df3$precip_cv_year_sd_space[zidx_p_year] <- 0.001


# fit model ---------------------------------------------------------------

DATA <- list(N = NROW(bird_df3),
             Y_obs = bird_df3$lGL * y_scalar,
             sd_Y = sd_Y * y_scalar,
             # K = NCOL(tt),
             K = 5,
             J = length(unique(bird_df3$niche_idx)),
             lMass = bird_df3$lMass * lMass_scalar,
             temp_sd_season = bird_df3$temp_sd_season *
               temp_sd_season_scalar,
             temp_sd_year = bird_df3$temp_sd_year *
               temp_sd_year_scalar,
             precip_cv_season = bird_df3$precip_cv_season *
               precip_cv_season_scalar,
             precip_cv_year = bird_df3$precip_cv_year *
               precip_cv_year_scalar,
             temp_sd_season_sd_space = bird_df3$temp_sd_season_sd_space *
               temp_sd_season_scalar,
             temp_sd_year_sd_space = bird_df3$temp_sd_year_sd_space *
               temp_sd_year_scalar,
             precip_cv_season_sd_space = bird_df3$precip_cv_season_sd_space *
               precip_cv_season_scalar,
             precip_cv_year_sd_space = bird_df3$precip_cv_year_sd_space *
               precip_cv_year_scalar,
             niche_idx = bird_df3$niche_idx,
             mu_kappa = 1,
             sigma_kappa = 1,
             Rho = Rho, #corr matrix
             pro_data = bird_df3)

# summary(lm(DATA$Y ~ DATA$X[,1] +
#              DATA$X[,2] +
#              DATA$X[,3] +
#              DATA$X[,4] +
#              DATA$X[,5]))

options(mc.cores = parallel::detectCores())

# DELTA <- 0.92
# TREE_DEPTH <- 10
STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 4000

#compile model
mod <- cmdstanr::cmdstan_model(paste0(sc_dir, 
                                      'Scripts/Model_files/5-phylo-vint-berk-oe.stan'),
                               force_recompile = TRUE)

print("STARTING TO FIT MODEL")

#sample
fit <- mod$sample(
  data = DATA,
  chains = CHAINS,
  iter_sampling = ITER / 2,
  iter_warmup = ITER / 2,
  parallel_chains = CHAINS,
  refresh = 500,
  step_size = STEP_SIZE)
# max_treedepth = TREE_DEPTH
# adapt_delta = DELTA


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fit,
                  round = 4,
                  file_name = paste0('bird-gl-phylo-vint-berk-oe-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('bird-gl-phylo-vint-berk-oe-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('bird-gl-phylo-vint-berk-oe-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('bird-gl-phylo-vint-berk-oe-data-', run_date),
                  cp_file = c(paste0(sc_dir, 'Scripts/Model_files/5-phylo-vint-berk-oe.stan'),
                              paste0(sc_dir, 'Scripts/5-model/5-gl-phylo-vint-berk-oe.R')),
                  cp_file_names = c(paste0('5-phylo-vint-berk-oe-', run_date, '.stan'),
                                    paste0('5-gl-phylo-vint-berk-oe-', run_date, '.R')))

fig_dir <- paste0(dir, 'Results/bird-gl-phylo-vint-berk-oe-', run_date, '/')

# fit <- readRDS(paste0(dir, '/Results/bird-gl-phylo-vint-berk-oe-', run_date,
#                       '/bird-gl-phylo-vint-berk-oe-fit-', run_date, '.rds'))
# library(shinystan)
# shinystan::launch_shinystan(fit)


# residuals ---------------------------------------------------------------

kappa_mn <- MCMCvis::MCMCpstr(fit, params = 'kappa')[[1]]
gamma_mn <- MCMCvis::MCMCpstr(fit, params = 'gamma')[[1]]
alpha_mn <- MCMCvis::MCMCpstr(fit, params = 'alpha')[[1]]
beta_mn <- MCMCvis::MCMCpstr(fit, params = 'beta')[[1]]
yt_mn <- MCMCvis::MCMCpstr(fit, params = 'Y_t')[[1]]

#extract latent cov
tr_temp_sd_season_mn <- MCMCvis::MCMCpstr(fit, params = 'tr_temp_sd_season',
                                exact = FALSE, ISB = TRUE)[[1]]
tr_temp_sd_year_mn <- MCMCvis::MCMCpstr(fit, params = 'tr_temp_sd_year',
                                          exact = FALSE, ISB = TRUE)[[1]]
tr_precip_cv_season_mn <- MCMCvis::MCMCpstr(fit, params = 'tr_precip_cv_season',
                                          exact = FALSE, ISB = TRUE)[[1]]
tr_precip_cv_year_mn <- MCMCvis::MCMCpstr(fit, params = 'tr_precip_cv_year',
                                          exact = FALSE, ISB = TRUE)[[1]]

#get linear predictor
mu_mn <- kappa_mn + 
  gamma_mn[DATA$niche_idx] + 
  alpha_mn + 
  beta_mn[1] * DATA$lMass +
  beta_mn[2] * tr_temp_sd_season_mn +
  beta_mn[3] * tr_temp_sd_year_mn +
  beta_mn[4] * tr_precip_cv_season_mn +
  beta_mn[5] * tr_precip_cv_year_mn 

#resids
resid <- DATA$Y_obs - mu_mn


# residuals ch ---------------------------------------------------------------

kappa_ch <- MCMCvis::MCMCchains(fit, params = 'kappa')
gamma_ch <- MCMCvis::MCMCchains(fit, params = 'gamma')
alpha_ch <- MCMCvis::MCMCchains(fit, params = 'alpha')
beta_ch <- MCMCvis::MCMCchains(fit, params = 'beta')
yt_ch <- MCMCvis::MCMCchains(fit, params = 'Y_t')

#extract latent cov
tr_temp_sd_season_ch <- MCMCvis::MCMCchains(fit, params = 'tr_temp_sd_season',
                                          exact = FALSE, ISB = TRUE)
tr_temp_sd_year_ch <- MCMCvis::MCMCchains(fit, params = 'tr_temp_sd_year',
                                        exact = FALSE, ISB = TRUE)
tr_precip_cv_season_ch <- MCMCvis::MCMCchains(fit, params = 'tr_precip_cv_season',
                                            exact = FALSE, ISB = TRUE)
tr_precip_cv_year_ch <- MCMCvis::MCMCchains(fit, params = 'tr_precip_cv_year',
                                          exact = FALSE, ISB = TRUE)

#get linear predictor
mu_ch <- matrix(NA, nrow = NROW(kappa_ch), ncol = length(mu_mn))
for (i in 1:NROW(kappa_ch))
{
  #i <- 1
  print(paste0('processing ', i, ' of ', NROW(kappa_ch)))
  mu_ch[i,] <- kappa_ch[i,] + 
    gamma_ch[i,DATA$niche_idx] + 
    alpha_ch[i,] + 
    beta_ch[i,1] * DATA$lMass +
    beta_ch[i,2] * tr_temp_sd_season_ch[i,] +
    beta_ch[i,3] * tr_temp_sd_year_ch[i,] +
    beta_ch[i,4] * tr_precip_cv_season_ch[i,] +
    beta_ch[i,5] * tr_precip_cv_year_ch[i,]
}

#resids - sweep out y from linear predictor and multiple by 1 (equivalent to y - mu)
resid_ch <- sweep(mu_ch, 2, DATA$Y_obs) * -1
resid_mn <- apply(resid_ch, 2, mean)
resid_sd <- apply(resid_ch, 2, sd)


# Summary -----------------------------------------------------------------

#model summary
# MCMCvis::MCMCsummary(fit, round = 3,
#                      params = c('beta', 'kappa', 'gamma', 'sigma', 'sigma_phylo'),
#                      pg0 = TRUE)


# covariate effect on LH trait -----------------------------------------------

#INTERPRETATION
#((e^param) - 1) * 100 = percent change in trait for every one unit change in covariate
#((e^(param * L)) - 1) * 100 = percent change in trait for every L unit change in covariate
beta1_ch <- MCMCvis::MCMCchains(fit, params = 'beta[1]',
                                exact = TRUE, ISB = FALSE) *
  lMass_scalar / y_scalar
beta2_ch <- MCMCvis::MCMCchains(fit, params = 'beta[2]',
                                exact = TRUE, ISB = FALSE) *
  temp_sd_season_scalar / y_scalar
beta3_ch <- MCMCvis::MCMCchains(fit, params = 'beta[3]',
                                exact = TRUE, ISB = FALSE) *
  temp_sd_year_scalar / y_scalar
beta4_ch <- MCMCvis::MCMCchains(fit, params = 'beta[4]',
                                exact = TRUE, ISB = FALSE) *
  precip_cv_season_scalar / y_scalar
beta5_ch <- MCMCvis::MCMCchains(fit, params = 'beta[5]',
                                exact = TRUE, ISB = FALSE) *
  precip_cv_year_scalar / y_scalar

#scaling cov to measured scale bc transformed param est
#% change in LH trait for 1 sd change in covariate
beta1_rs_ch <- (exp(beta1_ch * sd(DATA$lMass / lMass_scalar)) - 1) * 100
beta2_rs_ch <- (exp(beta2_ch * sd(tr_temp_sd_season_mn / temp_sd_season_scalar)) - 1) * 100
beta3_rs_ch <- (exp(beta3_ch * sd(tr_temp_sd_year_mn / temp_sd_year_scalar)) - 1) * 100
beta4_rs_ch <- (exp(beta4_ch * sd(tr_precip_cv_season_mn / precip_cv_season_scalar)) - 1) * 100
beta5_rs_ch <- (exp(beta5_ch * sd(tr_precip_cv_year_mn / precip_cv_year_scalar)) - 1) * 100

mean(beta1_rs_ch)
quantile(beta1_rs_ch, probs = c(0.055, 0.5, 0.945))
sum(beta1_rs_ch < 0) / NROW(beta1_rs_ch)

mean(beta2_rs_ch)
quantile(beta2_rs_ch, probs = c(0.055, 0.5, 0.945))
sum(beta2_rs_ch < 0) / NROW(beta2_rs_ch)

mean(beta3_rs_ch)
quantile(beta3_rs_ch, probs = c(0.055, 0.5, 0.945))
sum(beta3_rs_ch > 0) / NROW(beta3_rs_ch)

mean(beta4_rs_ch)
quantile(beta4_rs_ch, probs = c(0.055, 0.5, 0.945))
sum(beta4_rs_ch < 0) / NROW(beta4_rs_ch)

mean(beta5_rs_ch)
quantile(beta5_rs_ch, probs = c(0.055, 0.5, 0.945))
sum(beta5_rs_ch > 0) / NROW(beta5_rs_ch)

median((exp(beta1_ch * diff(range(DATA$lMass / lMass_scalar))) - 1) * 100)
median((exp(beta2_ch * diff(range(tr_temp_sd_season_mn / temp_sd_season_scalar))) - 1) * 100)
median((exp(beta3_ch * diff(range(tr_temp_sd_year_mn / temp_sd_year_scalar))) - 1) * 100)
median((exp(beta4_ch * diff(range(tr_precip_cv_season_mn / precip_cv_season_scalar))) - 1) * 100)
median((exp(beta5_ch * diff(range(tr_precip_cv_year_mn / precip_cv_year_scalar))) - 1) * 100)


# multiplicative effect of niche ------------------------------------------

gamma_ch <- MCMCvis::MCMCchains(fit, params = 'gamma')
gamma_rs_ch <- (exp(gamma_ch / y_scalar) - 1) * 100
apply(gamma_rs_ch, 2, mean)
apply(gamma_rs_ch, 2, function(x) quantile(x, probs = c(0.055, 0.5, 0.945)))
apply(gamma_rs_ch, 2, function(x) sum(x > 0) / NROW(gamma_rs_ch))

#pairwise differences niche categories
#k(k-1)/2
k <- (NCOL(gamma_ch)*(NCOL(gamma_ch) -1)) / 2
pw_df <- data.frame(nc1 = rep(NA, k),
           nc2 = NA,
           mean = NA,
           LCI = NA,
           UCI = NA,
           Pg0 = NA)
mn_mat <- matrix(NA, nrow = NCOL(gamma_ch), ncol = NCOL(gamma_ch))
Pg0_mat <- matrix(NA, nrow = NCOL(gamma_ch), ncol = NCOL(gamma_ch))
counter <- 1
for (i in 1:(NCOL(gamma_ch) - 1))
{
  #i <- 1
  for (j in (i+1):NCOL(gamma_ch))
  {
    #j <- 2
    d_t <- gamma_ch[,i] - gamma_ch[,j]
    d_rs <- (exp(d_t / y_scalar) - 1) * 100

    pw_df$nc1[counter] <- i
    pw_df$nc2[counter] <- j
    pw_df$mean[counter] <- mean(d_rs)
    pw_df$LCI[counter] <- quantile(d_rs, probs = 0.055)
    pw_df$UCI[counter] <- quantile(d_rs, probs = 0.945)
    pw_df$Pg0[counter] <- sum(d_rs > 0) / NROW(d_rs)

    #stack into upper triangular mat
    mn_mat[i, j] <- mean(d_rs)
    Pg0_mat[i, j] <- sum(d_rs > 0) / NROW(d_rs)
    counter <- counter + 1
  }
}

#join with niche category names
ref_df <- data.frame(num = 1:NCOL(gamma_ch),
           niche_names)

pw_df2 <- dplyr::left_join(pw_df, ref_df, by = c('nc1' = 'num')) %>%
  dplyr::left_join(ref_df, by = c('nc2' = 'num'))

saveRDS(pw_df2, paste0(fig_dir, '/pairwise_gamma.rds'))
# pw_df2 <- readRDS(paste0(fig_dir, '/pairwise_gamma.rds'))

#largest pairwise difference between categories
pw_df2[which.max(abs(pw_df2$mean)),]


# partial resid plots ------------------------------------------------

# https://www.wikiwand.com/en/Partial_residual_plot
pr_fun <- function(num, nm, YL = NULL, XL = NULL)
{
  names <- c('Mass', 'Temp seasonality', 'Temp interannual',
             'Precip seasonality', 'Precip interannual')
  
  #partial residuals
  # pr <- (resid / y_scalar) + (beta_mn[num] * DATA_X_sc[,num])
  pr <- matrix(NA, nrow = NROW(resid_ch), ncol = NCOL(resid_ch))
  if (num == 1)
  {
    for (i in 1:NROW(resid_ch))
    {
      # print(paste0('processing ', i, ' of ', NROW(resid_ch)))
      pr[i,] <- (resid_ch[i,]) + 
        (beta_ch[i,1] * DATA$lMass)
    }
    
    cov_sim <- seq(min(DATA$lMass), max(DATA$lMass), length.out = 15)
    out_sim <- matrix(NA, nrow = NROW(beta_ch), ncol = length(cov_sim))
    for (i in 1:NROW(beta_ch))
    {
      #i <- 1
      out_sim[i,] <- beta_ch[i,1] %*% cov_sim
    }
    pp_df <- data.frame(pr_mn = apply(pr, 2, mean) / y_scalar,
                        pr_sd = apply(pr, 2, sd) / y_scalar,
                        cov_mn = DATA$lMass / lMass_scalar)
    p2_df <- data.frame(x = cov_sim / lMass_scalar,
                        LCI = apply(out_sim, 2, function(x) quantile(x, probs = 0.055)) / 
                          y_scalar,
                        UCI = apply(out_sim, 2, function(x) quantile(x, probs = 0.945)) /
                          y_scalar,
                        mn = apply(out_sim, 2, mean) / y_scalar)
    
    pps <- ggplot() +
      # geom_errorbar(data = pp_df, 
      #               aes(x = cov_mn,
      #                   ymin = pr_mn - pr_sd, 
      #                   ymax = pr_mn + pr_sd),
      #               alpha = 0.05) +
      geom_point(data = pp_df, 
                 aes(cov_mn, pr_mn),
                 alpha = 0.15) +
      # geom_ribbon(data = p2_df,
      #             aes(x = x,
      #                 ymin = LCI,
      #                 ymax = UCI),
      #             fill = 'red',
      #             alpha = 0.2) +
      geom_line(data = p2_df,
                aes(x = x,
                    y = mn),
                col = 'red',
                alpha = 0.5,
                linewidth = 1.1) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(linewidth = 2),
            axis.ticks = element_line(linewidth = 1.5),
            axis.text.x = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 18))
    
    if (!is.null(YL))
    {
      pps <- pps +
        coord_cartesian(ylim =  YL)
    }
    if (!is.null(XL))
    {
      pps <- pps +
        coord_cartesian(xlim =  XL)
    }
    if (!is.null(YL) & !is.null(XL))
    {
      pps <- pps +
        coord_cartesian(ylim =  YL, xlim =  XL)
    }
  }
  
  plt_fun <- function(p_df, p2_df, num2)
  {
    #function for plotting
    pr_plt <- ggplot() +
      # geom_errorbar(data = p_df,
      #               aes(x = cov_mn,
      #                   ymin = pr_mn - pr_sd,
      #                   ymax = pr_mn + pr_sd),
      #               alpha = 0.05) +
      # geom_errorbarh(data = p_df,
      #                aes(y = pr_mn,
      #                    xmin = cov_mn - cov_sd,
      #                    xmax = cov_mn + cov_sd),
      #                alpha = 0.05) +
      geom_point(data = p_df, 
                 aes(cov_mn, pr_mn),
                 alpha = 0.2,
                 size = 1.6) +
      # geom_ribbon(data = p2_df,
      #             aes(x = x,
      #                 ymin = LCI,
      #                 ymax = UCI),
      #             fill = 'red',
      #             alpha = 0.2) +
      geom_line(data = p2_df,
                aes(x = x,
                    y = mn),
                col = 'red',
                alpha = 0.8,
                linewidth = 1.7) +
      geom_hline(yintercept = 0,
                 linetype = 'dashed',
                 # alpha = 0.3,
                 col = 'grey70',
                 linewidth = 1.7) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(linewidth = 2),
            axis.ticks = element_line(linewidth = 1.5),
            axis.text.x = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 18))
    
    if (!is.null(YL))
    {
      pr_plt <- pr_plt +
        coord_cartesian(ylim =  YL)
    }
    if (!is.null(XL))
    {
      pr_plt <- pr_plt +
        coord_cartesian(xlim =  XL)
    }
    if (!is.null(YL) & !is.null(XL))
    {
      pr_plt <- pr_plt +
        coord_cartesian(ylim =  YL, xlim =  XL)
    }
    return(pr_plt)
  }
  
  if (num == 2)
  {
    for (i in 1:NROW(resid_ch))
    {
      #i <- 1
      # print(paste0('processing ', i, ' of ', NROW(resid_ch)))
      pr[i,] <- (resid_ch[i,]) + 
        (beta_ch[i,2] * (tr_temp_sd_season_ch[i,]))
    }
    
    #get point (mean and sd)
    pp_df <- data.frame(pr_mn = apply(pr, 2, mean) / y_scalar,
                        pr_sd = apply(pr, 2, sd) / y_scalar,
                        cov_mn = apply(tr_temp_sd_season_ch, 2, mean) /
                          temp_sd_season_scalar,
                        cov_sd = apply(tr_temp_sd_season_ch, 2, sd) /
                          temp_sd_season_scalar)
    
    #sim beta over cov
    cov_sim <- seq(min(pp_df$cov_mn), max(pp_df$cov_mn), length.out = 15)
    out_sim <- matrix(NA, nrow = NROW(beta_ch), ncol = length(cov_sim))
    for (i in 1:NROW(beta_ch))
    {
      #i <- 1
      out_sim[i,] <- beta_ch[i,2] %*% (cov_sim * temp_sd_season_scalar)
    }
    p2_df <- data.frame(x = cov_sim,
                        LCI = apply(out_sim, 2, function(x) quantile(x, probs = 0.055)) / 
                          y_scalar,
                        UCI = apply(out_sim, 2, function(x) quantile(x, probs = 0.945)) /
                          y_scalar,
                        mn = apply(out_sim, 2, mean) /
                          y_scalar)
    
    #run function to plot
    pps <- plt_fun(pp_df, p2_df, num2 = 2)
  }
  if (num == 3)
  {
    for (i in 1:NROW(resid_ch))
    {
      #i <- 1
      # print(paste0('processing ', i, ' of ', NROW(resid_ch)))
      pr[i,] <- (resid_ch[i,]) + 
        (beta_ch[i,3] * (tr_temp_sd_year_ch[i,]))
    }
    
    #get point (mean and sd)
    pp_df <- data.frame(pr_mn = apply(pr, 2, mean) / y_scalar,
                        pr_sd = apply(pr, 2, sd) / y_scalar,
                        cov_mn = apply(tr_temp_sd_year_ch, 2, mean) /
                          temp_sd_year_scalar,
                        cov_sd = apply(tr_temp_sd_year_ch, 2, sd) /
                          temp_sd_year_scalar)
    #sim beta over cov
    cov_sim <- seq(min(pp_df$cov_mn), max(pp_df$cov_mn), length.out = 15)
    out_sim <- matrix(NA, nrow = NROW(beta_ch), ncol = length(cov_sim))
    for (i in 1:NROW(beta_ch))
    {
      #i <- 1
      out_sim[i,] <- beta_ch[i,3] %*% (cov_sim * temp_sd_year_scalar)
    }
    p2_df <- data.frame(x = cov_sim,
                        LCI = apply(out_sim, 2, function(x) quantile(x, probs = 0.055)) / 
                          y_scalar,
                        UCI = apply(out_sim, 2, function(x) quantile(x, probs = 0.945)) / 
                          y_scalar,
                        mn = apply(out_sim, 2, mean) / 
                          y_scalar)
    
    #run function to plot
    pps <- plt_fun(pp_df, p2_df, num2 = 3)
  }
  if (num == 4)
  {
    for (i in 1:NROW(resid_ch))
    {
      #i <- 1
      # print(paste0('processing ', i, ' of ', NROW(resid_ch)))
      pr[i,] <- (resid_ch[i,]) + 
        (beta_ch[i,4] * (tr_precip_cv_season_ch[i,]))
    }
    
    #get point (mean and sd)
    pp_df <- data.frame(pr_mn = apply(pr, 2, mean) / y_scalar,
                        pr_sd = apply(pr, 2, sd) / y_scalar,
                        cov_mn = apply(tr_precip_cv_season_ch, 2, mean) /
                          precip_cv_season_scalar,
                        cov_sd = apply(tr_precip_cv_season_ch, 2, sd) /
                          precip_cv_season_scalar)
    #sim beta over cov
    cov_sim <- seq(min(pp_df$cov_mn), max(pp_df$cov_mn), length.out = 15)
    out_sim <- matrix(NA, nrow = NROW(beta_ch), ncol = length(cov_sim))
    for (i in 1:NROW(beta_ch))
    {
      #i <- 1
      out_sim[i,] <- beta_ch[i,4] %*% (cov_sim * precip_cv_season_scalar)
    }
    p2_df <- data.frame(x = cov_sim,
                        LCI = apply(out_sim, 2, function(x) quantile(x, probs = 0.055)) / 
                          y_scalar,
                        UCI = apply(out_sim, 2, function(x) quantile(x, probs = 0.945)) / 
                          y_scalar,
                        mn = apply(out_sim, 2, mean) / y_scalar)
    
    #run function to plot
    pps <- plt_fun(pp_df, p2_df, num2 = 4)
  }
  if (num == 5)
  {
    for (i in 1:NROW(resid_ch))
    {
      #i <- 1
      # print(paste0('processing ', i, ' of ', NROW(resid_ch)))
      pr[i,] <- (resid_ch[i,]) + 
        (beta_ch[i,5] * (tr_precip_cv_year_ch[i,]))
    }
    #get point (mean and sd)
    pp_df <- data.frame(pr_mn = apply(pr, 2, mean) / y_scalar,
                        pr_sd = apply(pr, 2, sd) / y_scalar,
                        cov_mn = apply(tr_precip_cv_year_ch, 2, mean) /
                          precip_cv_year_scalar,
                        cov_sd = apply(tr_precip_cv_year_ch, 2, sd) /
                          precip_cv_year_scalar)
    #sim beta over cov
    cov_sim <- seq(min(pp_df$cov_mn), max(pp_df$cov_mn), length.out = 15)
    out_sim <- matrix(NA, nrow = NROW(beta_ch), ncol = length(cov_sim))
    for (i in 1:NROW(beta_ch))
    {
      #i <- 1
      out_sim[i,] <- beta_ch[i,5] %*% (cov_sim * precip_cv_year_scalar)
    }
    p2_df <- data.frame(x = cov_sim,
                        LCI = apply(out_sim, 2, function(x) quantile(x, probs = 0.055)) / 
                          y_scalar,
                        UCI = apply(out_sim, 2, function(x) quantile(x, probs = 0.945)) / 
                          y_scalar,
                        mn = apply(out_sim, 2, mean) / y_scalar)
    
    #run function to plot
    pps <- plt_fun(pp_df, p2_df, num2 = 5)
  }
  
  ggsave(filename = paste0(fig_dir, nm, '-pr-', run_date, '.pdf'),
         pps,
         width = 6,
         height = 4.5)
  # print(pps)
}

pr_fun(num = 1, nm = 'mass-neb') #Mass
pr_fun(num = 2, nm = 'temp-season-neb', #temp season
       YL = c(-0.6, 0.65)) #excludes one outlier
pr_fun(num = 3, nm = 'temp-year-neb', #temp year
       YL = c(-0.6, 0.65)) #excludes one outlier
pr_fun(num = 4, nm = 'precip-season-neb') #precip season
pr_fun(num = 5, nm = 'precip-year-neb') #precip year


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
MCMCvis::MCMCplot(cbind(beta3_rs_ch,
                        beta2_rs_ch,
                        beta5_rs_ch,
                        beta4_rs_ch),
                  labels = c('T interannual var',
                             'T seasonality',
                             'P interannual var',
                             'P seasonality'),
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  xlim = c(-3,2),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = '% change GL for 1 sd change in cov',
                  guide_lines = TRUE)
dev.off()


pdf(paste0(fig_dir, 'gamma-cat-rs-', run_date, '.pdf'),
    height = 5, width = 5)
MCMCvis::MCMCplot(gamma_rs_ch,
                  labels = niche_names,
                  rank = TRUE,
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  xlim = c(-25, 25),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = '% change GL by niche group',
                  guide_lines = TRUE)
dev.off()


# PPC ---------------------------------------------------------------------

# PPC - normal model
sigma_ch <- MCMCvis::MCMCchains(fit, params = 'sigma')

#500 iterations
sidx <- sample(1:NROW(sigma_ch), size = 500)
mu_rep <- matrix(NA, nrow = length(sidx), ncol = length(DATA$Y_obs))
y_rep <- matrix(NA, nrow = length(sidx), ncol = length(DATA$Y_obs))
yt_rep <- matrix(NA, nrow = length(sidx), ncol = length(DATA$Y_obs))
for (i in 1:length(sidx))
{
  #i <- 1
  print(paste0('iter: ', i, ' of ', length(sidx)))
  for (j in 1:length(DATA$Y_obs))
  {
    #j <- 1
    mu_rep[i, j] <- kappa_ch[sidx[i], 1] + 
      gamma_ch[sidx[i], DATA$niche_idx[j]] +
      alpha_ch[sidx[i], j] + 
      beta_ch[i,1] * DATA$lMass[j] +
      beta_ch[i,2] * tr_temp_sd_season_ch[i,j] +
      beta_ch[i,3] * tr_temp_sd_year_ch[i,j] +
      beta_ch[i,4] * tr_precip_cv_season_ch[i,j] +
      beta_ch[i,5] * tr_precip_cv_year_ch[i,j]
    
    yt_rep[i,j] <- rnorm(1, mu_rep[i, j], sigma_ch[sidx[i], 1])
    y_rep[i,j] <- rnorm(1, yt_rep[i,j], DATA$sd_Y)
  }
}

pdf(paste0(fig_dir, 'PPC-', run_date, '.pdf'), height = 5, width = 5)
plot(density(DATA$Y / y_scalar), col = 'black', lwd = 3, ylim = c(0, 1.5))#, xlim = c(0, 3.5))
for (i in 1:500)
{
  lines(density(y_rep[i,] / y_scalar), col = rgb(1,0,0,0.05))
}
dev.off()


# VIF ---------------------------------------------------------------------

# #covariates as a function of other covariates
# tf1 <- lm(DATA$lMass ~ tr_temp_sd_season_mn + tr_temp_sd_year_mn +
#             tr_precip_cv_season_mn + tr_precip_cv_year_mn)
# stf1 <- summary(tf1)
# 
# tf2 <- lm(tr_temp_sd_season_mn ~ DATA$lMass + tr_temp_sd_year_mn +
#             tr_precip_cv_season_mn + tr_precip_cv_year_mn)
# stf2 <- summary(tf2)
# 
# tf3 <- lm(tr_temp_sd_year_mn ~  DATA$lMass + tr_temp_sd_season_mn +
#             tr_precip_cv_season_mn + tr_precip_cv_year_mn)
# stf3 <- summary(tf3)
# 
# tf4 <- lm(tr_precip_cv_season_mn ~ DATA$lMass + tr_temp_sd_season_mn +
#             tr_temp_sd_year_mn + tr_precip_cv_year_mn)
# stf4 <- summary(tf4)
# 
# tf5 <- lm(tr_precip_cv_year_mn ~ DATA$lMass + tr_temp_sd_season_mn +
#             tr_temp_sd_year_mn + tr_precip_cv_season_mn)
# stf5 <- summary(tf5)
# 
# 
# #calc VIF per covariate
# 1 / (1 - stf1$r.squared) #lMass
# 1 / (1 - stf2$r.squared) #temp_sd_year
# 1 / (1 - stf3$r.squared) #temp_sd_season
# 1 / (1 - stf4$r.squared) #precip_cv_year
# 1 / (1 - stf5$r.squared) #precip_cv_season
