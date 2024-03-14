# TITLE:            Fit bayesian model 
# PROJECT:          IBEEM Environmental Variability and Life History
# AUTHORS:          Casey Youngflesh, Peter Williams, Jeff Dozer, Adriana Uscanga, Lala Kounta, Kelly Kapsar, Phoebe Zarnetske, Pat Bills
# DATA INPUT:       Bird consensus tree (4b), 
# DATA OUTPUT:      Model output summary and figures (caterpillar plots, partial residual plots, etc.)
# DATE:             October 2023 
# OVERVIEW:         Fits bayesian model to data (gl ~ env + phylo + vint)


# specify dir -------------------------------------------------------------

dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
sc_dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
# dir <- '/mnt/research/ibeem/variability/'
# sc_dir <- '/mnt/home/ccy/variability/'
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

# #exclude Sicalis auriventris - outlier precip_cv_season
# s_auriventris <- read.csv(paste0(dir, 'data/L3/main-bird-data-birdtree2.csv')) %>%
#   dplyr::filter(Migration == 1, precip_cv_season >= 2.5) %>%
#   dplyr::select(precip_cv_season)
# 
# all_precip_cv_season <- read.csv(paste0(dir, 'data/L3/main-bird-data-birdtree2.csv')) %>%
#   dplyr::filter(Order %ni% or_excl,
#                 Migration == 1) %>%
#   dplyr::select(precip_cv_season) %>%
#   as.vector()
# 
# #11.6 MAD
# (s_auriventris$precip_cv_season - median(all_precip_cv_season$precip_cv_season)) / mad(all_precip_cv_season$precip_cv_season)


'%ni%' <- Negate('%in%')
bird_df <- read.csv(paste0(dir, 'data/L3/main-bird-data-birdtree2.csv')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::filter(Order %ni% or_excl,
                Migration == 1,
                #filter out precip outlier
                precip_cv_season < 2.5) %>% 
  dplyr::mutate(lMass = log(Mass),
                lGL = log(GenLength),
                species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  dplyr::rename(Trophic_niche = Trophic.Niche) %>%
  #drop duplicated species
  dplyr::group_by(species) %>%
  dplyr::slice_head() %>%
  dplyr::ungroup() %>%
  #how much var will change in 1 generation (in sds)
  #((degrees / year) * (year / generation) * (1 / degrees) = degrees (sd) / generation
  dplyr::mutate(temp_delta = (temp_slope * GenLength) / temp_sd_year, 
                precip_delta = (precip_slope * GenLength) / precip_sd_year)

#subset out just traits of interest
bird_df2 <- dplyr::select(bird_df, 
                          ID,
                          species,
                          Order,
                          Family,
                          lGL,
                          Trophic_niche,
                          lMass,
                          temp_sd_year,
                          temp_sd_season,
                          precip_cv_year,
                          precip_cv_season,
                          temp_delta,
                          precip_delta)


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

#get estimate of Pagel's kappa to scale phylogeny
#https://lukejharmon.github.io/pcm/chapter6_beyondbm/
fit_ka <- geiger::fitContinuous(pr_tree2, dd[,'lGL'], model = "kappa")

#other transformations - lambda, delta, OU, Brownian motion
fit_la <- geiger::fitContinuous(pr_tree2, dd[,'lGL'], model = "lambda")
fit_de <- geiger::fitContinuous(pr_tree2, dd[,'lGL'], model = "delta")
fit_ou <- geiger::fitContinuous(pr_tree2, dd[,'lGL'], model = "OU")
fit_bm <- geiger::fitContinuous(pr_tree2, dd[,'lGL'], model = "BM")

aic_df <- data.frame(method = c('kappa', 'lambda', 'delta', 'ou', 'bm'),
                     AIC = c(fit_ka$opt$aic, fit_la$opt$aic, fit_de$opt$aic, 
                             fit_ou$opt$aic, fit_bm$opt$aic))

#rescale tree
pr_tree_k <- phytools::rescale(pr_tree, 'kappa', 
                     kappa = fit_ka$opt$kappa, sigsq = fit_ka$opt$sigsq)

#get corr matrix of rescaled tree
Rho <- ape::vcv.phylo(pr_tree_k, corr = TRUE)


# niche levels ------------------------------------------------------------

#assign missing species (owls) to Vertivore
bird_df3$Trophic_niche[which(is.na(bird_df3$Trophic_niche))] <- "Vertivore"

bird_df3$niche_idx <- as.numeric(factor(bird_df3$Trophic_niche))
niche_names <- levels(factor(bird_df3$Trophic_niche))


# summary metrics ---------------------------------------------------------

#7476 species
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
temp_sd_year_scalar <- 0.1
precip_cv_season_scalar <- 0.1
precip_cv_year_scalar <- 0.5
y_scalar <- 2

#center predictors
tt <- data.frame(lMass = bird_df3$lMass * 
                   lMass_scalar,
                 temp_sd_season = bird_df3$temp_sd_season * 
                   temp_sd_season_scalar,
                 temp_sd_year = bird_df3$temp_sd_year * 
                   temp_sd_year_scalar,
                 precip_cv_season = bird_df3$precip_cv_season * 
                   precip_cv_season_scalar,
                 precip_cv_year = bird_df3$precip_cv_year * 
                   precip_cv_year_scalar) %>%
  apply(2, function(x) scale(x, scale = FALSE)[,1])



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


# fit model ---------------------------------------------------------------

DATA <- list(N = NROW(bird_df3),
             Y = bird_df3$lGL * y_scalar,
             K = NCOL(tt),
             J = length(unique(bird_df3$niche_idx)),
             X = tt,
             niche_idx = bird_df3$niche_idx,
             mu_kappa = 2,
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
# STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 4000

#compile model
mod <- cmdstanr::cmdstan_model(paste0(sc_dir, 'Scripts/Model_files/5-phylo-vint.stan'),
                               force_recompile = TRUE)

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
                  file_name = paste0('bird-gl-phylo-vint-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('bird-gl-phylo-vint-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('bird-gl-phylo-vint-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('bird-gl-phylo-vint-data-', run_date),
                  cp_file = c(paste0(sc_dir, 'Scripts/Model_files/5-phylo-vint.stan'), 
                              paste0(sc_dir, 'Scripts/5-model/5-gl-phylo-vint.R')),
                  cp_file_names = c(paste0('5-phylo-vint-', run_date, '.stan'),
                                    paste0('5-gl-phylo-vint-', run_date, '.R')))

fig_dir <- paste0(dir, 'Results/bird-gl-phylo-vint-', run_date, '/')

# fit <- readRDS(paste0(dir, '/Results/bird-gl-phylo-vint-', run_date,
#                       '/bird-gl-phylo-vint-fit-', run_date, '.rds'))
# library(shinystan)
# shinystan::launch_shinystan(fit)


# residuals ---------------------------------------------------------------

alpha_mn <- MCMCvis::MCMCpstr(fit, params = 'alpha')[[1]]
kappa_mn <- MCMCvis::MCMCpstr(fit, params = 'kappa')[[1]]
gamma_mn <- MCMCvis::MCMCpstr(fit, params = 'gamma')[[1]]
beta_mn <- MCMCvis::MCMCpstr(fit, params = 'beta')[[1]]

mu_mn <- kappa_mn + gamma_mn[DATA$niche_idx] + alpha_mn + (DATA$X %*% beta_mn)[,1]
resid <- DATA$Y - mu_mn


# Summary -----------------------------------------------------------------

#model summary
MCMCvis::MCMCsummary(fit, round = 3, 
                     params = c('beta', 'kappa', 'gamma', 'sigma', 'sigma_phylo'),
                     pg0 = TRUE)


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

beta_mn_sc <- apply(cbind(beta1_ch, beta2_ch, beta3_ch, beta4_ch, beta5_ch), 2, mean)
DATA_X_sc <- DATA$X
DATA_X_sc[,] <- NA
DATA_X_sc[,1] <- DATA$X[,1] / lMass_scalar
DATA_X_sc[,2] <- DATA$X[,2] / temp_sd_season_scalar
DATA_X_sc[,3] <- DATA$X[,3] / temp_sd_year_scalar
DATA_X_sc[,4] <- DATA$X[,4] / precip_cv_season_scalar
DATA_X_sc[,5] <- DATA$X[,5] / precip_cv_year_scalar


#scaling cov to measured scale bc transformed param est
#% change in LH trait for 1 sd change in covariate
beta1_rs_ch <- (exp(beta1_ch * sd(tt[,1] / lMass_scalar)) - 1) * 100
beta2_rs_ch <- (exp(beta2_ch * sd(tt[,2] / temp_sd_season_scalar)) - 1) * 100
beta3_rs_ch <- (exp(beta3_ch * sd(tt[,3] / temp_sd_year_scalar)) - 1) * 100
beta4_rs_ch <- (exp(beta4_ch * sd(tt[,4] / precip_cv_season_scalar)) - 1) * 100
beta5_rs_ch <- (exp(beta5_ch * sd(tt[,5] / precip_cv_year_scalar)) - 1) * 100

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
sum(beta5_rs_ch < 0) / NROW(beta5_rs_ch)

median((exp(beta1_ch * diff(range(tt[,1] / lMass_scalar))) - 1) * 100)
median((exp(beta2_ch * diff(range(tt[,2] / temp_sd_season_scalar))) - 1) * 100)
median((exp(beta3_ch * diff(range(tt[,3] / temp_sd_year_scalar))) - 1) * 100)
median((exp(beta4_ch * diff(range(tt[,4] / precip_cv_season_scalar))) - 1) * 100)
median((exp(beta5_ch * diff(range(tt[,5] / precip_cv_year_scalar))) - 1) * 100)


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

# saveRDS(pw_df2, paste0(fig_dir, '/pairwise_gamma.rds'))
# pw_df2 <- readRDS(paste0(fig_dir, '/pairwise_gamma.rds'))

#largest pairwise difference between categories
pw_df2[which.max(abs(pw_df2$mean)),]


# partial resid plots ------------------------------------------------

# https://www.wikiwand.com/en/Partial_residual_plot
pr_fun <- function(num, nm, YL = NULL)
{
  names <- c('Mass', 'Temp seasonality', 'Temp interannual',
             'Precip seasonality', 'Precip interannual')
  
  #partial residuals
  pr <- (resid / y_scalar) + (beta_mn_sc[num] * DATA_X_sc[,num])
  
  #range(pr) #mass: -0.698  1.710
  #range(pr) #temp season: -0.0438  0.3780
  #range(pr) #temp year: -0.4113  0.3652
  #range(pr) #precip season: -0.4112  0.3708
  #range(pr) #precip year: -0.4170  0.3757
  #range(c(pr2, pr5)): -0.4170  0.3780
  
  pdf(paste0(fig_dir, nm, '-pr-', run_date, '.pdf'),
      height = 5, width = 5)
  plot(DATA_X_sc[,num], pr, col = rgb(0,0,0,0.2), pch = 19,
       xlab = 'Predictor',
       ylab = 'Partial residual',
       main = names[num],
       ylim = YL)
  abline(h = 0, col = 'grey', lwd = 4, lty = 2)
  abline(a = 0, b = beta_mn_sc[num], col = rgb(1,0,0,0.5), lwd = 4)
  dev.off()
}

pr_fun(num = 1, nm = 'mass') #Mass
pr_fun(num = 2, nm = 'temp-season', YL = c(-0.4170, 0.3780)) #temp season
pr_fun(num = 3, nm = 'temp-year') #temp year
pr_fun(num = 4, nm = 'precip-season') #precip season
pr_fun(num = 5, nm = 'precip-year', YL = c(-0.4170, 0.3780)) #precip year


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
                  sz_thick = 3,
                  sz_thin = 3,
                  main = '% change GL by niche group',
                  guide_lines = TRUE)
dev.off()


# PPC ---------------------------------------------------------------------

# PPC - normal model
sigma_ch <- MCMCvis::MCMCchains(fit, params = 'sigma')
alpha_ch <- MCMCvis::MCMCchains(fit, params = 'alpha')
kappa_ch <- MCMCvis::MCMCchains(fit, params = 'kappa')
gamma_ch <- MCMCvis::MCMCchains(fit, params = 'gamma')
beta_ch <- MCMCvis::MCMCchains(fit, params = 'beta')

#500 iterations
sidx <- sample(1:NROW(sigma_ch), size = 500)
mu_rep <- matrix(NA, nrow = length(sidx), ncol = length(DATA$Y))
y_rep <- matrix(NA, nrow = length(sidx), ncol = length(DATA$Y))
for (i in 1:length(sidx))
{
  #i <- 1
  print(paste0('iter: ', i, ' of ', length(sidx)))
  for (j in 1:length(DATA$Y))
  {
    #j <- 1
    #t-dis
    # eps <- rt(n = 1, df = nu_ch[sidx[i],1]) * sigma_ch[sidx[i],1]
    # y_rep[i,j] <- mu_ch[sidx[i],j] + eps
    
    mu_rep[i, j] <- kappa_ch[sidx[i], 1] + gamma_ch[sidx[i], DATA$niche_idx[j]] + 
      alpha_ch[sidx[i], j] + (DATA$X[j,] %*% beta_ch[sidx[i],])[,1]
    
    y_rep[i,j] <- rnorm(1, mu_rep[i, j], sigma_ch[sidx[i], 1])
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

#covariates as a function of other covariates
# tf1 <- lm(lMass ~ temp_sd_year + temp_sd_season +
#             precip_cv_year + precip_cv_season, 
#           data = bird_df3)
# stf1 <- summary(tf1)
# 
# tf2 <- lm(temp_sd_year ~ lMass + temp_sd_season +
#             precip_cv_year + precip_cv_season, data = 
#             bird_df3)
# stf2 <- summary(tf2)
# 
# tf3 <- lm(temp_sd_season ~ lMass + temp_sd_year +
#             precip_cv_year + precip_cv_season, 
#           data = bird_df3)
# stf3 <- summary(tf3)
# 
# tf4 <- lm(precip_cv_year ~ lMass + temp_sd_year + temp_sd_season + 
#             precip_cv_season, 
#           data = bird_df3)
# stf4 <- summary(tf4)
# 
# tf5 <- lm(precip_cv_season ~ lMass + temp_sd_year + temp_sd_season + 
#             precip_cv_year, 
#           data = bird_df3)
# stf5 <- summary(tf5)
# 
# 
# #calc VIF per covariate
# 1 / (1 - stf1$r.squared) #lMass
# 1 / (1 - stf2$r.squared) #temp_sd_year
# 1 / (1 - stf3$r.squared) #temp_sd_season
# 1 / (1 - stf4$r.squared) #precip_cv_year
# 1 / (1 - stf5$r.squared) #precip_cv_season
