####################
# Fit Bayes model - clutch size ~ env + phylo + oe + varyiny int by niche
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
  #sample a subset of species
  # dplyr::slice_sample(n = 2000) %>%
  dplyr::mutate(fac_Family = factor(Family),
                fac_Order = factor(Order),
                lMass = log(Mass),
                lGL = log(GenLength),
                f_id = as.numeric(fac_Family)) %>%
  #drop duplicated species (for now)
  dplyr::group_by(Birdtree_name) %>%
  dplyr::slice_head() %>%
  dplyr::ungroup()

#survival etc from Bird et al.
bs <- read.csv(paste0(dir, 'data/L1/trait/bird-et-al-data-with-id.csv')) %>%
  dplyr::select(-Order, -Family, -GenLength)

#clutch size from Bird et al.
bcs <- read.csv(paste0(dir, 'data/L1/trait/bird_et_al_clutch_size.csv'))

#join
bird_df2 <- dplyr::left_join(bird_df, bs, by = c('Birdtree_name' = 'Sci_name')) %>%
  dplyr::arrange(Birdtree_name) %>%
  dplyr::mutate(Scientific.name = stringr::str_to_sentence(Birdtree_name),
                #correct for incorrect scaling in data
                Measured_survival = Measured_survival * 0.01) %>%
  dplyr::left_join(dplyr::select(bcs, -Order, -Family, -Genus), 
                   by = 'Scientific.name') %>%
  dplyr::select(ID, Birdtree_name, Avonet_name, Family, Order, Trophic.Niche,
                temp_mean, temp_sd_year, temp_sd_season, temp_sp_color_month,
                precip_mean, precip_cv_year, precip_cv_season, precip_sp_color_month,
                Mass, GenLength, lMass, lGL, Mean.clutch.size, Measured_survival, 
                Measured_age_first_breeding, Measured_max_longevity, 
                Modeled_survival, Modeled_age_first_breeding, Modeled_max_longevity)

#subset out just traits of interest
tri <- dplyr::mutate(bird_df2,
                     species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  dplyr::mutate(Measured_log_age_first_breeding = log(Measured_age_first_breeding),
                Measured_log_max_longevity = log(Measured_max_longevity),
                Measured_log_clutch_size = log(Mean.clutch.size)) %>%
  dplyr::select(species, 
                lMass,
                Measured_survival,
                Measured_log_age_first_breeding,
                Measured_log_max_longevity,
                Measured_log_clutch_size)


# trait imputation --------------------------------------------------------

#load consensus tree - bird.phylo
load(paste0(dir, 'data/L3/bird-consensus-tree.rda'))

#df with names and idx
idx_df <- data.frame(idx = 1:NROW(tri),
                     name = tri$species)

#species not found in both datasets (species to drop from tree)
nm <- setdiff(bird.phylo$tip.label, tri$species)

#prune specified tips from tree
pr_tree <- ape::drop.tip(bird.phylo, nm)

#get index for name order on tips
j_idx <- dplyr::left_join(data.frame(name = pr_tree$tip.label), idx_df, 
                          by = 'name')


# #run phylo imputation - not working properly on cluster - read in bird_df3 instead
# ir <- Rphylopars::phylopars(trait_data = tri[j_idx$idx,],
#                             tree = pr_tree,
#                             phylo_correlated = TRUE,
#                             # model = "BM") # AIC = ???
#                             # model = "OU") # AIC = ???
#                             model = "lambda") # AIC = 10106.53
# # model = "mvOU") # threw an error
# # model = "delta") # threw an error
# # model = "EB") # AIC = ???
# # model = "star") # ???
# 
# # AIC(ir)
# 
# #returns species means as well as some internal nodes
# #get just species rows
# ridx <- which(row.names(ir$anc_rec) %in% tri$species)
# 
# #variance (uncertainty)
# ir_unc <- data.frame(species = tri$species[j_idx$idx],
#                      ir$anc_var[ridx,]) %>%
#   dplyr::mutate(SD_survival = sqrt(Measured_survival),
#                 SD_log_age_first_breeding = sqrt(Measured_log_age_first_breeding),
#                 SD_log_max_longevity = sqrt(Measured_log_max_longevity),
#                 SD_log_clutch_size = sqrt(Measured_log_clutch_size)) %>%
#   dplyr::arrange(species) %>%
#   dplyr::select(species,
#                 SD_survival,
#                 SD_log_age_first_breeding,
#                 SD_log_max_longevity,
#                 SD_log_clutch_size)
# 
# #merge imputed values with unc and env data
# ir_mrg <- data.frame(species = tri$species[j_idx$idx],
#                      ir$anc_rec[ridx,]) %>%
#   dplyr::arrange(species) %>%
#   dplyr::rename(Phylo_survival = Measured_survival,
#                 Phylo_log_age_first_breeding = Measured_log_age_first_breeding,
#                 Phylo_log_max_longevity = Measured_log_max_longevity,
#                 Phylo_log_clutch_size = Measured_log_clutch_size) %>%
#   #join with measured traits
#   dplyr::left_join(dplyr::select(tri, -lMass), by = 'species') %>%
#   #join with var
#   dplyr::left_join(ir_unc, by = 'species')
# row.names(ir_mrg) <- NULL
# 
# #join trait and env data
# bird_df3 <- dplyr::mutate(bird_df2,
#                           species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
#   dplyr::select(species,
#                 temp_mean,
#                 temp_sd_year,
#                 temp_sd_season,
#                 temp_sp_color_month,
#                 precip_cv_year,
#                 precip_cv_season,
#                 precip_sp_color_month,
#                 Trophic_niche = Trophic.Niche,
#                 lGL,
#                 Modeled_survival,
#                 Modeled_age_first_breeding,
#                 Modeled_max_longevity) %>%
#   dplyr::left_join(ir_mrg, by = 'species')
# 
# #assign missing species (owls) to Vertivore
# bird_df3$Trophic_niche[which(is.na(bird_df3$Trophic_niche))] <- "Vertivore"

# saveRDS(bird_df3, paste0(dir, 'Scripts/5-model/bird_df3.rds'))
bird_df3 <- readRDS(paste0(sc_dir, 'Scripts/5-model/bird_df3.rds'))


# phylo -------------------------------------------------------------------

#subset of imp
# set.seed(1)
# Nsel <- 500
# stidx <- sample(x = 1:NROW(bird_df3), size = Nsel)
# bird_df4 <- bird_df3[stidx,]
bird_df4 <- bird_df3

# #prune tree
# #df with names and idx
# idx_df2 <- data.frame(idx = 1:NROW(bird_df4),
#                       name = stringr::str_to_title(gsub(' ', '_', bird_df4$species)))
# 
# #species not found in both datasets (species to drop from tree)
# nm2 <- setdiff(pr_tree$tip.label, bird_df4$species)
# 
# #prune specified tips from tree
# pr_tree2 <- ape::drop.tip(pr_tree, nm2)

#get idx
j_idx3 <- dplyr::left_join(data.frame(species = pr_tree$tip.label), 
                           data.frame(idx = 1:NROW(bird_df4), bird_df4), 
                           by = 'species')

#apply
bird_df5 <- bird_df4[j_idx3$idx,]

#niche levels
bird_df5$niche_idx <- as.numeric(factor(bird_df5$Trophic_niche))
niche_names <- levels(factor(bird_df5$Trophic_niche))

#separate obs and imp values
obs_idx <- which(bird_df5$SD_log_clutch_size == 0)
imp_idx <- which(bird_df5$SD_log_clutch_size != 0)

#get corr matrix
# V <- ape::vcv.phylo(pr_tree2, corr = TRUE)
V <- ape::vcv.phylo(pr_tree, corr = TRUE)


# scale/prep data ---------------------------------------------------------

#scalars for data
lMass_scalar <- 0.5
temp_sd_season_scalar <- 0.1
temp_sd_year_scalar <- 0.05
precip_cv_season_scalar <- 0.05
precip_cv_year_scalar <- 0.05
y_scalar <- 10

#split predictors into obs and imp
tt_obs <- data.frame(lMass = bird_df5$lMass[obs_idx] * 
                       lMass_scalar,
                     temp_sd_season = bird_df5$temp_sd_season[obs_idx] * 
                       temp_sd_season_scalar,
                     temp_sd_year = bird_df5$temp_sd_year[obs_idx] * 
                       temp_sd_year_scalar,
                     precip_cv_season = bird_df5$precip_cv_season[obs_idx] * 
                       precip_cv_season_scalar,
                     precip_cv_year = bird_df5$precip_cv_year[obs_idx] * 
                       precip_cv_year_scalar)

tt_imp <- data.frame(lMass = bird_df5$lMass[imp_idx] * 
                       lMass_scalar,
                     temp_sd_season = bird_df5$temp_sd_season[imp_idx] * 
                       temp_sd_season_scalar,
                     temp_sd_year = bird_df5$temp_sd_year[imp_idx] * 
                       temp_sd_year_scalar,
                     precip_cv_season = bird_df5$precip_cv_season[imp_idx] * 
                       precip_cv_season_scalar,
                     precip_cv_year = bird_df5$precip_cv_year[imp_idx] * 
                       precip_cv_year_scalar)

#subtract off mean to center vars
tt_mns <- rbind(tt_obs, tt_imp) %>%
  apply(2, mean)

tt_obs2 <- sweep(tt_obs, 2, tt_mns)
tt_imp2 <- sweep(tt_imp, 2, tt_mns)


# fit model ---------------------------------------------------------------

DATA <- list(N = NROW(bird_df5),
             No = length(obs_idx),
             Ni = length(imp_idx),
             y_obs = bird_df5$Phylo_log_clutch_size[obs_idx] * y_scalar,
             y_imp = bird_df5$Phylo_log_clutch_size[imp_idx] * y_scalar,
             sd_y = bird_df5$SD_log_clutch_size[imp_idx] * y_scalar,
             K = NCOL(tt_obs),
             J = length(unique(bird_df5$niche_idx)),
             X_obs = tt_obs2,
             X_imp = tt_imp2,
             obs_idx = obs_idx,
             imp_idx = imp_idx,
             niche_obs_idx = bird_df5$niche_idx[obs_idx],
             niche_imp_idx = bird_df5$niche_idx[imp_idx],
             LRho = chol(V)) #cholesky factor of corr matrix

options(mc.cores = parallel::detectCores())

# DELTA <- 0.92
# TREE_DEPTH <- 12
# STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 2000

#compile model
mod <- cmdstanr::cmdstan_model(paste0(sc_dir, 'Scripts/Model_files/5-phylo-oe-vint.stan'))

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
                  file_name = paste0('bird-clutch-phylo-oe-vint-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('bird-clutch-phylo-oe-vint-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('bird-clutch-phylo-oe-vint-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('bird-clutch-phylo-oe-vint-data-', run_date),
                  cp_file = c(paste0(sc_dir, 'Scripts/Model_files/5-phylo-oe-vint.stan'), 
                              paste0(sc_dir, 'Scripts/5-model/5-clutch-phylo-oe-vint.R')),
                  cp_file_names = c(paste0('5-phylo-oe-vint-', run_date, '.stan'),
                                    paste0('5-clutch-phylo-oe-vint-', run_date, '.R')))

fig_dir <- paste0(dir, 'Results/bird-clutch-phylo-oe-vint-', run_date)

# library(shinystan)
# shinystan::launch_shinystan(fit)
# fit <- readRDS(paste0(dir, '/Results/bird-surv-phylo-oe-', run_date,
#                       '/bird-surv-phylo-oe-fit-', run_date, '.rds'))
# DATA <- readRDS(paste0(dir, '/Results/bird-surv-phylo-oe-', run_date,
#                       '/bird-surv-phylo-oe-data-', run_date, '.rds'))

# 
# # residuals ---------------------------------------------------------------
# 
# # extract residuals and calc phylo signal
# mu_obs_mn <- MCMCvis::MCMCpstr(fit, params = 'mu_obs')[[1]]
# 
# #calc resid
# resid <- DATA$y_obs - mu_obs_mn

beta_mn <- MCMCvis::MCMCpstr(fit, params = 'beta')[[1]]
kappa_mn <- MCMCvis::MCMCpstr(fit, params = 'kappa')[[1]]
gamma_mn <- MCMCvis::MCMCpstr(fit, params = 'gamma')[[1]]
alpha_mn <- MCMCvis::MCMCpstr(fit, params = 'alpha')[[1]]

mu_obs <- kappa_mn + gamma_mn[niche_obs_idx] + alpha_mn[DATA$obs_idx] + as.matrix(DATA$X_obs) %*% beta_mn
mu_imp <- kappa_mn + gamma_mn[niche_imp_idx] + alpha_mn[DATA$imp_idx] + as.matrix(DATA$X_imp) %*% beta_mn

resid_obs <- DATA$y_obs - mu_obs
resid_imp <- DATA$y_imp - mu_imp


# Summary -----------------------------------------------------------------

#model summary
MCMCvis::MCMCsummary(fit, round = 3, 
                     params = c('beta'),
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

#combine data (scaled)
tt_comb2 <- rbind(tt_obs2, tt_imp2)

#scaling cov to measured scale bc transformed param est
#% change in LH trait for 1 sd change in covariate
beta1_rs_ch <- (exp(beta1_ch * sd(tt_comb2[,1] / lMass_scalar)) - 1) * 100
beta2_rs_ch <- (exp(beta2_ch * sd(tt_comb2[,2] / temp_sd_season_scalar)) - 1) * 100
beta3_rs_ch <- (exp(beta3_ch * sd(tt_comb2[,3] / temp_sd_year_scalar)) - 1) * 100
beta4_rs_ch <- (exp(beta4_ch * sd(tt_comb2[,4] / precip_cv_season_scalar)) - 1) * 100
beta5_rs_ch <- (exp(beta5_ch * sd(tt_comb2[,5] / precip_cv_year_scalar)) - 1) * 100


# added variable and partial resid plots ------------------------------------------------

# # https://www.wikiwand.com/en/Partial_residual_plot
# pr_fun <- function(num, nm)
# {
#   tm <- cbind(tt_comb2[,1],
#               c(DATA$temp_sd_season_obs, DATA$temp_sd_season_imp),
#               c(DATA$temp_sd_year_obs, DATA$temp_sd_year_imp),
#               c(DATA$precip_cv_season_obs, DATA$precip_cv_season_imp),
#               c(DATA$precip_cv_year_obs, DATA$precip_cv_year_imp))
# 
#   pch <- cbind(beta_ch,
#                gamma1_ch, gamma2_ch,
#                theta1_ch, theta2_ch)
# 
#   names <- c('Mass', 'Temp seasonality', 'Temp interannual',
#              'Precip seasonality', 'Precip interannual')
# 
#   #partial residuals
#   pr <- resid_comb + (median(pch[,num]) * tm[,num])
# 
#   pdf(paste0(fig_dir, nm, '-pr-', run_date, '.pdf'),
#       height = 5, width = 5)
#   plot(tm[,num], pr, col = rgb(0,0,0,0.2), pch = 19,
#        xlab = 'Predictor',
#        ylab = 'Partial residual',
#        main = names[num])
#   abline(h = 0, col = 'grey', lwd = 4, lty = 2)
#   abline(a = 0, b = median(pch[,num]), col = rgb(1,0,0,0.5), lwd = 4)
#   dev.off()
# }
# 
# pr_fun(num = 1, nm = 'mass') #Mass
# pr_fun(num = 2, nm = 'temp-season') #temp season
# pr_fun(num = 3, nm = 'temp-year') #temp year
# pr_fun(num = 4, nm = 'precip-season') #precip season
# pr_fun(num = 5, nm = 'precip-year') #precip year


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
                  main = '% change Surv for 1 sd change in cov',
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

# PPC - t
# mu_ch <- MCMCvis::MCMCchains(fit, params = 'mu')
# sigma_ch <- MCMCvis::MCMCchains(fit, params = 'sigma')
# nu_ch <- MCMCvis::MCMCchains(fit, params = 'nu')

# PPC - normal model
# mu_obs_ch <- MCMCvis::MCMCchains(fit, params = 'mu_obs')
# mu_imp_ch <- MCMCvis::MCMCchains(fit, params = 'mu_imp')
# sigma_ch <- MCMCvis::MCMCchains(fit, params = 'sigma')
# mu_comb_ch <- cbind(mu_obs_ch, mu_imp_ch)
# 
# #500 iterations
# sidx <- sample(1:NROW(mu_comb_ch), size = 500)
# y_rep <- matrix(NA, nrow = length(sidx), ncol = NCOL(mu_comb_ch))
# for (i in 1:length(sidx))
# {
#   #i <- 1
#   print(paste0('iter: ', i, ' of ', length(sidx)))
#   for (j in 1:NCOL(mu_comb_ch))
#   {
#     #t-dis
#     # eps <- rt(n = 1, df = nu_ch[sidx[i],1]) * sigma_ch[sidx[i],1]
#     # y_rep[i,j] <- mu_ch[sidx[i],j] + eps
#     y_rep[i,j] <- rnorm(1, mu_comb_ch[sidx[i],j], sigma_ch[sidx[i],1])
#   }
# }
# 
# pdf(paste0(fig_dir, 'PPC-', run_date, '.pdf'), height = 5, width = 5)
# plot(density(y_comb), col = 'black', lwd = 3, ylim = c(0, 5))#, xlim = c(0, 3.5))
# for (i in 1:500)
# {
#   lines(density(y_rep[i,]), col = rgb(1,0,0,0.2))
# }
# dev.off()


# R^2 ---------------------------------------------------------------------

#Gelman et al. 2019 - American Statistician
#var predicted / (var predicted + var residuals)
#in other words:
#explained var / (explained var + resid var)

#with mass - 0.16
# var_pred <- apply(mu_comb_ch, 1, var)
# var_resid <- apply(sweep(mu_comb_ch, 2, y_comb), 1, var)
# r2_ch <- var_pred / (var_pred + var_resid)
# hist(r2_ch)

# #no mass - ???
# mu_nm_ch <- MCMCvis::MCMCchains(fit, params = 'mu_nm')
# var_pred_nm <- apply(mu_nm_ch, 1, var)
# var_resid_nm <- apply(sweep(mu_nm_ch, 2, DATA$y), 1, var)
# r2_ch_nm <- var_pred_nm / (var_pred_nm + var_resid_nm)
# hist(r2_ch_nm)


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


