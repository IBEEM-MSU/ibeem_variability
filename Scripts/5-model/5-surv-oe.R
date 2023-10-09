####################
# Fit Bayes model - survival ~ env
# Estimate surv from traits and phylo
# Fit model incorporating obs error
####################


# specify dir -------------------------------------------------------------

dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
run_date <- '2023-10-09'


# load packages -----------------------------------------------------------

library(tidyverse)
library(rstan)
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

#how many species in each family
dplyr::group_by(bird_df, Family) %>%
  dplyr::count() %>%
  dplyr::arrange(desc(n))

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

#run phylo imputation
ir <- Rphylopars::phylopars(trait_data = tri[j_idx$idx,], 
                            tree = pr_tree, 
                            phylo_correlated = TRUE,
                            # model = "BM") # AIC = 20804
                            # model = "OU") # AIC = 1942203
                            model = "lambda") # AIC = 10755
# model = "mvOU") # threw an error
# model = "delta") # threw an error
# model = "EB") # AIC = 20697
# model = "star") # ???

# AIC(ir)

#returns species means as well as some internal nodes
#get just species rows
ridx <- which(row.names(ir$anc_rec) %in% tri$species)

#variance (uncertainty)
ir_unc <- data.frame(species = tri$species[j_idx$idx], 
                 ir$anc_var[ridx,]) %>%
  dplyr::mutate(SD_survival = sqrt(Measured_survival),
                SD_log_age_first_breeding = sqrt(Measured_log_age_first_breeding),
                SD_log_max_longevity = sqrt(Measured_log_max_longevity),
                SD_log_clutch_size = sqrt(Measured_log_clutch_size)) %>%
  dplyr::arrange(species) %>%
  dplyr::select(species, 
                SD_survival, 
                SD_log_age_first_breeding,
                SD_log_max_longevity,
                SD_log_clutch_size)

#merge imputed values with unc and env data
ir_mrg <- data.frame(species = tri$species[j_idx$idx], 
                  ir$anc_rec[ridx,]) %>%
  dplyr::arrange(species) %>%
  dplyr::rename(Phylo_survival = Measured_survival,
                Phylo_log_age_first_breeding = Measured_log_age_first_breeding,
                Phylo_log_max_longevity = Measured_log_max_longevity,
                Phylo_log_clutch_size = Measured_log_clutch_size) %>%
  #join with measured traits
  dplyr::left_join(dplyr::select(tri, -lMass), by = 'species') %>%
  #join with var
  dplyr::left_join(ir_unc, by = 'species')
row.names(ir_mrg) <- NULL

#join trait and env data
bird_df3 <- dplyr::mutate(bird_df2,
                      species = stringr::str_to_title(gsub(' ', '_', Birdtree_name))) %>%
  dplyr::select(species,
                temp_mean,
                temp_sd_year,
                temp_sd_season,
                temp_sp_color_month,
                precip_cv_year,
                precip_cv_season,
                precip_sp_color_month,
                Modeled_survival,
                Modeled_age_first_breeding,
                Modeled_max_longevity) %>%
  dplyr::left_join(ir_mrg, by = 'species')

# saveRDS(bird_df3, paste0(dir, 'Scripts/5-model/bird_df3.rds'))


# Run Stan model --------------------------------------------------------------

#idx observed data
obs_idx <- which(bird_df3$SD_survival == 0)
imp_idx <- which(bird_df3$SD_survival != 0)

DATA <- list(No = length(obs_idx),
             Ni = length(imp_idx),
             y_obs = bird_df3$Phylo_survival[obs_idx],
             y_imp = bird_df3$Phylo_survival[imp_idx],
             sd_y = bird_df3$SD_survival[imp_idx],
             lMass_obs = bird_df3$lMass[obs_idx],
             temp_sd_season_obs = bird_df3$temp_sd_season[obs_idx],
             temp_sd_year_obs = bird_df3$temp_sd_year[obs_idx],
             precip_cv_season_obs = bird_df3$precip_cv_season[obs_idx],
             precip_cv_year_obs = bird_df3$precip_cv_year[obs_idx],
             lMass_imp = bird_df3$lMass[imp_idx],
             temp_sd_season_imp = bird_df3$temp_sd_season[imp_idx],
             temp_sd_year_imp = bird_df3$temp_sd_year[imp_idx],
             precip_cv_season_imp = bird_df3$precip_cv_season[imp_idx],
             precip_cv_year_imp = bird_df3$precip_cv_year[imp_idx],
             pro_data = bird_df3)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.95
TREE_DEPTH <- 12
STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 2000

#survival - ~25 min to fit
#obs err model
fit <- rstan::stan(paste0(dir, 'Scripts/Model_files/5-oe.stan'),
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'gamma1',
                            'gamma2',
                            # 'gamma3',
                            'theta1',
                            'theta2',
                            # 'theta3',
                            'sigma',
                            # 'nu',
                            'mu_obs',
                            'mu_imp',
                            'y_iv'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fit, 
                  round = 4,
                  file_name = paste0('bird-surv-oe-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('bird-surv-oe-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('bird-surv-oe-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('bird-surv-oe-data-', run_date),
                  cp_file = c(paste0(dir, 'Scripts/Model_files/5-oe.stan'), 
                              paste0(dir, 'Scripts/5-model/5-surv-oe.R')),
                  cp_file_names = c(paste0('5-oe-', run_date, '.stan'),
                                    paste0('5-surv-oe-', run_date, '.R')))

# library(shinystan)
# shinystan::launch_shinystan(fit)
# fit <- readRDS(paste0(dir, '/Results/se-bird-novar-oe-sep-', run_date,
#                       '/se-bird-novar-oe-sep-fit-', run_date, '.rds'))


# residuals ---------------------------------------------------------------

# extract residuals and calc phylo signal
mu_obs_mn <- MCMCvis::MCMCpstr(fit, params = 'mu_obs')[[1]]
mu_imp_mn <- MCMCvis::MCMCpstr(fit, params = 'mu_imp')[[1]]

#combine and calc resid
mu_comb <- c(mu_obs_mn, mu_imp_mn)
y_comb <- c(DATA$y_obs, DATA$y_imp)
resid_comb <- y_comb - mu_comb


# phylo signal in resids --------------------------------------------------

#NOTE: picante::phylosignal is quite slow. At least 30min for one tree...
#Seems to be some memory pressure
#https://cran.r-project.org/web/packages/picante/vignettes/picante-intro.pdf
#no strong phylo signal detected in tree 1 (K ~ 0.05). "K values closer to zero correspond to a random or convergent pattern of evolution, while K values greater than 1 indicate strong phylogenetic signal and conservatism of traits" (picante vignette)

#df with names and idx
idx_df <- data.frame(idx = 1:length(y_comb), 
                     name = stringr::str_to_title(gsub(' ', '_', 
                                                       c(bird_df3$species[obs_idx],
                                                         bird_df3$species[imp_idx]))))

#get index for name order on tips
j_idx <- dplyr::left_join(data.frame(name = tree_n$tip.label), idx_df, 
                          by = 'name')
resid_srt <- resid_comb[j_idx$idx]

#apply to residuals
# phy_res <- picante::phylosignal(resid_srt, tree_n) #quite slow (10s of minutes)

#K ~ 0.3
#lambda ~ 
# library(phytools)
# phytools::phylosig(tree_n, resid_srt, method = 'K') #quit slow
# phytools::phylosig(tree_n, resid_srt, method = 'lambda') #wouldn't finish??? (many hours)


#fit alternative model - no measures err - using nlme
j_idx2 <- dplyr::left_join(data.frame(species = tree_n$tip.label), 
                           data.frame(idx = 1:NROW(bird_df3), bird_df3), 
                          by = 'species')
#apply
bird_gls <- bird_df3[j_idx2$idx,]
library(nlme)
pgls_fit <- nlme::gls(Phylo_survival ~ lMass +
                         temp_sd_season +
                         temp_sd_year +
                         precip_cv_season +
                         precip_cv_year,
                       correlation = ape::corPagel(1, tree_n, fixed = FALSE),
                       data = bird_gls,
                       method = "ML")
summary(pgls_fit)
# car::vif(pgls_fit)
#Pagel (1)
# Coefficients:
#   Value  Std.Error   t-value p-value
# (Intercept)       0.3577705 0.03975689   8.99896  0.0000
# lMass             0.0458616 0.00084637  54.18605  0.0000
# temp_sd_season   -0.0028178 0.00026847 -10.49553  0.0000
# temp_sd_year     -0.0041413 0.00403653  -1.02596  0.3049
# precip_cv_season -0.0033795 0.00232325  -1.45464  0.1458
# precip_cv_year    0.0142654 0.01024590   1.39230  0.1639

#just measured
bird_ms <- dplyr::filter(bird_df3, !is.na(Measured_survival))

#prune tree
#df with names and idx
idx_df2 <- data.frame(idx = 1:NROW(bird_ms), 
                     name = stringr::str_to_title(gsub(' ', '_', bird_ms$species)))

#species not found in both datasets (species to drop from tree)
nm2 <- setdiff(tree[[1]]$tip.label, bird_ms$species)

#prune specified tips from all trees
pr_tree2 <- lapply(tree, ape::drop.tip, tip = nm2)
class(pr_tree2) <- "multiPhylo"
tree_n2 <- pr_tree2[[1]]

#get idx
j_idx3 <- dplyr::left_join(data.frame(species = tree_n2$tip.label), 
                           data.frame(idx = 1:NROW(bird_ms), bird_ms), 
                           by = 'species')
#apply
bird_gls2 <- bird_ms[j_idx3$idx,]

library(nlme)
pgls_fit2 <- nlme::gls(Measured_survival ~ lMass +
                        temp_sd_season +
                        temp_sd_year +
                        precip_cv_season +
                        precip_cv_year,
                      # correlation = ape::corBrownian(phy = tree_n),
                      correlation = ape::corPagel(1, tree_n2, fixed = FALSE),
                      data = bird_gls2,
                      method = "ML")
summary(pgls_fit2)
car::vif(pgls_fit2)
#Pagel
# Coefficients:
#   Value  Std.Error   t-value p-value
# (Intercept)       0.3459467 0.08343573  4.146266  0.0000
# lMass             0.0499811 0.01056317  4.731644  0.0000
# temp_sd_season   -0.0115342 0.00558003 -2.067054  0.0395
# temp_sd_year      0.0268342 0.11241383  0.238709  0.8115
# precip_cv_season  0.0246000 0.06265611  0.392619  0.6949
# precip_cv_year    0.2150230 0.23048542  0.932914  0.3515


# Summary -----------------------------------------------------------------

#model summary
MCMCvis::MCMCsummary(fit, round = 3, 
                     params = c('beta', #effect log mass
                                'gamma1', #effect temp seasonality
                                'gamma2', #effect temp interannual var
                                'theta1', #effect precip seasonality
                                'theta2'), #effect precip interranvual var
                     pg0 = TRUE)


# covariate effect on LH trait -----------------------------------------------

#INTERPRETATION
#((e^param) - 1) * 100 = percent change in trait for every one unit change in covariate
#((e^(param * L)) - 1) * 100 = percent change in trait for every L unit change in covariate
beta_ch <- MCMCvis::MCMCchains(fit, params = 'beta')
gamma1_ch <- MCMCvis::MCMCchains(fit, params = 'gamma1')
gamma2_ch <- MCMCvis::MCMCchains(fit, params = 'gamma2')
theta1_ch <- MCMCvis::MCMCchains(fit, params = 'theta1')
theta2_ch <- MCMCvis::MCMCchains(fit, params = 'theta2')

# median((exp(beta_ch * diff(range(DATA$lMass))) - 1) * 100)
# median((exp(gamma1_ch * diff(range(DATA$temp_sd_season))) - 1) * 100)
# median((exp(gamma2_ch * diff(range(DATA$temp_sd_year))) - 1) * 100)
# median((exp(theta1_ch * diff(range(DATA$precip_cv_season))) - 1) * 100)
# median((exp(theta2_ch * diff(range(DATA$precip_cv_year))) - 1) * 100)

#% change in LH trait for 1 sd change in covariate
beta_rs_ch <- (exp(beta_ch * sd(c(DATA$lMass_obs, 
                                  DATA$lMass_imp))) - 1) * 100
gamma1_rs_ch <- (exp(gamma1_ch * sd(c(DATA$temp_sd_season_obs,
                                    DATA$temp_sd_season_imp))) - 1) * 100
gamma2_rs_ch <- (exp(gamma2_ch * sd(c(DATA$temp_sd_year_obs,
                                    DATA$temp_sd_year_imp))) - 1) * 100
theta1_rs_ch <- (exp(theta1_ch * sd(c(DATA$precip_cv_season_obs,
                                    DATA$precip_cv_season_imp))) - 1) * 100
theta2_rs_ch <- (exp(theta2_ch * sd(c(DATA$precip_cv_year_obs,
                                      DATA$precip_cv_year_imp))) - 1) * 100


# added variable and partial resid plots ------------------------------------------------

fig_dir <- paste0('bird-surv-oe-', run_date, '/')

# https://www.wikiwand.com/en/Partial_residual_plot
pr_fun <- function(num, nm)
{
  tm <- cbind(c(DATA$lMass_obs, DATA$lMass_imp),
              c(DATA$temp_sd_season_obs, DATA$temp_sd_season_imp),
              c(DATA$temp_sd_year_obs, DATA$temp_sd_year_imp),
              c(DATA$precip_cv_season_obs, DATA$precip_cv_season_imp),
              c(DATA$precip_cv_year_obs, DATA$precip_cv_year_imp))
  
  pch <- cbind(beta_ch, 
               gamma1_ch, gamma2_ch,
               theta1_ch, theta2_ch)
  
  names <- c('Mass', 'Temp seasonality', 'Temp interannual',
             'Precip seasonality', 'Precip interannual')
  
  #partial residuals
  pr <- resid_comb + (median(pch[,num]) * tm[,num])
  
  pdf(paste0(fig_dir, nm, '-pr-', run_date, '.pdf'),
      height = 5, width = 5)
  plot(tm[,num], pr, col = rgb(0,0,0,0.2), pch = 19,
       xlab = 'Predictor',
       ylab = 'Partial residual',
       main = names[num])
  abline(h = 0, col = 'grey', lwd = 4, lty = 2)
  abline(a = 0, b = median(pch[,num]), col = rgb(1,0,0,0.5), lwd = 4)
  dev.off()
}

pr_fun(num = 1, nm = 'mass') #Mass
pr_fun(num = 2, nm = 'temp-season') #temp season
pr_fun(num = 3, nm = 'temp-year') #temp year
pr_fun(num = 4, nm = 'precip-season') #precip season
pr_fun(num = 5, nm = 'precip-year') #precip year


# cat plots ---------------------------------------------------------------

#temp
pdf(paste0(fig_dir, 'param-cat-raw-', run_date, '.pdf'),
    height = 5, width = 5)
MCMCvis::MCMCplot(fit,
                  params = c('beta', 
                             'gamma1', 'gamma2',
                             'theta1', 'theta2'),
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
MCMCvis::MCMCplot(cbind(gamma1_rs_ch,
                        gamma2_rs_ch,
                        theta1_rs_ch,
                        theta2_rs_ch),
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


# PPC ---------------------------------------------------------------------

# PPC - t
# mu_ch <- MCMCvis::MCMCchains(fit, params = 'mu')
# sigma_ch <- MCMCvis::MCMCchains(fit, params = 'sigma')
# nu_ch <- MCMCvis::MCMCchains(fit, params = 'nu')

# PPC - normal model
mu_obs_ch <- MCMCvis::MCMCchains(fit, params = 'mu_obs')
mu_imp_ch <- MCMCvis::MCMCchains(fit, params = 'mu_imp')
sigma_ch <- MCMCvis::MCMCchains(fit, params = 'sigma')
mu_comb_ch <- cbind(mu_obs_ch, mu_imp_ch)

#500 iterations
sidx <- sample(1:NROW(mu_comb_ch), size = 500)
y_rep <- matrix(NA, nrow = length(sidx), ncol = NCOL(mu_comb_ch))
for (i in 1:length(sidx))
{
  #i <- 1
  print(paste0('iter: ', i, ' of ', length(sidx)))
  for (j in 1:NCOL(mu_comb_ch))
  {
    #t-dis
    # eps <- rt(n = 1, df = nu_ch[sidx[i],1]) * sigma_ch[sidx[i],1]
    # y_rep[i,j] <- mu_ch[sidx[i],j] + eps
    y_rep[i,j] <- rnorm(1, mu_comb_ch[sidx[i],j], sigma_ch[sidx[i],1])
  }
}

pdf(paste0(fig_dir, 'PPC-', run_date, '.pdf'), height = 5, width = 5)
plot(density(y_comb), col = 'black', lwd = 3, ylim = c(0, 5))#, xlim = c(0, 3.5))
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

#with mass - 0.16
var_pred <- apply(mu_comb_ch, 1, var)
var_resid <- apply(sweep(mu_comb_ch, 2, y_comb), 1, var)
r2_ch <- var_pred / (var_pred + var_resid)
hist(r2_ch)

# #no mass - ???
# mu_nm_ch <- MCMCvis::MCMCchains(fit, params = 'mu_nm')
# var_pred_nm <- apply(mu_nm_ch, 1, var)
# var_resid_nm <- apply(sweep(mu_nm_ch, 2, DATA$y), 1, var)
# r2_ch_nm <- var_pred_nm / (var_pred_nm + var_resid_nm)
# hist(r2_ch_nm)


# VIF ---------------------------------------------------------------------

#covariates as a function of other covariates
tf1 <- lm(lMass ~ temp_sd_year + temp_sd_season +
            precip_cv_year + precip_cv_season, 
          data = bird_df)
stf1 <- summary(tf1)

tf2 <- lm(temp_sd_year ~ lMass + temp_sd_season +
            precip_cv_year + precip_cv_season, data = 
            bird_df)
stf2 <- summary(tf2)

tf3 <- lm(temp_sd_season ~ lMass + temp_sd_year +
            precip_cv_year + precip_cv_season, 
          data = bird_df)
stf3 <- summary(tf3)

tf4 <- lm(precip_cv_year ~ lMass + temp_sd_year + temp_sd_season + 
             precip_cv_season, 
          data = bird_df)
stf4 <- summary(tf4)

tf5 <- lm(precip_cv_season ~ lMass + temp_sd_year + temp_sd_season + 
            precip_cv_year, 
          data = bird_df)
stf5 <- summary(tf5)


#calc VIF per covariate
1 / (1 - stf1$r.squared) #lMass
1 / (1 - stf2$r.squared) #temp_sd_year
1 / (1 - stf3$r.squared) #temp_sd_season
1 / (1 - stf4$r.squared) #precip_cv_year
1 / (1 - stf5$r.squared) #precip_cv_season

#correlation
cor(as.matrix(dplyr::select(bird_df3,
                            lMass, temp_sd_year, temp_sd_season,
                            precip_cv_year, precip_cv_season)))
