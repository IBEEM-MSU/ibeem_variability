####################
# Fit Bayes model - gen length ~ env (varying intercepts or slopes by family)
# 
####################


# specify dir -------------------------------------------------------------

dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
run_date <- '2023-09-22'


# load packages -----------------------------------------------------------

library(tidyverse)
library(rstan)
library(MCMCvis)
library(ape)
library(picante)


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

#just one Fam
of <- dplyr::filter(bird_df, Family == 'Psittacidae')
summary(lm(lGL ~ lMass + temp_sp_color_month + temp_sd_season + temp_sd_year, 
           data = of))


# # load mammal data --------------------------------------------------------
# 
# mam_df <- read.csv(paste0(dir, 'Data/L3/main-mammal-data.csv')) %>%
#   dplyr::mutate(Family = PH_Family,
#                 Order = PH_Order,
#                 # LH_Mass = LH_AdultBodyMass_g, #Pacifici mass
#                 Mass = PH_Mass.g, #Phylacine mass
#                 GenLength = LH_GenerationLength_d,
#                 fac_Family = factor(Family),
#                 fac_Order = factor(Order),
#                 lMass = log(Mass),
#                 lGL = log(GenLength))
# 
# #one species with Inf precip_cv_season (due to precip_mean = 0 and using median)
# dplyr::filter(mam_df, !is.finite(precip_cv_season))
# mam_df$precip_cv_season[which(!is.finite(mam_df$precip_cv_season))] <- 0


# Run Stan model --------------------------------------------------------------

#data for model
DATA <- list(N = NROW(bird_df),
             Nf = length(unique(bird_df$f_id)),
             y = bird_df$lGL,
             f_id = bird_df$f_id, #family id for each data point
             lMass = bird_df$lMass,
             temp_sd_season = bird_df$temp_sd_season,
             temp_sd_year = bird_df$temp_sd_year,
             temp_sp_color_month = bird_df$temp_sp_color_month,
             precip_cv_season = bird_df$precip_cv_season,
             precip_cv_year = bird_df$precip_cv_year,
             precip_sp_color_month = bird_df$precip_sp_color_month,
             pro_data = bird_df)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.95
TREE_DEPTH <- 12
STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 2000

#varying intercepts or slopes
fit <- rstan::stan(paste0(dir, 'Scripts/Model_files/5-varfam.stan'),
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'gamma1',
                            'gamma2',
                            'gamma3',
                            'theta1',
                            'theta2',
                            'theta3',
                            'sigma',
                            'mu_alpha',
                            'mu_beta',
                            'mu_gamma1',
                            'mu_gamma2',
                            'mu_gamma3',
                            'mu_theta1',
                            'mu_theta2',
                            'mu_theta3',
                            'nu',
                            'sigma_abgt',
                            'Rho',
                            'mu',
                            'mu_nm'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fit, 
                  round = 4,
                  file_name = paste0('ge-bird-varfam-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('ge-bird-varfam-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('ge-bird-varfam-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('ge-bird-varfam-data-', run_date),
                  cp_file = c(paste0(dir, 'Scripts/Model_files/5-varfam.stan'), 
                              paste0(dir, 'Scripts/5-model/5-varfam.R')),
                  cp_file_names = c(paste0('5-varfam-', run_date, '.stan'),
                                    paste0('5-varfam-', run_date, '.R')))

# library(shinystan)
# shinystan::launch_shinystan(fit)


# phylo signal in resids --------------------------------------------------

#NOTE: picante::phylosignal is quite slow. At least 15min for one tree...
#Seems to be some memory pressure
#https://cran.r-project.org/web/packages/picante/vignettes/picante-intro.pdf
#no strong phylo signal detected in tree 1 (K ~ 0.05). "K values closer to zero correspond to a random or convergent pattern of evolution, while K values greater than 1 indicate strong phylogenetic signal and conservatism of traits" (picante vignette)

#extract resids
mu_mn <- MCMCvis::MCMCpstr(fit, params = 'mu')[[1]]
resids <- DATA$y - mu_mn

#read in phylo tree (birdtree.org - Ericson 0001-1000)
tr <- ape::read.tree(paste0(dir, 'data/L1/trait/AllBirdsEricson1.tre'))

#species not found in both datasets (species to drop from tree)
nm <- setdiff(tr[[1]]$tip.label, idx_df$name)

#prune specified tips from all trees
pr_tr <- lapply(tr, drop.tip, tip = nm)
class(pr_tr) <- "multiPhylo"

#df with names and idx
idx_df <- data.frame(idx = 1:NROW(bird_df), 
                     name = stringr::str_to_title(gsub(' ', '_', bird_df$Birdtree_name)))

#for each of 1000 trees, calculate phylo signal (Blomberg's K) in resids
out.df <- data.frame(K = rep(NA, length(pr_tr)), 
                     PIC.var.P = NA)
for (i in 1:100)
{
  #i <- 2
  print(paste0('tree: ', i, ' of ', length(pr_tr)))
  tree_n <- pr_tr[[i]]
  
  #get index for name order on tips
  j_idx <- dplyr::left_join(data.frame(name = tree_n$tip.label), idx_df, 
                            by = 'name')
  #apply to residuals
  resid_srt <- resids[j_idx$idx]
  phy_res <- picante::phylosignal(resid_srt, tree_n) #quite slow
  out.df$K[i] <- phy_res$K
  out.df$PIC.var.P[i] <- phy_res$PIC.variance.P
}

#summarize output
hist(out.df$PIC.var.P)
sum(out.df$PIC.var.P < 0.05)

mean(out.df$K)
sd(out.df$K)


# Summary -----------------------------------------------------------------

#model summary
MCMCvis::MCMCsummary(fit, round = 3, 
                     params = c('mu_beta', #effect log mass
                                'mu_gamma1', #effect temp seasonality
                                'mu_gamma2', #effect temp interannual var
                                'mu_gamma3', #effect temp spectra
                                'mu_theta1', #effect precip seasonality
                                'mu_theta2', #effect precip interranvual var
                                'mu_theta3'), #effect precip spectra
                     pg0 = TRUE)


# covariate effect on gen length -----------------------------------------------

#INTERPRETATION
#((e^param) - 1) * 100 = percent change in trait for every one unit change in covariate
#((e^(param * L)) - 1) * 100 = percent change in trait for every L unit change in covariate
mu_beta_ch <- MCMCvis::MCMCchains(fit, params = 'mu_beta')
mu_gamma1_ch <- MCMCvis::MCMCchains(fit, params = 'mu_gamma1')
mu_gamma2_ch <- MCMCvis::MCMCchains(fit, params = 'mu_gamma2')
mu_gamma3_ch <- MCMCvis::MCMCchains(fit, params = 'mu_gamma3')
mu_theta1_ch <- MCMCvis::MCMCchains(fit, params = 'mu_theta1')
mu_theta2_ch <- MCMCvis::MCMCchains(fit, params = 'mu_theta2')
mu_theta3_ch <- MCMCvis::MCMCchains(fit, params = 'mu_theta3')

# median((exp(mu_beta_ch * diff(range(DATA$lMass))) - 1) * 100)
# median((exp(mu_gamma1_ch * diff(range(DATA$temp_sd_season))) - 1) * 100)
# median((exp(mu_gamma2_ch * diff(range(DATA$temp_sd_year))) - 1) * 100)
# median((exp(mu_gamma3_ch * diff(range(DATA$temp_sp_color_month))) - 1) * 100)
# median((exp(mu_theta1_ch * diff(range(DATA$precip_cv_season))) - 1) * 100)
# median((exp(mu_theta2_ch * diff(range(DATA$precip_cv_year))) - 1) * 100)
# median((exp(mu_theta3_ch * diff(range(DATA$precip_sp_color_month))) - 1) * 100)

#% change in gen length for 1 sd change in covariate
mu_beta_rs_ch <- (exp(mu_beta_ch * sd(DATA$lMass)) - 1) * 100
mu_gamma1_rs_ch <- (exp(mu_gamma1_ch * sd(DATA$temp_sd_season)) - 1) * 100
mu_gamma2_rs_ch <- (exp(mu_gamma2_ch * sd(DATA$temp_sd_year)) - 1) * 100
mu_gamma3_rs_ch <- (exp(mu_gamma3_ch * sd(DATA$temp_sp_color_month)) - 1) * 100
mu_theta1_rs_ch <- (exp(mu_theta1_ch * sd(DATA$precip_cv_season)) - 1) * 100
mu_theta2_rs_ch <- (exp(mu_theta2_ch * sd(DATA$precip_cv_year)) - 1) * 100
mu_theta3_rs_ch <- (exp(mu_theta3_ch * sd(DATA$precip_sp_color_month)) - 1) * 100


# cat plots ---------------------------------------------------------------

fig_dir <- paste0(dir, 'Results/ge-bird-varfam-', run_date, '/')

#temp
pdf(paste0(fig_dir, 'param-cat-raw-', run_date, '.pdf'),
    height = 14, width = 5)
MCMCvis::MCMCplot(fit,
                  params = c('mu_beta', 
                             'mu_gamma1', 'mu_gamma2', 'mu_gamma3',
                             'mu_theta1', 'mu_theta2', 'mu_theta3'),
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
MCMCvis::MCMCplot(cbind(mu_beta_rs_ch, 
                        mu_gamma1_rs_ch,
                        mu_gamma2_rs_ch,
                        mu_gamma3_rs_ch,
                        mu_theta1_rs_ch,
                        mu_theta2_rs_ch,
                        mu_theta3_rs_ch),
                  labels = c('log(Mass)',
                             'T season',
                             'T interannual var',
                             'T spectra',
                             'P seasonality',
                             'P interannual var',
                             'P spectra'),
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = '% change in Gen Length for 1 sd change in cov',
                  guide_lines = TRUE)
dev.off()


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

pdf(paste0(fig_dir, 'varfam_PPC-', run_date, '.pdf'), height = 5, width = 5)
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

#with mass
var_pred <- apply(mu_ch, 1, var)
var_resid <- apply(sweep(mu_ch, 2, DATA$y), 1, var)
r2_ch <- var_pred / (var_pred + var_resid)
hist(r2_ch)

#no mass
mu_nm_ch <- MCMCvis::MCMCchains(fit, params = 'mu_nm')
var_pred_nm <- apply(mu_nm_ch, 1, var)
var_resid_nm <- apply(sweep(mu_nm_ch, 2, DATA$y), 1, var)
r2_ch_nm <- var_pred_nm / (var_pred_nm + var_resid_nm)
hist(r2_ch_nm)


# VIF ---------------------------------------------------------------------

#covariates as a function of other covariates
tf1 <- lm(lMass ~ temp_sd_year + temp_sd_season + temp_sp_color_month +
            precip_cv_year + precip_cv_season + precip_sp_color_month, 
          data = bird_df)
stf1 <- summary(tf1)

tf2 <- lm(temp_sd_year ~ lMass + temp_sd_season + temp_sp_color_month +
            precip_cv_year + precip_cv_season + precip_sp_color_month, data = 
            bird_df)
stf2 <- summary(tf2)

tf3 <- lm(temp_sd_season ~ lMass + temp_sd_year + temp_sp_color_month +
            precip_cv_year + precip_cv_season + precip_sp_color_month, 
          data = bird_df)
stf3 <- summary(tf3)

tf4 <- lm(temp_sp_color_month ~ lMass + temp_sd_year + temp_sd_season +
            precip_cv_year + precip_cv_season + precip_sp_color_month, 
          data = bird_df)
stf4 <- summary(tf4)

tf5 <- lm(precip_cv_year ~ lMass + temp_sd_year + temp_sd_season + 
            temp_sp_color_month + precip_cv_season + precip_sp_color_month, 
          data = bird_df)
stf5 <- summary(tf5)

tf6 <- lm(precip_cv_season ~ lMass + temp_sd_year + temp_sd_season + 
            temp_sp_color_month + precip_cv_year + precip_sp_color_month, 
          data = bird_df)
stf6 <- summary(tf6)

tf7 <- lm(precip_sp_color_month ~ lMass + temp_sd_year + temp_sd_season + 
            temp_sp_color_month + precip_cv_year + precip_cv_season, 
          data = bird_df)
stf7 <- summary(tf7)

#calc VIF per covariate
1 / (1 - stf1$r.squared) #lMass
1 / (1 - stf2$r.squared) #temp_sd_year
1 / (1 - stf3$r.squared) #temp_sd_season
1 / (1 - stf4$r.squared) #temp_sp_color_month
1 / (1 - stf5$r.squared) #precip_cv_year
1 / (1 - stf6$r.squared) #precip_cv_season
1 / (1 - stf7$r.squared) #precip_sp_color_month

#correlation
cor(as.matrix(dplyr::select(bird_df,
                            lMass, temp_sd_year, temp_sd_season, temp_sp_color_month,
                            precip_cv_year, precip_cv_season, precip_sp_color_month)))
