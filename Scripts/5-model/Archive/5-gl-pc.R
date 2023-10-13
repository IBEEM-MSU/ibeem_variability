####################
# Fit Bayes model - gen length ~ PC
####################


# specify dir -------------------------------------------------------------

dir <- '~/Google_Drive/Research/Projects/IBEEM_variabilty/'
run_date <- '2023-10-10'


# load packages -----------------------------------------------------------

library(tidyverse)
library(rstan)
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
                Mass, GenLength, lMass, lGL) %>%
  dplyr::mutate(species = stringr::str_to_title(gsub(' ', '_', Birdtree_name)))


# phylo -------------------------------------------------------------------

#load consensus tree - bird.phylo
load(paste0(dir, 'data/L3/bird-consensus-tree.rda'))

#df with names and idx
idx_df <- data.frame(idx = 1:NROW(tri),
                     name = bird_df2$species)

#species not found in both datasets (species to drop from tree)
nm <- setdiff(bird.phylo$tip.label, bird_df2$species)

#prune specified tips from all trees
pr_tree <- ape::drop.tip(bird.phylo, nm)


# PCA ---------------------------------------------------------------------

#for full run
j_idx3 <- dplyr::left_join(data.frame(species = pr_tree$tip.label),
                           data.frame(idx = 1:NROW(bird_df2), bird_df2),
                           by = 'species')

bird_df5 <- bird_df2[j_idx3$idx,]

#same results as raw, essentially
tt_pca <- dplyr::select(bird_df5, 
                        temp_sd_year, 
                        temp_sd_season,
                        precip_cv_year,
                        precip_cv_season) %>%
  prcomp(center = TRUE, scale. = TRUE)
factoextra::fviz_pca_var(tt_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')
factoextra::fviz_pca_var(tt_pca,
                         axes = c(2,3),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')

#add PC to data.frame
bird_df5$PC1 <- tt_pca$x[,1]
bird_df5$PC2 <- tt_pca$x[,2]
bird_df5$PC3 <- tt_pca$x[,3]


# scale/prep data ---------------------------------------------------------

#split predictors into obs and imp
tt <- data.frame(rep(1, NROW(bird_df5)),
                 bird_df5$lMass,
                 bird_df5$PC1,
                 bird_df5$PC2,
                 bird_df5$PC3)


# Run Stan model --------------------------------------------------------------

DATA <- list(N = NROW(bird_df5),
             y = bird_df5$lGL,
             K = NCOL(tt),
             X = tt,
             pro_data = bird_df5)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.93
TREE_DEPTH <- 12
STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 2000

#gen length - <1 min to fit
fit <- rstan::stan(paste0(dir, 'Scripts/Model_files/5-novar.stan'),
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('beta',
                            'sigma',
                            'mu'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))


# save summary space ------------------------------------------------------------

#save out summary, model fit, data
MCMCvis::MCMCdiag(fit, 
                  round = 4,
                  file_name = paste0('bird-gl-pc-results-', run_date),
                  dir = paste0(dir, 'Results'),
                  mkdir = paste0('bird-gl-pc-', run_date),
                  probs = c(0.055, 0.5, 0.945),
                  pg0 = TRUE,
                  save_obj = TRUE,
                  obj_name = paste0('bird-gl-pc-fit-', run_date),
                  add_obj = list(DATA),
                  add_obj_names = paste0('bird-gl-pc-data-', run_date),
                  cp_file = c(paste0(dir, 'Scripts/Model_files/5-novar.stan'), 
                              paste0(dir, 'Scripts/5-model/5-gl-pc.R')),
                  cp_file_names = c(paste0('5-novar-', run_date, '.stan'),
                                    paste0('5-gl-pc-', run_date, '.R')))

# library(shinystan)
# shinystan::launch_shinystan(fit)
# fit <- readRDS(paste0(dir, '/Results/se-bird-novar-oe-sep-', run_date,
#                       '/se-bird-novar-oe-sep-fit-', run_date, '.rds'))


# residuals ---------------------------------------------------------------

# extract residuals and calc phylo signal
mu_mn <- MCMCvis::MCMCpstr(fit, params = 'mu')[[1]]

#calc resid
resid <- DATA$y - mu_mn


# phylo signal in resids --------------------------------------------------

#df with names and idx
idx_df <- data.frame(idx = 1:length(resid), 
                     name = stringr::str_to_title(gsub(' ', '_', 
                                                       c(bird_df5$species))))

#get index for name order on tips
j_idx <- dplyr::left_join(data.frame(name = pr_tree$tip.label), idx_df, 
                          by = 'name')
resid_srt <- resid[j_idx$idx]

#K ~ 0.19
library(phytools)
phytools::phylosig(pr_tree, resid_srt, method = 'K') #quite slow
# phytools::phylosig(pr_tree, resid_srt, method = 'lambda')

# #just measured
# bird_ms <- dplyr::filter(bird_df3, !is.na(Measured_survival))
# 
# #prune specified tips from tree
# pr_tree2 <- ape::drop.tip(pr_tree, nm)
# 
# #get index for name order on tips
# j_idx <- dplyr::left_join(data.frame(name = pr_tree2$tip.label), idx_df, 
#                           by = 'name')
# 
# #apply
# bird_gls2 <- bird_ms[j_idx3$idx,]
# 
# library(nlme)
# pgls_fit2 <- nlme::gls(Measured_survival ~ lMass +
#                         temp_sd_season +
#                         temp_sd_year +
#                         precip_cv_season +
#                         precip_cv_year,
#                       # correlation = ape::corBrownian(phy = tree_n),
#                       correlation = ape::corPagel(1, pr_tree2, fixed = FALSE),
#                       data = bird_gls2,
#                       method = "ML")
# summary(pgls_fit2)
# car::vif(pgls_fit2)
#Pagel


# Summary -----------------------------------------------------------------

#model summary
MCMCvis::MCMCsummary(fit, round = 3, 
                     params = 'beta',
                     pg0 = TRUE)


# covariate effect on LH trait -----------------------------------------------

#INTERPRETATION
#((e^param) - 1) * 100 = percent change in trait for every one unit change in covariate
#((e^(param * L)) - 1) * 100 = percent change in trait for every L unit change in covariate
beta1_ch <- MCMCvis::MCMCchains(fit, params = 'beta[1]', ISB = FALSE, exact = TRUE)
beta2_ch <- MCMCvis::MCMCchains(fit, params = 'beta[2]', ISB = FALSE, exact = TRUE)
beta3_ch <- MCMCvis::MCMCchains(fit, params = 'beta[3]', ISB = FALSE, exact = TRUE)
beta4_ch <- MCMCvis::MCMCchains(fit, params = 'beta[4]', ISB = FALSE, exact = TRUE)
beta5_ch <- MCMCvis::MCMCchains(fit, params = 'beta[5]', ISB = FALSE, exact = TRUE)

# median((exp(beta_ch * diff(range(DATA$lMass))) - 1) * 100)
# median((exp(gamma1_ch * diff(range(DATA$temp_sd_season))) - 1) * 100)
# median((exp(gamma2_ch * diff(range(DATA$temp_sd_year))) - 1) * 100)
# median((exp(theta1_ch * diff(range(DATA$precip_cv_season))) - 1) * 100)
# median((exp(theta2_ch * diff(range(DATA$precip_cv_year))) - 1) * 100)

#% change in LH trait for 1 sd change in covariate
beta2_rs_ch <- (exp(beta2_ch * sd(tt[,2])) - 1) * 100
beta3_rs_ch <- (exp(beta3_ch * sd(tt[,3])) - 1) * 100
beta4_rs_ch <- (exp(beta4_ch * sd(tt[,4])) - 1) * 100
beta5_rs_ch <- (exp(beta5_ch * sd(tt[,5])) - 1) * 100


# added variable and partial resid plots ------------------------------------------------

fig_dir <- paste0(dir, 'Results/bird-gl-pc-', run_date, '/')

# https://www.wikiwand.com/en/Partial_residual_plot
pr_fun <- function(num, nm)
{
  tm <- tt[,-1]
  
  pch <- cbind(beta2_ch,
               beta3_ch, beta4_ch,
               beta5_ch)
  
  names <- c('Mass', 'PC1', 'PC2',
             'PC3')
  
  #partial residuals
  pr <- resid + (median(pch[,num]) * tm[,num])
  
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

pr_fun(num = 1, nm = 'Mass') #Mass
pr_fun(num = 2, nm = 'PC1')
pr_fun(num = 3, nm = 'PC2')
pr_fun(num = 4, nm = 'PC3')


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
                        beta4_rs_ch,
                        beta5_rs_ch),
                  labels = c('PC1',
                             'PC2',
                             'PC3'),
                  sz_labels = 1.5,
                  ci = c(89, 89),
                  sz_thick = 3,
                  sz_thin = 3,
                  main = '% change GL for 1 sd change in cov',
                  guide_lines = TRUE)
dev.off()


# PCA plots ---------------------------------------------------------------

pdf(paste0(fig_dir, 'PC-plt-1-', run_date, '.pdf'),
    height = 5, width = 5)
factoextra::fviz_pca_var(tt_pca,
                         axes = c(1,2),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')
dev.off()

pdf(paste0(fig_dir, 'PC-plt-2-', run_date, '.pdf'),
    height = 5, width = 5)
factoextra::fviz_pca_var(tt_pca,
                         axes = c(2,3),
                         #geom = 'arrow',
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         #repel = FALSE,
                         title = 'PCA')
dev.off()


# PPC ---------------------------------------------------------------------

# PPC - t
# mu_ch <- MCMCvis::MCMCchains(fit, params = 'mu')
# sigma_ch <- MCMCvis::MCMCchains(fit, params = 'sigma')
# nu_ch <- MCMCvis::MCMCchains(fit, params = 'nu')

# # PPC - normal model
mu_ch <- MCMCvis::MCMCchains(fit, params = 'mu')
sigma_ch <- MCMCvis::MCMCchains(fit, params = 'sigma')

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

#with mass
var_pred <- apply(mu_ch, 1, var)
var_resid <- apply(sweep(mu_ch, 2, DATA$y), 1, var)
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
          data = bird_df5)
stf1 <- summary(tf1)

tf2 <- lm(temp_sd_year ~ lMass + temp_sd_season +
            precip_cv_year + precip_cv_season, data = 
            bird_df5)
stf2 <- summary(tf2)

tf3 <- lm(temp_sd_season ~ lMass + temp_sd_year +
            precip_cv_year + precip_cv_season, 
          data = bird_df5)
stf3 <- summary(tf3)

tf4 <- lm(precip_cv_year ~ lMass + temp_sd_year + temp_sd_season + 
            precip_cv_season, 
          data = bird_df5)
stf4 <- summary(tf4)

tf5 <- lm(precip_cv_season ~ lMass + temp_sd_year + temp_sd_season + 
            precip_cv_year, 
          data = bird_df5)
stf5 <- summary(tf5)


#calc VIF per covariate
1 / (1 - stf1$r.squared) #lMass
1 / (1 - stf2$r.squared) #temp_sd_year
1 / (1 - stf3$r.squared) #temp_sd_season
1 / (1 - stf4$r.squared) #precip_cv_year
1 / (1 - stf5$r.squared) #precip_cv_season

#correlation
cor(as.matrix(dplyr::select(bird_df5,
                            lMass, temp_sd_year, temp_sd_season,
                            precip_cv_year, precip_cv_season)))
