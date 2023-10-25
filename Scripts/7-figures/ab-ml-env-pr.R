############
# Plot partial resids for Ab and Ml (from phylo vint models), as a function of env
############


# load packages -----------------------------------------------------------

library(MCMCvis)
library(tidyverse)


# load model fits and data ------------------------------------------------

fit_ab <- readRDS(paste0(dir, '/Results/bird-ab-phylo-vint-', run_date,
                         '/bird-ab-phylo-vint-fit-', run_date, '.rds'))
data_ab <- readRDS(paste0(dir, '/Results/bird-ab-phylo-vint-', run_date,
                          '/bird-ab-phylo-vint-data-', run_date, '.rds'))
fit_ml <- readRDS(paste0(dir, '/Results/bird-ml-phylo-vint-', run_date,
                         '/bird-ml-phylo-vint-fit-', run_date, '.rds'))
data_ml <- readRDS(paste0(dir, '/Results/bird-ml-phylo-vint-', run_date,
                          '/bird-ml-phylo-vint-data-', run_date, '.rds'))


# extract posteriors, sim linear predictor, and get residuals -------------

alpha_mn_ab <- MCMCvis::MCMCpstr(fit_ab, params = 'alpha')[[1]]
kappa_mn_ab <- MCMCvis::MCMCpstr(fit_ab, params = 'kappa')[[1]]
gamma_mn_ab <- MCMCvis::MCMCpstr(fit_ab, params = 'gamma')[[1]]
beta_mn_ab <- MCMCvis::MCMCpstr(fit_ab, params = 'beta')[[1]]
mu_mn_ab <- kappa_mn_ab + gamma_mn_ab[data_ab$niche_idx] + alpha_mn_ab + (data_ab$X %*% beta_mn_ab)[,1]
resid_ab <- data_ab$Y - mu_mn_ab

alpha_mn_ml <- MCMCvis::MCMCpstr(fit_ml, params = 'alpha')[[1]]
kappa_mn_ml <- MCMCvis::MCMCpstr(fit_ml, params = 'kappa')[[1]]
gamma_mn_ml <- MCMCvis::MCMCpstr(fit_ml, params = 'gamma')[[1]]
beta_mn_ml <- MCMCvis::MCMCpstr(fit_ml, params = 'beta')[[1]]
mu_mn_ml <- kappa_mn_ml + gamma_mn_ml[data_ml$niche_idx] + alpha_mn_ml + (data_ml$X %*% beta_mn_ml)[,1]
resid_ml <- data_ml$Y - mu_mn_ml


# calculate partial resids for temp_sd_season -----------------------------

names <- c('Mass', 'Temp seasonality', 'Temp interannual',
           'Precip seasonality', 'Precip interannual')

NUM <- 3

pr_ab <- resid_ab + (beta_mn_ab[NUM] * data_ab$X[,NUM])
pr_ml <- resid_ml + (beta_mn_ml[NUM] * data_ml$X[,NUM])


# plot --------------------------------------------------------------------

tplt <- data.frame(var = data_ab$X[,NUM],
                   pr_ab = pr_ab, 
                   pr_ml = pr_ml)

#plot
ggplot(tplt, aes(pr_ml, pr_ab, color = var)) +
  geom_point(size = 2.5, alpha = 0.5) +
  theme_bw() +
  labs(color = names[NUM]) +
  xlab('Parital residuals Max Longevity') +
  ylab('Partial residuals Age First Breeding')

