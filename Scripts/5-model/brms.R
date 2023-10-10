library(ape)
library(brms)
library(cmdstanr)

# # https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html
# phylo <- ape::read.nexus("https://paul-buerkner.github.io/data/phylo.nex")
# data_simple <- read.table(
#   "https://paul-buerkner.github.io/data/data_simple.txt",
#   header = TRUE)
# 
# A <- ape::vcv.phylo(phylo)
# 
# model_simple <- brms::brm(
#   phen ~ cofactor + (1|gr(phylo, cov = A)),
#   data = data_simple,
#   family = gaussian(),
#   data2 = list(A = A),
#   prior = c(
#     prior(normal(0, 10), "b"),
#     prior(normal(0, 50), "Intercept"),
#     prior(student_t(3, 0, 20), "sd"),
#     prior(student_t(3, 0, 20), "sigma")))
# 
# summary(model_simple)
# plot(model_simple, N = 2, ask = FALSE)
# 
# stancode(model_simple)

DATA <- list(N = NROW(bird_df5),
             R = R,
             I = diag(1, nrow = NROW(bird_df5), ncol = NROW(bird_df5)),
             y_obs = bird_df5$Phylo_survival,
             lMass_obs = bird_df5$lMass,
             temp_sd_season_obs = bird_df5$temp_sd_season,
             temp_sd_year_obs = bird_df5$temp_sd_year,
             precip_cv_season_obs = bird_df5$precip_cv_season,
             precip_cv_year_obs = bird_df5$precip_cv_year,
             pro_data = bird_df5)

#use corr matrix: https://discourse.mc-stan.org/t/covariance-matrix-phylogenetic-models/20477/4
V <- ape::vcv.phylo(tree_n2, corr = TRUE)

bird_df5$sc_PS <- bird_df5$Phylo_survival * 100

tt <- data.frame(lMass = bird_df5$lMass,
                 temp_sd_season = bird_df5$temp_sd_season,
                 temp_sd_year = bird_df5$temp_sd_year,
                 precip_cv_season = bird_df5$precip_cv_season,
                 precip_cv_year = bird_df5$precip_cv_year) %>%
  apply(2, function(x) scale(x, scale = TRUE))

dd <- data.frame(species = bird_df5$species,
  sc_PS = bird_df5$sc_PS,
  tt)

#https://discourse.mc-stan.org/t/running-speed-phylogenetic-models/30898
model_formula <- "sc_PS ~ lMass + 
    temp_sd_season + 
    temp_sd_year +
    precip_cv_season +
    precip_cv_year + 
    (1|gr(species, cov = V))"

# brms formula.
brms_formula <- brmsformula(model_formula)#, 
                            # family = cumulative(threshold = "equidistant"), 
                            # decomp = "QR")

brms_model <- brm(
  brms_formula,
  data = dd,
  data2 = list(V = V),
  prior = c(
    prior(normal(0, 1), "b"),
    prior(normal(0, 1), "Intercept"),
    prior(normal(0, 1), "sd"),
    prior(normal(0, 1), "sigma")),
  iter = 4000,
  chains = 4,
  cores = 4,
  normalize = FALSE,
  control = list(adapt_delta = 0.999,
                 stepsize = 0.001))#,
  # backend = "cmdstanr")#,
  # control = list(max_treedepth = 5),
  # threads = threading(25))

summary(brms_model)
# library(shinystan)
shinystan::launch_shinystan(brms_model$fit)

# 
# m1 <- brms::brm(
#   sc_PS ~ lMass + 
#     temp_sd_season + 
#     temp_sd_year +
#     precip_cv_season +
#     precip_cv_year + 
#     (1|gr(species, cov = V)),
#   data = bird_df5,
#   family = gaussian(),
#   decomp = 'QR',
#   data2 = list(V = V),
#   prior = c(
#     prior(normal(0, 1), "b"),
#     prior(normal(0, 1), "Intercept"),
#     prior(normal(0, 1), "sd"),
#     prior(normal(0, 1), "sigma")),
#     # prior(student_t(3, 0, 1), "sd"),
#     # prior(student_t(3, 0, 1), "sigma")),
#   cores = 4,
#   chains = 4,
#   iter = 5000)

summary(brms_model)
NROW(bird_df5)
stancode(brms_model)
brms::standata(brms_model)

str(brms_model)


# matrix with ones in first col
# center predictors
# center response


tt <- data.frame(bird_df5$lMass,
           bird_df5$temp_sd_season,
           bird_df5$temp_sd_year,
           bird_df5$precip_cv_season,
           bird_df5$precip_cv_year) %>%
  apply(2, function(x) scale(x, scale = TRUE))

tt2 <- cbind(1, tt)

DATA <- list(N = NROW(bird_df5),
             K = NCOL(tt2),
             y_obs = scale(bird_df5$Phylo_survival * 10, scale = FALSE)[,1],
             x = tt2,
             L_Rho = chol(R), #cholesky factor of corr (distance matrix)
             zeros = rep(0, NROW(bird_df5)),
             pro_data = bird_df5)

# DELTA <- 0.95
# TREE_DEPTH <- 12
# STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 2000

# XXXX sec N = 100 - slow, also doesn't fit well
fit <- rstan::stan(paste0(dir, 'Scripts/Model_files/5-novar-phylo.stan'),
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('beta',
                            'sigma',
                            'alpha',
                            'sigma_phy'))#,
                   # control = list(adapt_delta = DELTA,
                   #                max_treedepth = TREE_DEPTH,
                   #                stepsize = STEP_SIZE))

MCMCsummary(fit)
shinystan::launch_shinystan(fit)




#brms modified model

# data {
#   int<lower=1> N;  // total number of observations
#   vector[N] Y;  // response variable
#   int<lower=1> K;  // number of population-level effects
#   matrix[N, K] X;  // population-level design matrix
#   matrix[N, N] Lcov;  // cholesky factor of known covariance matrix
# }

V <- ape::vcv.phylo(tree_n2, corr = TRUE)

tt <- data.frame(bird_df5$lMass / 2,
                 bird_df5$temp_sd_season / 20,
                 bird_df5$temp_sd_year / 10,
                 bird_df5$precip_cv_season / 10,
                 bird_df5$precip_cv_year / 10) %>%
  apply(2, function(x) scale(x, scale = FALSE))

DATA <- list(N = NROW(bird_df5),
             # Y = scale(bird_df5$Phylo_survival * 10, scale = FALSE)[,1],
             Y = bird_df5$Phylo_survival * 10,
             K = NCOL(tt),
             X = tt,
             LRho = chol(V)) #cholesky factor of corr matrix

DELTA <- 0.92
TREE_DEPTH <- 12
STEP_SIZE <- 0.03
CHAINS <- 4
ITER <- 3000

# N = 500 - 550.415 sec (9 min) - no fit problems
# N = 2000 - 11971.2 sec (200 min) - no fit problems

#rstan with cor mat fastest
# N = 250; scale y, corr mat - 32 sec, 29 sec, 31 sec
# N = 250; no center y, corr mat - 32 sec, 28 sec, 30 sec
# N = 250; scale y, cov mat - 36 sec, 32 sec, 28 sec
# N = 250; no center y, cov mat - 31 sec, 37 sec, 50 sec
# N = 250; no optimized lm fun - 48, 32, 35 
# N = 250; no mu - 28, 30
# N = 250; no scale or center y - 44
# N = 250; scale 100 y - 23, 32, 52
# N = 250; no y scale w/ app priors - 58, 59, 62
# N = 250; scale 10 w/ app priors - 46 sec, 50
# N = 500; scale 10 w/ app priors - 400 sec (6-7 min)

#brms
# N = 250; scale y, corr mat - 342 sec, 243 sec
# N = 250; no center y, corr mat - 265 sec
# N = 250; scale y, cov mat - 
# N = 250; no center y, cov mat - 193 sec
fit <- rstan::stan(paste0(dir, 'Scripts/Model_files/brms_mod.stan'),
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('beta',
                            'sigma',
                            'sigma_phylo',
                            'kappa'),
                   control = list(adapt_delta = DELTA,
                   max_treedepth = TREE_DEPTH,
                   stepsize = STEP_SIZE))

shinystan::launch_shinystan(fit)
saveRDS(fit, '~/Desktop/brms_mod_stan.rds')

tt <- MCMCvis::MCMCpstr(fit, c('sigma', 'sigma_phylo'))

#lambda ~ 0.7
tt$sigma_phylo / (tt$sigma + tt$sigma_phylo)

# dd <- data.frame(species = bird_df5$species,
#                  # sc_PS = scale(bird_df5$Phylo_survival * 10, scale = FALSE)[,1],
#                  sc_PS = bird_df5$Phylo_survival * 10,
#                  lMass = tt[,1],
#                  temp_sd_season = tt[,2],
#                  temp_sd_year = tt[,3],
#                  precip_cv_season = tt[,4],
#                  precip_cv_year = tt[,5])
#            
# model_formula <- "sc_PS ~ lMass + 
#     temp_sd_season + 
#     temp_sd_year +
#     precip_cv_season +
#     precip_cv_year + 
#     (1|gr(species, cov = V))"
# 
# # brms formula.
# brms_formula <- brmsformula(model_formula)#, 
# # family = cumulative(threshold = "equidistant"), 
# # decomp = "QR")
# 
# brms_model <- brm(
#   brms_formula,
#   data = dd,
#   data2 = list(V = V),
#   prior = c(
#     prior(normal(0, 1), "b"),
#     prior(normal(0, 1), "Intercept"),
#     prior(normal(0, 1), "sd"),
#     prior(normal(0, 1), "sigma")),
#   iter = ITER,
#   chains = CHAINS,
#   cores = CHAINS,
#   normalize = FALSE,
#   control = list(adapt_delta = DELTA,
#                  max_treedepth = TREE_DEPTH,
#                  stepsize = STEP_SIZE))
# 
# str(brms_model)


bird_df5

obs_idx <- which(bird_df5$SD_survival == 0)
imp_idx <- which(bird_df5$SD_survival != 0)

V <- ape::vcv.phylo(tree_n2, corr = TRUE)

tt_obs <- data.frame(bird_df5$lMass[obs_idx] / 2,
                 bird_df5$temp_sd_season[obs_idx] / 15,
                 bird_df5$temp_sd_year[obs_idx] / 20,
                 bird_df5$precip_cv_season[obs_idx] / 20,
                 bird_df5$precip_cv_year[obs_idx] / 20) %>%
  apply(2, function(x) scale(x, scale = FALSE))

tt_imp <- data.frame(bird_df5$lMass[imp_idx] / 2,
                     bird_df5$temp_sd_season[imp_idx] / 15,
                     bird_df5$temp_sd_year[imp_idx] / 20,
                     bird_df5$precip_cv_season[imp_idx] / 20,
                     bird_df5$precip_cv_year[imp_idx] / 20) %>%
  apply(2, function(x) scale(x, scale = FALSE))

DATA <- list(N = NROW(bird_df5),
             No = length(obs_idx),
             Ni = length(imp_idx),
             y_obs = bird_df5$Phylo_survival[obs_idx] * 10,
             y_imp = bird_df5$Phylo_survival[imp_idx] * 10,
             sd_y = bird_df5$SD_survival[imp_idx] * 10,
             K = NCOL(tt_obs),
             X_obs = tt_obs,
             X_imp = tt_imp,
             imp_idx = imp_idx,
             obs_idx = obs_idx,
             LRho = chol(V)) #cholesky factor of cov matrix

#scales as cube?
#N = 500 - 148 sec (2.5 min)
#N = 1000 - 1141 sec (19 min)
fit <- rstan::stan(paste0(dir, 'Scripts/Model_files/brms_mod_oe.stan'),
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('beta',
                            'sigma',
                            'sigma_phylo',
                            'kappa'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))




library(brms)

idx <- sample(1:NROW(bird_df5), 100)
bird_df6 <- bird_df5[idx,]

#load consensus tree - bird.phylo
load(paste0(dir, 'data/L3/bird-consensus-tree.rda'))

#df with names and idx
idx_df <- data.frame(idx = 1:NROW(bird_df6),
                     name = bird_df6$species)

#species not found in both datasets (species to drop from tree)
nm <- setdiff(bird.phylo$tip.label, bird_df6$species)

#prune specified tips from all trees
pr_tree2 <- ape::drop.tip(bird.phylo, nm)

V <- ape::vcv.phylo(pr_tree2, corr = TRUE)

model_formula <- "lGL ~ lMass + 
    temp_sd_season + 
    temp_sd_year +
    precip_cv_season +
    precip_cv_year + 
    (1|gr(species, cov = V))"

# brms formula.
brms_formula <- brmsformula(model_formula)#, 
# family = cumulative(threshold = "equidistant"), 
# decomp = "QR")

brms_model <- brm(
  brms_formula,
  data = bird_df6,
  data2 = list(V = V),
  iter = 2000,
  chains = 4,
  cores = 4,
  normalize = FALSE,
  backend = 'cmdstanr')


fit_parallel <- update(
  brms_model, 
  chains = 2, 
  cores = 2,
  backend = "cmdstanr", 
  threads = threading(2))


stancode(brms_model)
stancode(fit_parallel)

brms::standata(brms_model)
brms::standata(fit_parallel)
