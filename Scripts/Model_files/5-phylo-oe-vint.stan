// phylo intercepts, observation model, varying int by niche
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> No;
  int<lower=1> Ni;
  vector[No] y_obs;  // response variable
  vector[Ni] y_imp;  // response variable
  vector[Ni] sd_y;
  int<lower=1> K;  // number of population-level effects
  int<lower=1> J;  // number of niche groups
  matrix[No, K] X_obs;  // population-level design matrix
  matrix[Ni, K] X_imp;  // population-level design matrix
  array[No] int<lower=1> obs_idx; // idx of alpha that corresponds to obs value
  array[Ni] int<lower=1> imp_idx; // idx of alpha that corresponds to imp value
  array[No] int<lower=1> niche_obs_idx; // niche id of obs value
  array[Ni] int<lower=1> niche_imp_idx; // niche id of imp value
  matrix[N, N] LRho;  // cholesky factor of known correlation matrix
}
parameters {
  vector[K] beta;  // regression coefficients
  real kappa;  // global int
  vector[J] gamma; // intercept for niche
  real mu_gamma; // mean gammas
  real<lower=0> sigma_gamma; // mean gammas
  real<lower=0> sigma;  // process error
  real<lower=0> sigma_phylo;  // phylo sd
  vector[N] z;  // standardized group-level effects
  vector[Ni] y_iv_raw;
}
transformed parameters {
  vector[N] alpha;  // phylo intercepts (per species)
  vector[Ni] y_iv;
  
  // implies alpha ~ MVN(0, Rho) * sigma_phylo
  alpha = (sigma_phylo * (LRho * z));
  
  // non-centered - implies y_iv ~ normal(X_imp * beta, sigma);
  y_iv = y_iv_raw * sigma + kappa + gamma[niche_imp_idx] + alpha[imp_idx] + X_imp * beta;
}
model {
  real lprior = 0;  // prior contributions to the log posterior
  // likelihood not including constants
  target += normal_lupdf(gamma | mu_gamma, sigma_gamma);
  
  // imputed data - observation model (y_iv equiv to y_obs)
  target += normal_lupdf(y_imp | y_iv, sd_y);

  // optimized call for glm
  // observed data - process model - could be t-dist
  target += normal_id_glm_lupdf(y_obs | X_obs, kappa + gamma[niche_obs_idx] + alpha[obs_idx], beta, sigma);

  // priors not including constants
  lprior += std_normal_lupdf(beta);
  lprior += normal_lupdf(kappa | 5, 2);
  lprior += std_normal_lupdf(sigma);
  lprior += std_normal_lupdf(sigma_phylo);
  // lprior += normal_lupdf(sigma_phylo | 0, 2);
  lprior += std_normal_lupdf(mu_gamma);
  lprior += std_normal_lupdf(sigma_gamma);
  target += lprior;
  target += std_normal_lupdf(z);
  target += std_normal_lupdf(y_iv_raw);
}
