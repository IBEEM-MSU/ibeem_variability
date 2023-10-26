// phylo intercepts, intercepts by niche group, oe
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> No;
  int<lower=1> Nm;
  vector[No] Y_obs;  // response variable - observed
  vector[Nm] Y_mod;  // response variable - modeled
  real sd_Y; // uncertainty
  int<lower=1> K;  // number of population-level effects
  int<lower=1> J;  // number of niche groups
  matrix[No, K] X_obs;  // design matrix - for observed data
  matrix[Nm, K] X_mod;  // design matrix - for modeled data
  array[No] int<lower=1> obs_idx; // idx of alpha that corresponds to obs value
  array[Nm] int<lower=1> mod_idx; // idx of alpha that corresponds to mod value
  array[No] int<lower=1> niche_obs_idx; // niche id of obs value
  array[Nm] int<lower=1> niche_mod_idx; // niche id of imp value
  matrix[N, N] Rho;  // known correlation matrix
  real mu_kappa;  // to set prior mu_kappa
  real sigma_kappa;  // to set prior sigma_kappa
}
transformed data {
  matrix[N, N] LRho = cholesky_decompose(Rho); // get cholesky factor of known corr matrix
}
parameters {
  vector[K] beta;  // regression coefficients
  real kappa;  // global int
  vector[J] gamma; // intercept for niche
  real<lower=0> sigma_gamma; // mean gammas
  real<lower=0> sigma;  // process error
  real<lower=0> sigma_phylo;  // phylo sd
  vector[N] z;  // standardized group-level effects
  vector[Nm] Y_mv_raw;
}
transformed parameters {
  vector[N] alpha;  // phylo intercepts (per species)
  vector[Nm] Y_mv;
  
  // implies alpha ~ MVN(0, Rho) * sigma_phylo^2
  alpha = (sigma_phylo^2 * (LRho * z));
  
  // implies y_iv ~ normal(kappa + gamma + ..., sigma)
  Y_mv = Y_mv_raw * sigma + kappa + gamma[niche_mod_idx] + alpha[mod_idx] + X_mod * beta;

}
model {
  // priors
  beta ~ std_normal();
  kappa ~ normal(mu_kappa, sigma_kappa); // hard coded mu_kappa and sigma_kappa
  sigma ~ std_normal();
  sigma_phylo ~ std_normal();
  sigma_gamma ~ std_normal();
  Y_mv_raw ~ std_normal();
  z ~ std_normal();
  
  // niche intercepts
  gamma ~ normal(0, sigma_gamma);
  
  // modeled data - observation model (Y_mv equiv to Y_obs)
  Y_mod ~ normal(Y_mv, sd_Y);
  
  // optimized call for glm
  // observed data - process model
  Y_obs ~ normal_id_glm(X_obs, kappa + gamma[niche_obs_idx] + alpha[obs_idx], beta, sigma);
}
