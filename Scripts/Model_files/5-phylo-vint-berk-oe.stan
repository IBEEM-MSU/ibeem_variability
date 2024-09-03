
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y_obs;  // response variable
  real<lower=0> sd_Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  vector[N] lMass;
  vector[N] temp_sd_season;
  vector[N] temp_sd_year;
  vector[N] precip_cv_season;
  vector[N] precip_cv_year;
  vector<lower=0>[N] temp_sd_season_sd_space;
  vector<lower=0>[N] temp_sd_year_sd_space;
  vector<lower=0>[N] precip_cv_season_sd_space;
  vector<lower=0>[N] precip_cv_year_sd_space;
  int<lower=1> J;  // number of niche groups
  array[N] int<lower=1> niche_idx; // niche id of obs value
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
  vector<lower=0>[N] tr_temp_sd_season; // true covariate (Berkson)
  vector<lower=0>[N] tr_temp_sd_year;
  vector<lower=0>[N] tr_precip_cv_season;
  vector<lower=0>[N] tr_precip_cv_year;
  vector[N] Y_t_raw;
}
transformed parameters {
  vector[N] alpha;  // phylo intercepts (per species)
  vector[N] Y_t;
  
  // implies alpha ~ MVN(0, Rho) * sigma_phylo^2
  alpha = (sigma_phylo^2 * (LRho * z));
  
  // implies y_t ~ normal(kappa + gamma + ..., sigma)
  Y_t = Y_t_raw * sigma + 
  kappa + 
  gamma[niche_idx] + 
  alpha + 
  beta[1] * lMass +
  beta[2] * tr_temp_sd_season +
  beta[3] * tr_temp_sd_year +
  beta[4] * tr_precip_cv_season +
  beta[5] * tr_precip_cv_year;
}
model {
  // priors
  beta ~ std_normal();
  kappa ~ normal(mu_kappa, sigma_kappa); // hard coded mu_kappa and sigma_kappa
  sigma ~ std_normal();
  sigma_phylo ~ std_normal();
  sigma_gamma ~ std_normal();
  z ~ std_normal();
  Y_t_raw ~ std_normal();
  
  // niche intercepts
  gamma ~ normal(0, sigma_gamma);
  
  // true ~ N(observed, known_spatial_var)
  // to account for aggregation over species range
  tr_temp_sd_season ~ normal(temp_sd_season, temp_sd_season_sd_space);
  tr_temp_sd_year ~ normal(temp_sd_year, temp_sd_year_sd_space);
  tr_precip_cv_season ~ normal(precip_cv_season, precip_cv_season_sd_space);
  tr_precip_cv_year ~ normal(precip_cv_year, precip_cv_year_sd_space);
  
  Y_obs ~ normal(Y_t, sd_Y);
}
