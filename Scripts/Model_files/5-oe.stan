// LH trait ~ env
// no varying intercepts or slopes - incorporate obs error
// separate obs and imputed vales

data {
int<lower=0> No;                        // number of observed vals
int<lower=0> Ni;                        // number of imputed vals
real y_obs[No];                          // observed vals
real y_imp[Ni];                          // imputed vals
real<lower=0> sd_y[Ni];                  // uncertaity of imputed vals
int<lower=1> K;  // number of population-level effects
matrix[No, K] X_obs;  // population-level design matrix
matrix[Ni, K] X_imp;  // population-level design matrix
int<lower=1> obs_idx[No]; // idx of alpha that corresponds to obs value
int<lower=1> imp_idx[Ni]; // idx of alpha that corresponds to imp value
}

parameters {
real<lower=0> sigma;
vector[K] beta;  // regression coefficients
vector[Ni] y_iv_raw;
}

transformed parameters {
vector[No] mu_obs;
vector[Ni] mu_imp;
vector[Ni] y_iv;

// predictors obs data
mu_imp = X_imp * beta;

// predictors imp data - same params as for obs data
mu_obs = X_obs * beta;

// non-centered for imputed data - same process error as obs
// equiv to y_iv ~ normal(mu_imp, sigma)
y_iv = y_iv_raw * sigma + mu_imp;
}

model {
// priors
beta ~ std_normal();
sigma ~ std_normal();
y_iv_raw ~ std_normal();
// nu ~ gamma(2, 0.1);             // degrees of freedom parameter

// centered imputed data - same process error as obs
// y_iv ~ normal(mu_imp, sigma)

// imputed data - observation model (y_iv equiv to y_obs)
y_imp ~ normal(y_iv, sd_y);

// observed data - process model - could be t-dist
y_obs ~ normal(mu_obs, sigma);
}

// generated quantities {
// vector[N] mu_nm;
// 
// mu_nm = alpha + gamma1 * temp_sd_season + gamma2 * temp_sd_year + gamma3 * temp_sp_color_month + theta1 * precip_cv_season + theta2* precip_cv_year + theta3 * precip_sp_color_month;
// }
