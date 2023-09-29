// gen lenth ~ env
// no varying intercepts or slopes - incorporate obs error
// separate obs and imputed vales

data {
int<lower=0> No;                        // number of observed vals
int<lower=0> Ni;                        // number of imputed vals
real y_obs[No];                          // observed vals
real y_imp[Ni];                          // imputed vals
real<lower=0> sd_y[Ni];                  // uncertaity of imputed vals
vector[No] lMass_obs;                   // predictors associated with observed vals
vector[No] temp_sd_season_obs;
vector[No] temp_sd_year_obs;
// vector[N] temp_sp_color_month;
vector[No] precip_cv_season_obs;
vector[No] precip_cv_year_obs;
// vector[N] precip_sp_color_month;
vector[Ni] lMass_imp;
vector[Ni] temp_sd_season_imp;
vector[Ni] temp_sd_year_imp;
// vector[N] temp_sp_color_month;
vector[Ni] precip_cv_season_imp;
vector[Ni] precip_cv_year_imp;
}

parameters {
real<lower=0> sigma;
// real<lower=0> nu;
real alpha;
real beta;
real gamma1;
real gamma2;
// real gamma3;
real theta1;
real theta2;
// real theta3;
// vector[N] y;
vector[Ni] y_iv_raw;
}

transformed parameters {
vector[No] mu_obs;
vector[Ni] mu_imp;
vector[Ni] y_iv;

// predictors obs data
mu_obs = alpha + 
beta * lMass_obs + 
gamma1 * temp_sd_season_obs + 
gamma2 * temp_sd_year_obs + 
theta1 * precip_cv_season_obs + 
theta2 * precip_cv_year_obs;

// predictors imp data - same params as for obs data
mu_imp = alpha + 
beta * lMass_imp + 
gamma1 * temp_sd_season_imp + 
gamma2 * temp_sd_year_imp + 
theta1 * precip_cv_season_imp + 
theta2 * precip_cv_year_imp;

// non-centered for imputed data - same process error as obs
// equiv to y_iv ~ normal(mu_imp, sigma)
y_iv = y_iv_raw * sigma + mu_imp;
}

model {
// priors
alpha ~ std_normal();
beta ~ std_normal();
gamma1 ~ std_normal();
gamma2 ~ std_normal();
// gamma3 ~ std_normal();
theta1 ~ std_normal();
theta2 ~ std_normal();
// theta3 ~ std_normal();
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
