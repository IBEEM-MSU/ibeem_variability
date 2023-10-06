// LH trait ~ env
// no varying intercepts or slopes - incorporate phylo
// no obs err

data {
int<lower=0> N;                        // number of observed vals
real y_obs[N];                          // observed vals
vector[No] lMass_obs;                   // predictors associated with observed vals
vector[No] temp_sd_season_obs;
vector[No] temp_sd_year_obs;
// vector[N] temp_sp_color_month;
vector[No] precip_cv_season_obs;
vector[No] precip_cv_year_obs;
// vector[N] precip_sp_color_month;
matrix[N, N] R;                   // scaled phylo distance matrix (from phylogeny)
matrix[N, N] I;                   // identity matrix
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
real<lower=0, upper=1> lambda;
}

transformed parameters {
vector[N] mu_obs;
matrix[N, N] R_lambda;
matrix[N, N] S;

// predictors obs data
mu_obs = alpha + 
beta * lMass_obs + 
gamma1 * temp_sd_season_obs + 
gamma2 * temp_sd_year_obs + 
theta1 * precip_cv_season_obs + 
theta2 * precip_cv_year_obs;

// Pagel's lambda
R_lambda = lambda * R + (1 - lambda) * I;
S = R_lambda * sigma;

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
lambda ~ uniform(0, 1);
// nu ~ gamma(2, 0.1);             // degrees of freedom parameter

y_obs ~ multi_normal(mu_obs, S); 
}

