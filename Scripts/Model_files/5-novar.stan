// LH trait ~ env

data {
int<lower=0> N;                        // number of observed vals
real y[N];                          // observed vals
int<lower=1> K;  // number of population-level effects
matrix[N, K] X;  // population-level design matrix
}

parameters {
real<lower=0> sigma;
vector[K] beta;  // regression coefficients
}

transformed parameters {
vector[N] mu;

// predictors
mu = X * beta;
}

model {
// priors
beta ~ std_normal();
sigma ~ std_normal();

// observed data - process model - could be t-dist
y ~ normal(mu, sigma);
}

// generated quantities {
// vector[N] mu_nm;
// 
// mu_nm = alpha + gamma1 * temp_sd_season + gamma2 * temp_sd_year + gamma3 * temp_sp_color_month + theta1 * precip_cv_season + theta2* precip_cv_year + theta3 * precip_sp_color_month;
// }
