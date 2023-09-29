// gen lenth ~ env
// no varying intercepts or slopes - incorporate obs error

data {
int<lower=0> N;                        // number of data points
real y_obs[N];
real<lower=0> sd_y[N];
vector[N] lMass;
vector[N] temp_sd_season;
vector[N] temp_sd_year;
// vector[N] temp_sp_color_month;
vector[N] precip_cv_season;
vector[N] precip_cv_year;
// vector[N] precip_sp_color_month;
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
vector[N] y_raw;
}

transformed parameters {
vector[N] mu;
vector[N] y;

// intercept and slope for each species/site
mu = alpha + beta * lMass + gamma1 * temp_sd_season + gamma2 * temp_sd_year + theta1 * precip_cv_season + theta2 * precip_cv_year;

// non-centered
y = y_raw * sigma + mu;
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
y_raw ~ std_normal();
// nu ~ gamma(2, 0.1);             // degrees of freedom parameter

// process model
// y ~ student_t(nu, mu, sigma);

// observation model
y_obs ~ normal(y, sd_y);
}

// generated quantities {
// vector[N] mu_nm;
// 
// mu_nm = alpha + gamma1 * temp_sd_season + gamma2 * temp_sd_year + gamma3 * temp_sp_color_month + theta1 * precip_cv_season + theta2* precip_cv_year + theta3 * precip_sp_color_month;
// }
