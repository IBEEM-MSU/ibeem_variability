// gen lenth ~ env
// varying intercepts

data {
int<lower=0> N;                        // number of data points
int<lower=0> Nf;                      // number of groups
real<lower=0> y[N];
int<lower=0> f_id[N];
vector[N] lMass;
vector[N] temp_sd_season;
vector[N] temp_sd_year;
vector[N] temp_sp_color_month;
vector[N] precip_cv_season;
vector[N] precip_cv_year;
vector[N] precip_sp_color_month;
}

parameters {
real<lower=0> sigma;
real<lower=0> nu;
vector[Nf] alpha;
real mu_alpha;
real<lower=0> sigma_alpha;
real beta;
real gamma1;
real gamma2;
real gamma3;
real theta1;
real theta2;
real theta3;
}

transformed parameters {
vector[N] mu;

// intercept and slope for each species/site
mu = alpha[f_id] + beta * lMass + gamma1 * temp_sd_season + gamma2 * temp_sd_year + gamma3 * temp_sp_color_month + theta1 * precip_cv_season + theta2* precip_cv_year + theta3 * precip_sp_color_month;
}

model {
// priors
mu_alpha ~ std_normal();
sigma_alpha ~ std_normal();
beta ~ std_normal();
gamma1 ~ std_normal();
gamma2 ~ std_normal();
gamma3 ~ std_normal();
theta1 ~ std_normal();
theta2 ~ std_normal();
theta3 ~ std_normal();
sigma ~ std_normal();
nu ~ gamma(2, 0.1);             // degrees of freedom parameter

// centered
alpha ~ normal(mu_alpha, sigma_alpha);

y ~ student_t(nu, mu, sigma);
}

generated quantities {
vector[N] mu_nm;

mu_nm = gamma1 * temp_sd_season + gamma2 * temp_sd_year + gamma3 * temp_sp_color_month + theta1 * precip_cv_season + theta2* precip_cv_year + theta3 * precip_sp_color_month;
}
