// gen lenth ~ env
// i = obs
// j = family
// y_ij ~ N(mu_ij, sigma)
// mu_ij = alpha_j + beta_j * x1_ij ...


data {
int<lower=0> N;                        // number of data points
int<lower=0> Nf;                      // number of families
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
real mu_alpha;
real mu_beta;
real mu_gamma1;
real mu_gamma2;
real mu_gamma3;
real mu_theta1;
real mu_theta2;
real mu_theta3;
vector<lower=0>[8] sigma_abgt;
cholesky_factor_corr[8] L_Rho;             // cholesky factor of corr matrix
matrix[8, Nf] z;                           // z-scores
}

transformed parameters {
vector[N] mu;
matrix[Nf, 8] abgt;
matrix[8, 8] Rho;               // corr matrix
vector[Nf] alpha;
vector[Nf] beta;
vector[Nf] gamma1;
vector[Nf] gamma2;
vector[Nf] gamma3;
vector[Nf] theta1;
vector[Nf] theta2;
vector[Nf] theta3;

// cholesky factor of covariance matrix multiplied by z score
// implies gt ~ MVN(0, sigma)
abgt = (diag_pre_multiply(sigma_abgt, L_Rho) * z)';
// implies Rho = L_Rho * L_Rho';
Rho = multiply_lower_tri_self_transpose(L_Rho);

alpha = mu_alpha + abgt[,1];
beta = mu_beta + abgt[,2];
gamma1 = mu_gamma1 + abgt[,3];
gamma2 = mu_gamma2 + abgt[,4];
gamma3 = mu_gamma3 + abgt[,5];
theta1 = mu_theta1 + abgt[,6];
theta2 = mu_theta2 + abgt[,7];
theta3 = mu_theta3 + abgt[,8];

// intercept and slope for each species/site
mu = alpha[f_id] + beta[f_id] .* lMass + gamma1[f_id] .* temp_sd_season + gamma2[f_id] .* temp_sd_year + gamma3[f_id] .* temp_sp_color_month + theta1[f_id] .* precip_cv_season + theta2[f_id] .* precip_cv_year + theta3[f_id] .* precip_sp_color_month;
}

model {
// priors
mu_alpha ~ std_normal();
mu_beta ~ std_normal();
mu_gamma1 ~ std_normal();
mu_gamma2 ~ std_normal();
mu_gamma3 ~ std_normal();
mu_theta1 ~ std_normal();
mu_theta2 ~ std_normal();
mu_theta3 ~ std_normal();
sigma_abgt ~ std_normal();
sigma ~ std_normal();
nu ~ gamma(2, 0.1);             // degrees of freedom parameter

to_vector(z) ~ std_normal();
L_Rho ~ lkj_corr_cholesky(2);

y ~ student_t(nu, mu, sigma);
}

generated quantities {
vector[N] mu_nm;

mu_nm = alpha + gamma1 * temp_sd_season + gamma2 * temp_sd_year + gamma3 * temp_sp_color_month + theta1 * precip_cv_season + theta2* precip_cv_year + theta3 * precip_sp_color_month;
}
