// gen lenth ~ env
// i = obs
// j = family
// y_ij ~ N(mu_ij, sigma)
// mu_ij = alpha_j + beta_j * x1_ij ...

data {
int<lower=0> N;                        // number of data points
int<lower=0> J;                      // number of families
vector[N] y;
array[N] int<lower=1, upper=J> niche_idx; // niche id of obs value
vector[N] lMass;
vector[N] temp_sd_season;
vector[N] temp_sd_year;
vector[N] precip_cv_season;
vector[N] precip_cv_year;
}

parameters {
real<lower=0> sigma;
real<lower=0> nu;
real mu_alpha;
real mu_beta;
real mu_gamma1;
real mu_gamma2;
real mu_theta1;
real mu_theta2;
vector<lower=0>[6] sigma_abgt;
cholesky_factor_corr[6] L_Rho;             // cholesky factor of corr matrix
matrix[6, J] z;                           // z-scores
}

transformed parameters {
vector[N] mu;
matrix[J, 6] abgt;
matrix[6, 6] Rho;               // corr matrix
vector[J] alpha;
vector[J] beta;
vector[J] gamma1;
vector[J] gamma2;
vector[J] theta1;
vector[J] theta2;

// cholesky factor of covariance matrix multiplied by z score
// implies gt ~ MVN(0, sigma)
abgt = (diag_pre_multiply(sigma_abgt, L_Rho) * z)';
// implies Rho = L_Rho * L_Rho';
Rho = multiply_lower_tri_self_transpose(L_Rho);

alpha = mu_alpha + abgt[,1];
beta = mu_beta + abgt[,2];
gamma1 = mu_gamma1 + abgt[,3];
gamma2 = mu_gamma2 + abgt[,4];
theta1 = mu_theta1 + abgt[,5];
theta2 = mu_theta2 + abgt[,6];

// intercept and slope for each species/site
mu = alpha[niche_idx] + beta[niche_idx] .* lMass + gamma1[niche_idx] .* temp_sd_season + gamma2[niche_idx] .* temp_sd_year + theta1[niche_idx] .* precip_cv_season + theta2[niche_idx] .* precip_cv_year;
}

model {
// priors
mu_alpha ~ std_normal();
mu_beta ~ std_normal();
mu_gamma1 ~ std_normal();
mu_gamma2 ~ std_normal();
mu_theta1 ~ std_normal();
mu_theta2 ~ std_normal();
sigma_abgt ~ std_normal();
sigma ~ std_normal();
nu ~ gamma(2, 0.1);             // degrees of freedom parameter

to_vector(z) ~ std_normal();
L_Rho ~ lkj_corr_cholesky(2);

y ~ student_t(nu, mu, sigma);
}

