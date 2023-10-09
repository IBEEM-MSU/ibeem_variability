// LH trait ~ env - TOO SLOW SEE 5-phylo.stan or 5-phylo-oe.stan for optimized scripts
// no varying intercepts or slopes - incorporate phylo
// no obs err -

data {
int<lower=0> N;                        // number of observed vals
int<lower=0> K;                        // number of predictors(5) + 1 
vector[N] y_obs;                          // observed vals
matrix[N, K] x;                          // predictor matrix
vector[N] zeros;
// vector[N] lMass_obs;                   // predictors associated with observed vals
// vector[N] temp_sd_season_obs;
// vector[N] temp_sd_year_obs;
// vector[N] temp_sp_color_month;
// vector[N] precip_cv_season_obs;
// vector[N] precip_cv_year_obs;
// vector[N] precip_sp_color_month;
// matrix[N, N] R;                   // scaled phylo distance matrix (from phylogeny)
// matrix[N, N] I;                   // identity matrix
matrix[N, N] L_Rho;                 // cholesky factor of known corr (scaled distance) matrix
// matrix[N, N] L_cov;    // cholesky factor of known cov
}

parameters {
real<lower=0> sigma;
// real<lower=0> nu;
vector[K] beta;
real<lower=0> sigma_phy;
vector[N] alpha;
// real alpha;
// real beta;
// real gamma1;
// real gamma2;
// // real gamma3;
// real theta1;
// real theta2;
// real theta3;
// vector[N] y;
// real<lower=0, upper=1> lambda;
// vector[N] eps;
// vector[N] z;
}

transformed parameters {
// vector[N] mu_obs;
// matrix[N, N] R_lambda;
matrix[N, N] S;
// cholesky_factor_corr[N] L_lambda;
// cholesky_factor_cov[N] L;

// Pagel's lambda - mixture of R (scaled distance [corr] mat) and I
// R_lambda = lambda * R + (1 - lambda) * I;

// https://archive.ph/ZeSim - Cholesky factors of covariance and correlation matrices in Stan
// approach one: get Cholesky factor corr (L_lambda) and compute Cholesky cov (L)
// L_lambda = cholesky_decompose(R_lambda);
// L = sigma * L_lambda;
// vector[N] alpha;

// species-specific intercepts
// alpha = ((sigma^2) * L_Rho) * z; // diag... gives cholesky of cov_mat whre L_Rho is cholesky of corr mat
// y = gamma + beta * X + (lambda * alpha + (1 - lambda) * eps)

// approach two: get cov and compute Cholesky cov
// multiply correlation matrix by variance to get covariance matrix
// S = R_lambda * sigma2;
// ^^ equiv to S = quad_form_diag(R_lambda, sqrt(sigma2));
// get Cholesky cov
// L = cholesky_decompose(S);

// predictors obs data
// mu_obs = 
// beta * lMass_obs + 
// gamma1 * temp_sd_season_obs + 
// gamma2 * temp_sd_year_obs + 
// theta1 * precip_cv_season_obs + 
// theta2 * precip_cv_year_obs;

S = L_Rho * sigma_phy;
}

model {
// priors
// alpha ~ std_normal();
beta ~ std_normal();
// gamma1 ~ std_normal();
// gamma2 ~ std_normal();
// gamma3 ~ std_normal();
// theta1 ~ std_normal();
// theta2 ~ std_normal();
// theta3 ~ std_normal();
sigma ~ std_normal();
sigma_phy ~ std_normal();
// lambda ~ uniform(0, 1);
// kappa ~ std_normal();
// z ~ std_normal();
// nu ~ gamma(2, 0.1);             // degrees of freedom parameter

// y_obs ~ multi_normal_cholesky(mu_obs, L);
// y_obs ~ multi_normal(x * beta, S);

// eps ~ normal(0, sigma);

// mixture of covariance with phylo and just sigma
// y_obs = x * beta + (lambda * alpha + (1 - lambda) * eps);
// y_obs ~ normal(x * beta + (lambda * alpha), (1 - lambda) * sigma);

alpha ~ multi_normal_cholesky(zeros, S);
y_obs ~ normal_id_glm(x, alpha, beta, sigma);
}

