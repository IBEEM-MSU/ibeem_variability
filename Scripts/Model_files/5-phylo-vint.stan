// phylo intercepts, intercepts by niche group
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
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
  vector[N] theta;  // unscaled phylo intecepts
  // vector[N] z;  // standardized group-level effects
}
transformed parameters {
  vector[N] alpha;  // phylo intercepts (per species)
  // implies alpha ~ MVN(0, Rho) * sigma_phylo^2
  // alpha = (sigma_phylo^2 * (LRho * z));
  
  // scale phylo intercepts
  alpha = theta * sigma_phylo;
}
model {
  // priors
  beta ~ std_normal();
  kappa ~ normal(mu_kappa, sigma_kappa); // hard coded mu_kappa and sigma_kappa
  sigma ~ std_normal();
  sigma_phylo ~ std_normal();
  sigma_gamma ~ std_normal();
  // z ~ std_normal();
  
  // niche intercepts
  gamma ~ normal(0, sigma_gamma);
  // phylo intercepts (unscaled)
  theta ~ multi_normal_cholesky(0, LRho)
  
  // optimized call for glm
  Y ~ normal_id_glm(X, kappa + gamma[niche_idx] + alpha, beta, sigma);
  // Y ~ normal(kappa + gamma[niche_idx] + alpha + X * beta, sigma);
}
