// phylo intercepts
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> J;  // number of niche groups
  array[N] int<lower=1> niche_idx; // niche id of obs value
  matrix[N, N] LRho;  // cholesky factor of known correlation matrix
}
parameters {
  vector[K] beta;  // regression coefficients
  real kappa;  // global int
  vector[J] gamma; // intercept for niche
  real mu_gamma; // mean gammas
  real<lower=0> sigma_gamma; // mean gammas
  real<lower=0> sigma;  // process error
  real<lower=0> sigma_phylo;  // phylo sd
  vector[N] z;  // standardized group-level effects
}
transformed parameters {
  vector[N] alpha;  // phylo intercepts (per species)
  // implies alpha ~ MVN(0, Rho) * sigma_phylo
  alpha = (sigma_phylo * (LRho * z));
}
model {
  real lprior = 0;  // prior contributions to the log posterior
  // likelihood not including constants
  target += normal_lupdf(gamma | mu_gamma, sigma_gamma);
  
  // optimized call for glm
  target += normal_id_glm_lupdf(Y | X, kappa + gamma[niche_idx] + alpha, beta, sigma);

  // priors not including constants
  lprior += std_normal_lupdf(beta);
  lprior += normal_lupdf(kappa | 6, 3);
  lprior += std_normal_lupdf(sigma);
  lprior += std_normal_lupdf(sigma_phylo);
  lprior += std_normal_lupdf(mu_gamma);
  lprior += std_normal_lupdf(sigma_gamma);
  target += lprior;
  target += std_normal_lupdf(z);
}
