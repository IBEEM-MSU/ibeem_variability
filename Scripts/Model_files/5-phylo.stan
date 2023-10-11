// phylo intercepts
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  matrix[N, N] LRho;  // cholesky factor of known correlation matrix
}
parameters {
  vector[K] beta;  // regression coefficients
  real kappa;  // global int
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
  // optimized call for glm
  target += normal_id_glm_lupdf(Y | X, kappa + alpha, beta, sigma);

  // priors not including constants
  lprior += std_normal_lupdf(beta);
  lprior += normal_lupdf(kappa | 6, 3);
  lprior += std_normal_lupdf(sigma);
  lprior += std_normal_lupdf(sigma_phylo);
  target += lprior;
  target += std_normal_lupdf(z);
}
