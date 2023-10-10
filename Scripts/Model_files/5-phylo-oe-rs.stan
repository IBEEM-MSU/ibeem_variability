// phylo intercepts w/ reduce_sum to parallelize within chains

functions {
  /* integer sequence of values
   * Args:
   *   start: starting integer
   *   end: ending integer
   * Returns:
   *   an integer sequence from start to end
   */
  array[] int sequence(int start, int end) {
    array[end - start + 1] int seq;
    for (n in 1 : num_elements(seq)) {
      seq[n] = n + start - 1;
    }
    return seq;
  }
  // compute partial sums of the log-likelihood
  real partial_log_lik_lpmf(array[] int seq, int start, int end,
                            data vector y_obs, data matrix X_obs, vector beta,
                            real kappa, vector alpha, real sigma) {
    real ptarget = 0;
    int N = end - start + 1;

    ptarget += normal_id_glm_lupdf(y_obs[start : end] | X_obs[start : end], kappa + alpha[start : end], beta, sigma);
    return ptarget;
  }
}

data {
  int<lower=1> N;  // total number of observations
  // int<lower=1> No;
  // int<lower=1> Ni;
  vector[N] y_obs;  // response variable
  // vector[Ni] y_imp;  // response variable
  // vector[Ni] sd_y;
  int grainsize; // grainsize for threading
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X_obs;  // population-level design matrix
  // matrix[Ni, K] X_imp;  // population-level design matrix
  // int<lower=1> obs_idx[No]; // idx of alpha that corresponds to obs value
  // int<lower=1> imp_idx[Ni]; // idx of alpha that corresponds to imp value
  matrix[N, N] LRho;  // cholesky factor of known correlation matrix
}
transformed data {
  array[N] int seq = sequence(1, N);
}
parameters {
  vector[K] beta;  // regression coefficients
  real kappa;  // global int
  real<lower=0> sigma;  // process error
  real<lower=0> sigma_phylo;  // phylo sd
  vector[N] z;  // standardized group-level effects
  // vector[Ni] y_iv_raw;
}
transformed parameters {
  vector[N] alpha;  // phylo intercepts (per species)
  
  // implies alpha ~ MVN(0, Rho) * sigma_phylo
  alpha = (sigma_phylo * (LRho * z));
}
model {
  real lprior = 0;  // prior contributions to the log posterior

  // optimized call for glm
  target += reduce_sum(partial_log_lik_lpmf, seq, grainsize, y_obs, X_obs, beta, 
  kappa, alpha, sigma);

  // priors not including constants
  lprior += std_normal_lupdf(beta);
  lprior += normal_lupdf(kappa | 5, 2);
  lprior += std_normal_lupdf(sigma);
  lprior += std_normal_lupdf(sigma_phylo);
  // lprior += normal_lupdf(sigma_phylo | 0, 2);
  target += lprior;
  target += std_normal_lupdf(z);
}
