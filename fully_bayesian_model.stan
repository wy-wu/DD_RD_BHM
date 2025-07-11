data {
  int<lower=1> N;           // total number of observations
  int<lower=1> K;           // number of groups
  int<lower=1, upper=K> group[N]; // group indicator for each observation
  vector[N] x;              // observed data

  real mu0;                 // hyperprior mean for u
  real<lower=0> tau0;       // hyperprior std for u
  real<lower=0> eta;        // scale for half-Cauchy prior on r
  real<lower=0> xi;         // scale for half-Cauchy prior on s
}

parameters {
  real<lower=0> s;          // std dev of likelihood
  real<lower=0> r;          // std dev of m_i prior
  real u;                   // global mean
  vector[K] m;              // group-specific means
}

model {
  // Hyperpriors
  u ~ normal(mu0, tau0);
  r ~ cauchy(0, eta);
  s ~ cauchy(0, xi);

  // Priors for group means
  m ~ normal(u, r);

  // Likelihood
  for (n in 1:N)
    x[n] ~ normal(m[group[n]], s);
}
