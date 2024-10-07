data {
  //Hyperparameters
  real a;
  real b;
  //real s;
  real alpha;
  //define the data
  int<lower=1> N; // sapmle size
  int<lower=1> K; // of clusters
  vector[N] X_obs;
  int<lower=1> miss_idx;
}

parameters {
  vector[K] mu;
  vector<lower=0>[K] sigma;
  simplex[K] pi;
  real X_miss;
}

transformed parameters {
  vector<lower=0>[K] sigma2;
  vector[N] X;
  sigma2=square(sigma);
  X = X_obs;
  X[miss_idx] = X_miss;
}


model {
//temporary vector for loop
real contributions[K];
//prior
//mean
mu ~ normal(0,sigma);
//variance
sigma2 ~ inv_gamma(a,b);
//cluster probability
pi ~ dirichlet(rep_vector(alpha,K));



//likelihood
for (i in 1:N){
  for (k in 1:K){
      contributions[k] = log(pi[k])+normal_lpdf(X[i] | mu[k],sigma[k]);
    }
    target += log_sum_exp(contributions);
  }
}
