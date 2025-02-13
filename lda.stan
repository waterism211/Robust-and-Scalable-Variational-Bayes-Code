data {
  int<lower=2> K;               // num topics
  int<lower=2> V;               // num words
  int<lower=1> M;               // num docs
  int<lower=1> N;               // total word instances
  int<lower=1,upper=V> w[N];    // word n
  int<lower=1,upper=M> doc[N];  // doc ID for word n
  vector<lower=0>[K] alpha;     // topic prior
  vector<lower=0>[V] beta;      // word prior
  simplex[V] truephi0[K];
}
parameters {
  simplex[K] theta[M];   // topic dist for doc m
  simplex[V] phi[K];     // word dist for topic k
}
model {
  for (m in 1:M)  
    theta[m] ~ dirichlet(alpha);  // prior
  for (k in 1:K)  
    phi[k] ~ dirichlet(beta);     // prior
  for (n in 1:N) {
    real gamma[K];
    for (k in 1:K) 
      gamma[k] = log(theta[doc[n],k]) + log(phi[k,w[n]]);
    increment_log_prob(log_sum_exp(gamma));  // likelihood
  }
}
generated quantities {
  vector<lower=0>[K] klphi[K];
  real temp;
  temp = 0;
  for (i in 1:K){
    for (j in 1:K){
      temp = 0;
      for (k in 1:V){
        temp = temp + truephi0[j][k] * log(truephi0[j][k]) - truephi0[j][k] * log(phi[i][k]);
      }
      klphi[i][j] = temp;
    }
  }
}