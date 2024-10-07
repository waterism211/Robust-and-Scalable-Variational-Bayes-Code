library(rstan)
library(tidyverse)
graphics.off()
rm(list=ls())
cat("\014")  

set.seed(0211)
this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
getwd()
gauss_mix<-stan_model(file = "code/gauss_mix.stan")

#generate the data
gauss_mix_simul<-function(n,mu,sd,pi)
{
  K<-length(mu)
  components <- sample(1:K,prob=pi,size=n,replace=TRUE)
  m_samples <- rnorm(n=n,mean=mu[components],sd=sd[components])
  return(m_samples)
}
gauss_data <- gauss_mix_simul(n=100,mu=c(-2,2),
                              sd=c(1,1),pi=c(0.5,0.5))

data = list(alpha = 0.1,a=0.1,b=0.1,N=100,K=2,X=gauss_data)

plot(density(gauss_data))

#for posterior predictive density
samples<-sampling(object = gauss_mix,data=data,iter=3000,warmup=1000,thin=1,
                  init="random",chains=1)
K=2
s<-summary(samples)
result<-s$summary[1:K,]

plot(samples)
stan_dens(samples,pars = "x_tilde")

#for missing data
gauss_data_missing = gauss_data
miss_idx = sample(1:100,size = 1)
gauss_data_missing[miss_idx] = 9999
data = list(alpha = 0.1,a=0.1,b=0.1,N=100,K=2,miss_idx = miss_idx,
            X_obs=gauss_data_missing)
gauss_mix_missing<-stan_model(file = "code/gauss_mix_missing.stan")


data = list(alpha = 0.1,a=0.1,b=0.1,N=100,K=2,X=gauss_data_missing)
#for posterior predictive density
samples<-sampling(object = gauss_mix_missing,data=data,iter=3000,
                  warmup=1000,thin=1,init="random",chains=1)
stan_dens(samples,pars = "X_miss")
