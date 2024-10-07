graphics.off()
rm(list=ls())
cat("\014")  

library(tidyverse)
set.seed(0211)
this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
source("code/gmm_code.R")


for (i in 1:25)
{
  result_vec=c()
  for (j in 1:10)
  {
    data = readRDS(file = paste("data/dataset",i,"_",j,".RDS",sep = ""))
    X <- as.matrix(data)
    #ggplot() + 
    #  geom_point(data = data.frame(X), mapping = aes(X1, X2))
    K <- 2      # Number of clusters
    
    # Run vb-gmm model model
    vb_gmm_model <- vb_gmm(X = X, K = K, alpha_0 = 1e-5, max_iter = 2001, 
                           is_animation = FALSE, is_verbose = FALSE)
    
    pi = vb_gmm_model$pi_k
    mu = vb_gmm_model$m
    sigma= list()
    c_sigma = c()
    for (k in 1:length(pi))
    {
      sigma_k = vb_gmm_model$W[,,k] * (vb_gmm_model$nu[k]-1)*vb_gmm_model$beta[k]/(1+vb_gmm_model$beta[k])
      sigma[[k]] = solve(sigma_k)
      c_sigma_k = c(solve(sigma_k))
      c_sigma = cbind(c_sigma,c_sigma_k)
    }
    result_vec_i = rbind(mu,pi,c_sigma)
    result_vec = rbind(result_vec,result_vec_i)
  }
  #load the data

  #obj <- structure(list(mu=mu,sigma=sigma,pi=pi), class = "vb_gmm_post")
  #result[[i]] = obj
  
  write_csv(data.frame(result_vec),paste("result/gmm_result",i,".csv",sep = ""))
  
}

result_vec=c()
for (i in 1:25)
{
  data = readRDS(file = paste("data/dataset",i,".RDS",sep = ""))
  X <- as.matrix(data)
  #ggplot() + 
  #  geom_point(data = data.frame(X), mapping = aes(X1, X2))
  K <- 2      # Number of clusters
  
  # Run vb-gmm model model
  vb_gmm_model <- vb_gmm(X = X, K = K, alpha_0 = 1e-5, max_iter = 2001, 
                         is_animation = FALSE, is_verbose = FALSE)
  
  pi = vb_gmm_model$pi_k
  mu = vb_gmm_model$m
  sigma= list()
  c_sigma = c()
  for (k in 1:length(pi))
  {
    sigma_k = vb_gmm_model$W[,,k] * (vb_gmm_model$nu[k]-1)*vb_gmm_model$beta[k]/(1+vb_gmm_model$beta[k])
    sigma[[k]] = solve(sigma_k)
    c_sigma_k = c(solve(sigma_k))
    c_sigma = cbind(c_sigma,c_sigma_k)
  }
  result_vec_i = rbind(mu,pi,c_sigma)
  result_vec = rbind(result_vec,result_vec_i)
  

}
write_csv(data.frame(result_vec),paste("result/gmm_result.csv",sep = ""))

