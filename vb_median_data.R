graphics.off()
rm(list=ls())
cat("\014")  

library(tidyverse)
set.seed(0211)
this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
source("code/gmm_code.R")

n = 200
mu = list(c(2, 2), c(-2, -2))
sigma = list(diag(1, 2), diag(1, 2))
pi = c(0.5, 0.5)

for (i in 1:25)
{
  simulated_data = gauss_mix_simul(n=n,mu=mu,sigma = sigma, pi=pi)
  simulated_data[200,] = i*apply(simulated_data,2,max)
  names(simulated_data) = c("X1","X2")
  simulated_data = data.frame(simulated_data)
  saveRDS(simulated_data,file = paste("data/dataset",i,".RDS", sep=""))
  perm <- sample(n)
  # Create indices for ten folds
  folds <- cut(seq(1, n), breaks=10, labels=FALSE)
  # Split the dataframe based on the folds
  split_df <- lapply(1:10, function(x) simulated_data[perm[folds == x], ])
  for (j in 1:10)
  {
    saveRDS(split_df[[j]],file = paste("data/dataset",i,"_",j,".RDS", sep=""))
  }
}


