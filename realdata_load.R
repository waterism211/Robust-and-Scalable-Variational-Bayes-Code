graphics.off()
rm(list=ls())
cat("\014")  

library(tidyverse)
set.seed(0211)
this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
source("code/gmm_code.R")

penguins = read.csv("data/penguins_size.csv")

penguins_new = na.omit(penguins)[1:300,c(1,3,4)]

max_value = apply(penguins_new[1:300,c(2,3)], MARGIN = 2, max)

for (i in 300:300)
{
  penguins_new[i,c(2,3)] = max_value*5
}

K <- 3    # Number of clusters
X = as.matrix(penguins_new[,c(2,3)])
# Run vb-gmm model model
vb_gmm_model <- vb_gmm(X = X, K = K, alpha_0 = 1e-5, max_iter = 10001, 
                       is_animation = FALSE, is_verbose = FALSE)


data.grid <- expand.grid(x = seq(from = min(X[,1]) - 2, 
                                 to = max(X[,1]) + 2, length.out = 100), 
                         y = seq(from = min(X[,2]) - 8, 
                                 to = max(X[,2]) + 2, length.out = 100))
q.samp <- cbind(data.grid, z = mixture_pdf_gaussian(vb_gmm_model,data.grid))



ggplot(data = penguins_new, mapping = aes(x= culmen_length_mm, y = culmen_depth_mm, color = species))+
  geom_contour(data = q.samp, mapping = aes(x = x,y = y, z = z, 
                                            colour = ..level..), binwidth = 0.001) + 
  gg_theme()

###split the data####
n = 300
perm <- sample(n)
# Create indices for ten folds
folds <- cut(seq(1, n), breaks=10, labels=FALSE)
# Split the dataframe based on the folds
split_penguins <- lapply(1:10, function(x) penguins_new[perm[folds == x], ])

ggplot(data = penguins_new, mapping = aes(x= culmen_length_mm, y = culmen_depth_mm, color = species))+
  geom_point()


result_vec=c()
for (i in 1:10)
{
  data = split_penguins[i][[1]][,2:3]
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

write_csv(data.frame(result_vec),paste("results_real_data/k3_penguins_result.csv",sep = ""))


result_vec=c()
for (i in 1:1)
{
  data = penguins_new[,2:3]
  X <- as.matrix(data)
  #ggplot() + 
  #  geom_point(data = data.frame(X), mapping = aes(X1, X2))
  K <- 3      # Number of clusters
  
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
write_csv(data.frame(result_vec),paste("results_real_data/k3_ori_penguins_result.csv",sep = ""))

