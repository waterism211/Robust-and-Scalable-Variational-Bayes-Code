library("gtools");  # for rdirichlet
library("rstan")
library(parallel)
options(mc.cores = detectCores() - 1)
library(ggplot2)
library(tidyverse)
library(SBmedian)
library(dplyr)
library(tidyr)
graphics.off()
rm(list=ls())
cat("\014")  
set.seed(0211)
this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
source("func.R")

# Set parameters
V <- 4  # Vocabulary size
K <- 2  # Number of topics
M <- 19 # Number of documents
alpha_val <- 2
beta_val <- 1

# Range of N1 values
N1_values <- c(5, 10, 20, 30,40, 50,70, 100,200)

# Initialize a data frame to store summary statistics for plotting
# Initialize data frames to store mean and timing information for plotting
mean_data <- data.frame()
time_data <- data.frame()

# Iterate over N1 values
for (N1 in N1_values) {
  cat("\nRunning simulation and model fitting for N1 =", N1, "\n")
  
  # Simulate LDA data
  result <- simulate_lda_data(V = V, K = K, M = M, N1 = N1, avg_doc_length = 10, alpha_val = alpha_val, beta_val = beta_val)
  
  # Time Bayesian fitting and calculate mean KL
  bayes_time <- system.time({
    fitlda <- stan(file = 'lda.stan', data = result, iter = 1000, chains = 2)
    klphifit <- rstan::extract(fitlda)$klphi
  })
  klphifit_mean <- mean(as.vector(klphifit))
  
  # Time Variational Bayes fitting and calculate mean KL
  vb_time <- system.time({
    ldam <- stan_model(file = 'lda.stan', verbose = TRUE)
    ldafit_vb <- vb(ldam, data = result, tol_rel_obj = 1e-3)
    klphifit_vb <- rstan::extract(ldafit_vb)$klphi
  })
  klphifit_vb_mean <- mean(as.vector(klphifit_vb))
  
  # Split data for geometric median calculation
  split_data <- split_lda_data(result, n_splits = 10)
  
  # Time Geometric Median for Bayesian KL and calculate mean
  geo_bayes_time <- system.time({
    klphifits <- extract_klphifits_bayes(split_data = split_data, result = result, stan_file = 'lda.stan', iter = 2000, chains = 2)
    d2 <- dim(klphifits[[1]])[2]
    d3 <- dim(klphifits[[1]])[3]
    kl_list <- vector("list", d2 * d3)
    pos <- 1
    for (i in 1:d2) {
      for (j in 1:d3) {
        mydata <- lapply(1:10, function(k) t(klphifits[[k]][, i, j]))
        myrun <- mpost.euc(mydata, show.progress = FALSE)
        id <- sample(1:nrow(myrun$med.atoms), 100, prob = myrun$med.weights, replace = TRUE)
        kl_list[[pos]] <- as.vector(myrun$med.atoms[id, ])
        pos <- pos + 1
      }
    }
  })
  kl_list_mean <- mean(unlist(kl_list))
  
  # Time Geometric Median for VB KL and calculate mean
  geo_vb_time <- system.time({
    klphifits_vb <- extract_klphifits_vb(split_data = split_data, result = result, stan_file = 'lda.stan', tol_rel_obj = 1e-3)
    kl_vb_list <- vector("list", d2 * d3)
    pos <- 1
    for (i in 1:d2) {
      for (j in 1:d3) {
        mydata <- lapply(1:10, function(k) t(klphifits_vb[[k]][, i, j]))
        myrun <- mpost.euc(mydata, show.progress = FALSE)
        id <- sample(1:nrow(myrun$med.atoms), 100, prob = myrun$med.weights, replace = TRUE)
        kl_vb_list[[pos]] <- as.vector(myrun$med.atoms[id, ])
        pos <- pos + 1
      }
    }
  })
  kl_vb_list_mean <- mean(unlist(kl_vb_list))
  
  # Store mean KL divergence data for plotting
  mean_data <- mean_data %>%
    bind_rows(data.frame(
      N1 = N1,
      Method = 'Bayesian KL',
      Mean = klphifit_mean
    )) %>%
    bind_rows(data.frame(
      N1 = N1,
      Method = 'VB KL',
      Mean = klphifit_vb_mean
    )) %>%
    bind_rows(data.frame(
      N1 = N1,
      Method = 'Median Bayesian KL',
      Mean = kl_list_mean
    )) %>%
    bind_rows(data.frame(
      N1 = N1,
      Method = 'Median VB KL',
      Mean = kl_vb_list_mean
    ))
  
  # Store timing data for plotting
  time_data <- time_data %>%
    bind_rows(data.frame(
      N1 = N1,
      Method = 'Bayesian KL',
      Time = bayes_time['elapsed']
    )) %>%
    bind_rows(data.frame(
      N1 = N1,
      Method = 'VB KL',
      Time = vb_time['elapsed']
    )) %>%
    bind_rows(data.frame(
      N1 = N1,
      Method = 'Median Bayesian KL',
      Time = geo_bayes_time['elapsed']
    )) %>%
    bind_rows(data.frame(
      N1 = N1,
      Method = 'Median VB KL',
      Time = geo_vb_time['elapsed']
    ))
}

# Plot the mean values
ggplot(mean_data, aes(x = factor(N1), y = Mean, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(x = "# of words in the outlier document",
       y = "Mean KL Divergence",
       color = "Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(time_data, aes(x = factor(N1), y = Time, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(x = "# of words in the outlier document",
       y = "Execution Time (seconds)",
       color = "Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
