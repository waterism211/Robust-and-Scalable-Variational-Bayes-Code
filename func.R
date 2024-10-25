library("gtools");  # for rdirichlet
library("rstan")
library(parallel)
options(mc.cores = detectCores() - 1)
library(ggplot2)
library(tidyverse)
library(SBmedian)
graphics.off()
rm(list=ls())
cat("\014")  

set.seed(0211)
this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)


simulate_lda_data <- function(V, K, M, N1, avg_doc_length = 10, alpha_val = 2, beta_val = 1) {
  # Set priors
  alpha <- rep(alpha_val, K)  # Dirichlet prior for topic distribution per document
  beta <- rep(beta_val, V)    # Dirichlet prior for word distribution per topic
  
  # Generate the word distribution per topic (phi)
  phi <- rdirichlet(K, beta)
  truephi0 <- phi
  
  # Generate the topic distribution per document (theta)
  theta <- rdirichlet(M, alpha)
  truetheta0 <- theta
  
  # Generate document lengths
  doc_length <- rpois(M, avg_doc_length)
  N <- sum(doc_length)  # Total number of words
  
  # Initialize vectors to store words and document indices
  w <- rep(NA, N)
  doc <- rep(NA, N)
  
  # Generate words for each document
  n <- 1
  for (m in 1:M) {
    for (i in 1:doc_length[m]) {
      z <- which(rmultinom(1, 1, theta[m,]) == 1)  # Sample topic for the word
      w[n] <- which(rmultinom(1, 1, phi[z,]) == 1)  # Sample word based on the topic
      doc[n] <- m  # Assign document index
      n <- n + 1
    }
  }
  
  # Add an additional document with N1 words all being the first word in the vocabulary
  doc_length <- c(doc_length, N1)
  M <- M + 1
  for (i in 1:N1) {
    w[n] <- 1  # All words are the first word
    doc[n] <- M  # Assign the additional document index
    n <- n + 1
  }
  
  # Update total number of words after adding the new document
  N <- N + N1
  
  # Return a list containing the generated data
  return(list(
    w = w,  # Word indices
    doc = doc,  # Document indices
    theta = theta,  # Topic distribution per document
    phi = phi,  # Word distribution per topic
    doc_length = doc_length,  # Document lengths
    M = M,  # Total number of documents
    N = N,  # Total number of words
    alpha = alpha,
    beta = beta,
    K = K,
    V = V,
    truephi0=truephi0
  ))
}



split_lda_data <- function(result, n_splits = 10) {
  # Randomly assign documents to n_splits groups
  w = result$w
  doc = result$doc
  doc_length = result$doc_length
  M = result$M
  random_split <- sample(cut(1:M, breaks = n_splits, labels = FALSE))
  
  # Initialize lists to hold the split data
  split_docs_id <- vector("list", n_splits)
  split_w <- vector("list", n_splits)
  split_doc_lengths <- vector("list", n_splits)
  split_doc_index <- vector("list", n_splits)
  
  # Document IDs
  doc_id <- seq(1:M)
  
  # Loop over the splits and extract the relevant data
  for (i in 1:n_splits) {
    # Find the documents assigned to this split
    subset_indices <- random_split == i
    split_docs_id[[i]] <- doc_id[subset_indices]
    
    # Get the document IDs for the current group
    doc_ids_in_group <- split_docs_id[[i]]
    
    # Find the corresponding indices in the `doc` array for these document IDs
    indices_in_group <- which(doc %in% doc_ids_in_group)
    
    # Extract the words and document lengths for the current group
    split_w[[i]] <- w[indices_in_group]
    split_doc_lengths[[i]] <- doc_length[doc_ids_in_group]
    
    # Create the re-indexed document index for Stan
    split_doc_index[[i]] <- integer()
    for (j in 1:length(split_doc_lengths[[i]])) {
      split_doc_index[[i]] <- c(split_doc_index[[i]], rep(j, split_doc_lengths[[i]][j]))
    }
  }
  
  # Return the split data as a list
  return(list(
    split_w = split_w,  # Words for each split
    split_doc_lengths = split_doc_lengths,  # Document lengths for each split
    split_doc_index = split_doc_index,  # Re-indexed document indices
    split_docs_id = split_docs_id  # Document IDs for each split
  ))
}


extract_klphifits_bayes <- function(split_data, result, stan_file, iter = 1000, chains = 2) {
  # Initialize a list to hold the KL divergence for each subset
  K = result$K
  V = result$V
  alpha = result$alpha
  beta = result$beta
  truephi0 = result$truephi0
  klphifits <- list()
  
  # Loop over the splits
  for (i in 1:length(split_data$split_w)) {
    # Prepare the subset data for Stan
    subset_data <- list(
      K = K,  # Number of topics
      V = V,  # Vocabulary size
      M = length(split_data$split_doc_lengths[[i]]),  # Number of documents in the current subset
      N = sum(split_data$split_doc_lengths[[i]]),  # Total number of words in the current subset
      w = split_data$split_w[[i]],  # Word indices for the current subset
      doc = split_data$split_doc_index[[i]],  # Document indices for the current subset (re-indexed)
      alpha = alpha,  # Dirichlet prior for topic distribution per document
      beta = beta,  # Dirichlet prior for word distribution per topic
      truephi0 = truephi0  # True word distribution for comparison (if needed)
    )
    
    # Fit the LDA model using Stan
    fitlda_subset <- stan(file = stan_file, data = subset_data, iter = iter, chains = chains)
    
    # Extract KL divergences (klphi) from the Stan output
    klphifit_subset <- rstan::extract(fitlda_subset)$klphi
    
    # Store the KL divergences for this subset
    klphifits[[i]] <- klphifit_subset
  }
  
  # Return the list of KL divergence values for each subset
  return(klphifits)
}

extract_klphifits_vb <- function(split_data, result, stan_file, tol_rel_obj = 1e-3) {
  # Initialize a list to hold the KL divergence for each subset
  ldam = stan_model(file = 'lda.stan', verbose = TRUE)
  K = result$K
  V = result$V
  alpha = result$alpha
  beta = result$beta
  truephi0 = result$truephi0
  klphifits_vb <- list()
  
  # Loop over the splits
  for (i in 1:length(split_data$split_w)) {
    # Prepare the subset data for Stan
    subset_data <- list(
      K = K,  # Number of topics
      V = V,  # Vocabulary size
      M = length(split_data$split_doc_lengths[[i]]),  # Number of documents in the current subset
      N = sum(split_data$split_doc_lengths[[i]]),  # Total number of words in the current subset
      w = split_data$split_w[[i]],  # Word indices for the current subset
      doc = split_data$split_doc_index[[i]],  # Document indices for the current subset (re-indexed)
      alpha = alpha,  # Dirichlet prior for topic distribution per document
      beta = beta,  # Dirichlet prior for word distribution per topic
      truephi0 = truephi0  # True word distribution for comparison (if needed)
    )
    
    # Fit the LDA model using Stan
    ldafit_subset = vb(ldam, data = subset_data, tol_rel_obj = 1e-3)
    
    # Extract KL divergences (klphi) from the Stan output
    klphifit_vb_subset <- rstan::extract(ldafit_subset)$klphi
    
    # Store the KL divergences for this subset
    klphifits_vb[[i]] <- klphifit_vb_subset
  }
  
  # Return the list of KL divergence values for each subset
  return(klphifits)
}


