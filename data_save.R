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
  data = readRDS(file = paste("data/dataset",i,".RDS",sep = ""))
  write_csv(data,paste("data_cvs/dataset",i,".csv",sep = ""))
}

