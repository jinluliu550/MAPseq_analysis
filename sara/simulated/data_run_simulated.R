#Fit on simulated data


setwd("~/Documents/GitHub/MAPseq_analysis")
load("~/data/simulated/simulated_data.RData")
source("~/R/main.R")

Y=lapply(Y,t)

mcmc_simulated = mcmc_run_all(Y,
                         J=10,
                         number_iter = 12000,
                         thinning = 5,
                         burn_in = 2000,
                         adaptive_prop = 0.001,
                         print_Z = TRUE,
                         iter_update = 100,
                         a_gamma=100,
                         b_gamma=2,
                         a_alpha = 1/5,
                         b_alpha = 1/2,
                         num.cores = 1)
