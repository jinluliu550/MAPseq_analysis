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

# Acceptance probabilities
plot(mcmc_simulated$acceptance_prob$omega, type = 'l')
plot(mcmc_simulated$acceptance_prob$alpha, type = 'l')
plot(mcmc_simulated$acceptance_prob$alpha_zero, type = 'l')
plot(mcmc_simulated$acceptance_prob$q_star, type = 'l')
plot(mcmc_simulated$acceptance_prob$gamma_star, type = 'l')
plot(mcmc_simulated$acceptance_prob$alpha_h, type = 'l')

# Trace plots
plot(mcmc_simulated$alpha_output, type = 'l')
plot(mcmc_simulated$alpha_zero_output , type = 'l')
alpha_h = matrix(unlist(mcmc_simulated$alpha_h_output),2001,R,byrow=TRUE)
plot(alpha_h[,1], type = 'l')
plot(alpha_h[,2], type = 'l')
plot(alpha_h[,3], type = 'l')

# PSM
psm_simulated <- similarity_matrix(mcmc_run_all_output = mcmc_simulated,
                                   num.cores = 1,
                                   run.on.pc = TRUE)
# Optimal clustering
simulated_Z <- opt.clustering.comb(mcmc_run_all_output = mcmc_simulated,
                                   post_similarity = psm_simulated)

# Posterior similarity matrix
plotpsm(psm.ind = psm_simulated$psm.within,
        psm.tot = psm_simulated$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)


plotpsm(psm.ind = psm_simulated$psm.within,
        psm.tot = psm_simulated$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')

# MCMC unique
mcmc_unique_simulated <- mcmc_run_post(mcmc_run_all_output = mcmc_simulated,
                                       Z = simulated_Z,
                                       thinning = 5,
                                       burn_in = 2000,
                                       number_iter = 12000,
                                       Y = Y,
                                       a_gamma = 100,
                                       b_gamma = 2,
                                       regions.name = paste('region', 1:mcmc_simulated$R))

mcmc_unique_simulated$estimated.pp.plot

simulated_Z_reordered <- mcmc_unique_simulated$Z

# w_jm
omega_JM_simulated <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_simulated,
                                        mcmc_run_post_output = mcmc_unique_simulated,
                                        thinning = 5,
                                        burn_in = 2000,
                                        number_iter = 12000)


omega_JM_simulated$omega_JM_plot


pp.standard.ordering(Y = Y,
                     Z = simulated_Z_reordered,
                     regions.name = paste('region', 1:3))




save.image("~/Documents/GitHub/MAPseq_analysis/data/simulated/mcmc_simulated_new.RData")