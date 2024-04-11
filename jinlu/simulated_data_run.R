# Load data
load('./data/simulated/simulated_data.RData')
load('./data/simulated/mcmc_simulated.RData')
load('./data/simulated/psm_simulated.RData')
load('./data/simulated/simulated_Z.RData')
load('./data/simulated/mcmc_unique_simulated.RData')
load('./data/simulated/simulated_Z_reordered.RData')
load('./data/simulated/omega_JM_simulated.RData')

Y <- lapply(1:4,
            function(m) t(Y[[m]]))

# MCMC run
mcmc_simulated <- mcmc_run_all(Y = Y,
                              J = 10,
                              number_iter = 3000,
                              thinning = 5,
                              burn_in = 1000,
                              adaptive_prop = 0.001,
                              print_Z = TRUE,
                              
                              a_gamma = 100,
                              b_gamma = 2,
                              a_alpha = 1/5,
                              b_alpha = 1/2,
                              num.cores = 1)

psm_simulated <- similarity_matrix(mcmc_run_all_output = mcmc_simulated,
                                   num.cores = 1,
                                   run.on.pc = TRUE)



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
                                       burn_in = 1500,
                                       number_iter = 3000,
                                       Y = Y,
                                       a_gamma = 100,
                                       b_gamma = 2,
                                       regions.name = paste('region', 1:3))

mcmc_unique_simulated$estimated.pp.plot

simulated_Z_reordered <- mcmc_unique_simulated$Z

# w_jm
omega_JM_simulated <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_simulated,
                                        mcmc_run_post_output = mcmc_unique_simulated,
                                        thinning = 5,
                                        burn_in = 1500,
                                        number_iter = 3000)
                                   

omega_JM_simulated$omega_JM_plot


pp.standard.ordering(Y = Y,
                     Z = simulated_Z_reordered,
                     regions.name = paste('region', 1:3))


#---------------------- tv -----------------------

mytv_dist = function(x,ind){
  xdim = dim(mcmc_simulated$omega_J_M_output[[1]])
  y = matrix(x[,ind],xdim[1],xdim[2])
  return(0.5*colSums(abs(x-y)))
}



# Posterior mean
tv_mean = matrix(0, 4, 4)
for (m in c(1:4)){
  tv_dist_m = lapply(mcmc_simulated$omega_J_M_output,mytv_dist, ind = m )
  tv_dist_m = data.frame(matrix(unlist(tv_dist_m), nrow=length(tv_dist_m), byrow=TRUE))
  tv_mean[m,] = colMeans(tv_dist_m)
}
tv_mean = data.frame(tv_mean, row.names = c("1", "2", "3", "4"))
names(tv_mean) = c("1", "2", "3", "4")
print(tv_mean)

tv_mean = data.frame(tv_mean, "Mouse 1" = c("1", "2", "3", "4"))
names(tv_mean)[1:4] = c("1", "2", "3", "4")
tv_mean <-  pivot_longer(tv_mean,
                         cols = !Mouse.1,
                         names_to = "Mouse.2", 
                         values_to = "TV"
)

ggplot(tv_mean, aes(x = Mouse.1, y = Mouse.2, fill = TV)) +
  geom_tile() +
  labs(x = "Mouse",
       y = "Mouse") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.4) +
  geom_text(aes(label = round(tv_mean$TV,3)), color = "black", size = 4) +
  coord_fixed()
