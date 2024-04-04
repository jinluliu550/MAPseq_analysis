Hans_data1 <- read.csv('./data/Han-data/han_brain4_gh.csv')
Hans_data2 <- read.csv('./data/Han-data/han_brain5_gh.csv')
Hans_data3 <- read.csv('./data/Han-data/han_brain6_gh.csv')

data_Hans <- list(t(Hans_data1),
                  t(Hans_data2),
                  t(Hans_data3))

load('./data/Han-data/mcmc_all_hans.RData')
load('./data/Han-data/psm_hans.RData')
load('./data/Han-data/mcmc_unique_Han.RData')

mcmc_all_hans <- mcmc_run_all(Y = data_Hans,
                             J = 20,
                             number_iter = 8000,
                             thinning = 5,
                             burn_in = 3000,
                             adaptive_prop = 0.1,
                             print_Z = TRUE,
                             
                             
                             a_gamma = 300,
                             b_gamma = 10,
                             a_alpha = 1/5,
                             b_alpha = 1/2,
                             num.cores = 10)



psm_hans <- similarity_matrix(mcmc_run_all_output = mcmc_all_hans,
                              num.cores = 10,
                              run.on.pc = FALSE)

# Optimal clustering
han_Z <- opt.clustering.comb(mcmc_run_all_output = mcmc_all_hans,
                             post_similarity = psm_hans)


# Plot of posterior similarity matrix
plotpsm(psm.ind = psm_hans$psm.within,
        psm.tot = psm_hans$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)


plotpsm(psm.ind = psm_hans$psm.within,
        psm.tot = psm_hans$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')

# MCMC unique
mcmc_unique_Han <- mcmc_run_post(mcmc_run_all_output = mcmc_all_hans,
                                 Z = han_Z,
                                 thinning = 5,
                                 burn_in = 1500,
                                 number_iter = 3000,
                                 Y = data_Hans,
                                 a_gamma = 500,
                                 b_gamma = 10,
                                 regions.name = rownames(data_Hans[[1]]))

# Number of neurons in each cluster
png(file = './plots/Hans/number_of_neuron.png',
    width = 1000,
    height = 700)

opt.clustering.frequency(clustering = mcmc_unique_Han$Z)

dev.off()

# Plot of estimated projection strength
png(file = './plots/Hans/estimated_pp.png',
    width = 1500,
    height = 800)

mcmc_unique_Han$estimated.pp.plot

dev.off()

# q_tilde
png(file = './plots/Hans/q_tilde.png',
    width = 1000,
    height = 300)

mcmc_unique_Han$q_tilde_plot

dev.off()



# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_hans,
                                   mcmc_run_post_output = mcmc_unique_Han,
                                   thinning = 5,
                                   burn_in = 1500,
                                   number_iter = 3000)

png(file = './plots/Hans/w_jm_EC.png',
    width = 1500,
    height = 800)

omega_JM_mcmc$omega_JM_plot

dev.off()


# Examining the difference of omega_jm between any 2 mice:
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc)


png(file = './plots/Hans/w_jm_difference.png',
    width = 1000,
    height = 300)

difference_omega_JM$probability_plot

dev.off()


# Heatmap of projection stength of neurons in each cluster
png(file = './plots/Hans/heatmap_neuron.png',
    width = 600,
    height = 600)

pp.standard.ordering(Y = data_Hans,
                     Z = mcmc_unique_Han$Z,
                     regions.name = rownames(data_Hans[[1]]))

dev.off()


# Distribution of N_{i,m} for neurons within the same cluster
data_EC_N <- lapply(1:length(data_Hans),
                    function(m) colSums(data_Hans[[m]]))

df <- data.frame(N = unlist(data_EC_N),
                 motif = unlist(mcmc_unique_Han$Z))

png(file = './plots/Hans/N_sum_by_cluster.png',
    width = 1000,
    height = 400)

ggplot(df)+
  geom_boxplot(mapping = aes(x = factor(motif, levels = 1:max(unlist(mcmc_unique_Han$Z))),
                             y = N))+
  theme_bw()+
  xlab('Cluster')

dev.off()


# Posterior predictive check with multiple replicates
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_hans,
                      Y = data_Hans,
                      N = 3,
                      regions.name = rownames(data_Hans[[1]]))


png(file = './plots/Hans/ppc_zero.png',
    width = 1200,
    height = 700)

ppc_multiple$zero.plot

dev.off()

png(file = './plots/Hans/ppc_nonzero.png',
    width = 1200,
    height = 700)

ppc_multiple$non.zero.plot

dev.off()


# Posterior predictive check with single replicated data
ppc_single <- ppc_single_f(mcmc_run_all_output = mcmc_all_hans,
                           Y = data_Hans,
                           regions.name = rownames(data_Hans[[1]]))

for(m in 1:length(data_Hans)){
  
  png(file = paste0('./plots/Hans/ppc_single_mouse_', m, '.png'),
      width = 1200,
      height = 700)
  
  print(ppc_single[[m]])
  
  dev.off()
}
