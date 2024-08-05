
# Hans data - also called the VC data
# Make sure data_Hans_5 and data_Hans are loaded before running the code below

M <- length(data_Hans)
C <- sapply(1:M, function(m) ncol(data_Hans[[m]]))


# Gel-plot
gel_plot_data_Hans <- gel_plot(Y = data_Hans)


png(file = './plots/Hans/gelplot_hans.png',
    width = 2000,
    height = 500)

ggarrange(gel_plot_data_Hans[[1]],
          gel_plot_data_Hans[[2]],
          gel_plot_data_Hans[[3]],
          gel_plot_data_Hans[[4]],
          nrow = 1,
          widths = c(1,1,1,1.3))

dev.off()


# Mouse index
mouse.index <- c(rep(1, C[1]),
                 rep(2, C[2]),
                 rep(3, C[3]),
                 rep(4, C[4]))

#---------------------------------------------------------------------------------------

data_Hans_cbind <- do.call(cbind, data_Hans)

df <- t(data_Hans_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))


C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(data_Hans[[m]]))))

# k-means
k_mean_clust_20 <- kmeans(df, 20, iter.max = 100, nstart = 25)$cluster

clust20 <- lapply(1:M,
                  function(m) k_mean_clust_20[(C_cumsum[m]+1):C_cumsum[m+1]])

# Run without threshold
mcmc_all_hans <- mcmc_run_all(Y = data_Hans,
                              J = 40,
                              number_iter = 20000,
                              thinning = 5,
                              burn_in = 5000,
                              adaptive_prop = 0.0001,
                              print_Z = TRUE,
                              a_gamma = 20,
                              b_gamma = 1,
                              a_alpha = 1/5,
                              b_alpha = 1/2,
                              Z.init = clust20)

#----------------------------------------------------------------------------------------
                              

Zmat = matrix(unlist(mcmc_all_hans$Z_output), length(mcmc_all_hans$Z_output), sum(C),byrow = TRUE)

# Number of occupied components
k = apply(Zmat,1,function(x){length(unique(x))})
plot(k, type = 'l')


# Posterior similarity matrix
psm_hans = similarity_matrix(mcmc_run_all_output = mcmc_all_hans)


# Reordered posterior samples of z
hans_z_reordered <- z_trace_updated(mcmc_run_all_output = mcmc_all_hans)


# optimal clustering
hans_Z <- opt.clustering.comb(z_trace = hans_z_reordered,
                              post_similarity = psm_hans,
                              max.k = max(k))


#-- Convert to a list
C_cumsum <- c(0, cumsum(C))

hans_Z <- lapply(1:M,
                function(m) hans_Z[(C_cumsum[m]+1):C_cumsum[m+1]])

#-------------------------------------------------------------------------------------------


# Plot of posterior similarity matrix
psm_hans_plot <- plotpsm(psm.ind = psm_hans$psm.within,
                         psm.tot = psm_hans$psm.combined)


png(file = './plots/Hans/heatmap_psm(1).png',
    width = 500,
    height = 400)

psm_hans_plot$plot.ind

dev.off()


png(file = './plots/Hans/heatmap_psm(2).png',
    width = 500,
    height = 400)

psm_hans_plot$plot.tot

dev.off()


# MCMC unique
mcmc_unique_hans <- mcmc_run_post(mcmc_run_all_output = mcmc_all_hans,
                                  Z = hans_Z,
                                  thinning = 5,
                                  burn_in = 2000,
                                  number_iter = 12000,
                                  Y = data_Hans,
                                  a_gamma = 20,
                                  b_gamma = 1,
                                  regions.name = rownames(data_Hans[[1]]))

                                 
# Neuron allocations
hans_Z_reordered <- mcmc_unique_hans$Z


# Clusters with at least 20 neurons
large.bayesian.cluster <- which(as.vector(table(unlist(hans_Z_reordered))) >= 20)




# Number of neurons by cluster and mosue
png(file = './plots/Hans/number_of_neuron_by_m.png',
    width = 600,
    height = 300)

opt.clustering.frequency(clustering = mcmc_unique_hans$Z)

dev.off()


# Plot of estimated projection strength
png(file = './plots/Hans/estimated_pp.png',
    width = 1200,
    height = 800)

mcmc_unique_hans$estimated.pp.plot

dev.off()



# q tilde
png(file = './plots/Hans/q_tilde.png',
    width = 600,
    height = 250)

mcmc_unique_hans$q_tilde_plot

dev.off()



# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_hans,
                                   mcmc_run_post_output = mcmc_unique_hans,
                                   thinning = 5,
                                   burn_in = 2000,
                                   number_iter = 12000)


png(file = './plots/Hans/w_jm.png',
    width = 1200,
    height = 800)

omega_JM_mcmc$omega_JM_plot

dev.off()



# Examining the difference of omega_jm between any 2 mice:
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc)


png(file = './plots/Hans/w_jm_difference.png',
    width = 800,
    height = 300)

difference_omega_JM$probability_plot

dev.off()


# data 1
significant_obs <- difference_omega_JM$significant_obs

data1_significant <- significant_obs %>%
  filter((data1 == 1) | (data2 == 1)) %>%
  pull(cluster)

length(data1_significant)

# data 2
data2_significant <- significant_obs %>%
  filter(data1 == 2 | data2 == 2) %>%
  pull(cluster)

length(data2_significant)

# data 3
data3_significant <- significant_obs %>%
  filter(data1 == 3 | data2 == 3) %>%
  pull(cluster)

length(data3_significant)

# data 4
data4_significant <- significant_obs %>%
  filter(data1 == 4 | data2 == 4) %>%
  pull(cluster)

length(data4_significant)

# Plot to show examples of differences
png(file = './plots/Hans/w_jm_difference_eg.png',
    width = 800,
    height = 400)

difference_in_omega_jm_plot(difference_in_omega_jm_output = difference_omega_JM,
                            mcmc_run_omega_output = omega_JM_mcmc,
                            N = 6)

dev.off()


# Heatmap of projection strength of neurons in each cluster
png(file = './plots/Hans/heatmap_neuron.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_Hans,
                      Z = mcmc_unique_hans$Z,
                      regions.name = rownames(data_Hans[[1]]),
                      mouse.index = mouse.index)

dev.off()


# Distribution of N_{i,m} for neurons within the same cluster

data_Hans_N <- lapply(1:length(data_Hans),
                      function(m) colSums(data_Hans[[m]]))

df <- data.frame(N = unlist(data_Hans_N),
                 motif = unlist(mcmc_unique_hans$Z))



png(file = './plots/Hans/N_sum_by_cluster.png',
    width = 800,
    height = 400)

ggplot(df)+
  geom_boxplot(mapping = aes(x = factor(motif, levels = 1:max(unlist(mcmc_unique_hans$Z))),
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
    width = 600,
    height = 300)

ppc_multiple$zero.plot

dev.off()


png(file = './plots/Hans/ppc_nonzero.png',
    width = 600,
    height = 300)

ppc_multiple$non.zero.plot

dev.off()



# Posterior predictive check with single replicated data
ppc_single <- ppc_single_f(mcmc_run_all_output = mcmc_all_hans,
                           Y = data_Hans,
                           regions.name = rownames(data_Hans[[1]]))


for(m in 1:M){
  
  png(file = paste0('./plots/Hans/ppc_single_mouse_', m, '.png'),
      width = 800,
      height = 400)
  
  print(ppc_single[[m]])
  
  dev.off()
}




#--------------------------------------------------- K-means -----------------------------

df <- t(apply(t(data_Hans_cbind), 1, function(x){return(x/sum(x))}))

# k = 20
clust20 <- lapply(1:M,
                  function(m) kmeans(df, 20, nstart = 25)$cluster[(C_cumsum[m]+1):C_cumsum[m+1]])

clust20_r <- k_means_reorder(Y = data_Hans,
                             Z = clust20)


# k = 30
clust30 <- lapply(1:M,
                  function(m) kmeans(df, 30, nstart = 25)$cluster[(C_cumsum[m]+1):C_cumsum[m+1]])

clust30_r <- k_means_reorder(Y = data_Hans,
                             Z = clust30)


# k = 40
clust40 <- lapply(1:M,
                  function(m) kmeans(df, 40, nstart = 25)$cluster[(C_cumsum[m]+1):C_cumsum[m+1]])


clust40_r <- k_means_reorder(Y = data_Hans,
                             Z = clust40)




# Heatmap of projection strength of neurons in each cluster
png(file = './plots/Hans/k_means_20.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_Hans,
                      Z = clust20_r,
                      regions.name = rownames(data_Hans[[1]]),
                      mouse.index = mouse.index)

dev.off()

png(file = './plots/Hans/k_means_30.png',
    width = 600,
    height = 600)


pp.standard.ordering2(Y = data_Hans,
                      Z = clust30_r,
                      regions.name = rownames(data_Hans[[1]]),
                      mouse.index = mouse.index)

dev.off()

png(file = './plots/Hans/k_means_40.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_Hans,
                     Z = clust40_r,
                     regions.name = rownames(data_Hans[[1]]),
                     mouse.index = mouse.index)

dev.off()



#--------------------------- Binomial clustering ----------------------------------

#------------------ Binomial Model 2 - Han's paper

# Assume less than 5 are noise
hans_binomial_5 <- binomial_model2(data = data_Hans,
                                   threshold = 5)


hans_binomial_reorder5 <- binom_cluster_reorder(Y = data_Hans,
                                                binomial_output = hans_binomial_5)




# Assume no noise
hans_binomial_1 <- binomial_model2(data = data_Hans,
                                   threshold = 1)


hans_binomial_reorder1 <- binom_cluster_reorder(Y = data_Hans,
                                                binomial_output = hans_binomial_1)


#----------------- Binomial Model 1 - Hippocampus paper

hans_binomial_5_pe <- binomial_model(data = data_Hans,
                                     threshold = 5)

hans_binomial_reorder5_pe <- binom_cluster_reorder(Y = data_Hans,
                                                   binomial_output = hans_binomial_5_pe)

hans_binomial_1_pe <- binomial_model(data = data_Hans,
                                     threshold = 1)

hans_binomial_reorder1_pe <- binom_cluster_reorder(Y = data_Hans,
                                                   binomial_output = hans_binomial_1_pe)



# hans_binomial_reorder5$cluster_summary %>%
#   filter(significant == 'significant',
#          cluster.type == 'over-represented')


png(file = './plots/Hans/binomial_cluster.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_Hans,
                      Z = hans_binomial_reorder1_pe$allocation,
                      regions.name = rownames(data_Hans[[1]]),
                      mouse.index = mouse.index)

dev.off()


png(file = './plots/Hans/binomial_cluster_5.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_Hans_5,
                      Z = hans_binomial_reorder5_pe$allocation,
                      regions.name = rownames(data_Hans[[1]]),
                      mouse.index = mouse.index)

dev.off()


#-------------------------- Variation of information between different estimates -------------------------------------



clusterings <- list(mcmc_unique_hans$Z,
                    hans_binomial_reorder1_pe$allocation,
                    hans_binomial_reorder5_pe$allocation,
                    clust20_r,
                    clust30_r,
                    clust40_r)


mx <- matrix(0, nrow = 6, ncol = 6)



colnames(mx) <- c('bayesian_motif',
                  'binomial_motif_1',
                  'binomial_motifs_5',
                  'k_means_20',
                  'k_means_30',
                  'k_means_40')

rownames(mx) <- c('bayesian_motif',
                  'binomial_motif_1',
                  'binomial_motifs_5',
                  'k_means_20',
                  'k_means_30',
                  'k_means_40')



for(i in 1:6){
  
  mx[i,] <- sapply(1:6,
                   function(v) variation_info(unlist(clusterings[[i]]),
                                              unlist(clusterings[[v]])))
}



mx <- mx/log(base = 2, x = sum(sapply(1:length(data_Hans),
                                      function(m) ncol(data_Hans[[m]]))))

mx <- round(mx, 2)


# Plot for the case with no threshold
ds <-  pivot_longer(as.data.frame(mx),
                    cols = 1:6,
                    names_to = "class.B", 
                    values_to = "vi"
)

ds$class.A <- rep(rownames(mx), each = 6)


png(file = './plots/Hans/comparison_vi.png',
    width = 700,
    height = 150)

ggplot(ds, aes(x = class.A, y = class.B, fill = vi)) +
  geom_tile() +
  labs(x = "Method A",
       y = "Method B") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "white", low="blue", midpoint = 0)

dev.off()







