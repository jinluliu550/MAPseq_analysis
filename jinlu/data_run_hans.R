
# Hans data - also called the VC data

M <- length(data_Hans)
C <- sapply(1:M, function(m) ncol(data_Hans[[m]]))

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


mouse.index <- c(rep(1, C[1]),
                 rep(2, C[2]),
                 rep(3, C[3]),
                 rep(4, C[4]))


#Initialize z
data_Hans_cbind <- do.call(cbind, data_Hans)

df <- t(data_Hans_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))


C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(data_Hans[[m]]))))

# k-means
k_mean_clust_20 <- kmeans(df, 20, iter.max = 100, nstart = 25)$cluster

clust20 <- lapply(1:M,
                  function(m) k_mean_clust_20[(C_cumsum[m]+1):C_cumsum[m+1]])


mcmc_all_hans <- mcmc_run_all(Y = data_Hans,
                              J = 40,
                              number_iter = 20000,
                              thinning = 5,
                              burn_in = 5000,
                              adaptive_prop = 0.0001,
                              print_Z = TRUE,
                              a_gamma = 30,
                              b_gamma = 1,
                              a_alpha = 1/5,
                              b_alpha = 1/2,
                              Z.init = clust20)


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


# Plot of posterior similarity matrix
psm_hans_plot <- plotpsm(psm.ind = psm_hans$psm.within,
                         psm.tot = psm_hans$psm.combined)


png(file = './plots/Hans/heatmap_psm_1.png',
    width = 500,
    height = 400)

psm_hans_plot$plot.ind

dev.off()


png(file = './plots/Hans/heatmap_psm_2.png',
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
                                  a_gamma = 30,
                                  b_gamma = 1,
                                  regions.name = rownames(data_Hans[[1]]))
                                 

hans_Z_reordered <- mcmc_unique_hans$Z



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

png(file = './plots/Hans/w_jm_difference_eg.png',
    width = 800,
    height = 400)

difference_in_omega_jm_plot(difference_in_omega_jm_output = difference_omega_JM,
                            mcmc_run_omega_output = omega_JM_mcmc,
                            N = 6)

dev.off()


# Heatmap of projection stength of neurons in each cluster
png(file = './plots/Hans/heatmap_neuron.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_Hans,
                      Z = mcmc_unique_hans$Z,
                      regions.name = rownames(data_Hans[[1]]),
                      mouse.index = mouse.index)

dev.off()


# Distribution of N_{i,m} for neurons within the same cluster
data_hans_N <- lapply(1:length(data_Hans),
                    function(m) colSums(data_Hans[[m]]))

df <- data.frame(N = unlist(data_hans_N),
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

#Initialize z
data_Hans_cbind <- do.call(cbind, data_Hans)

df <- t(data_Hans_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))

k_mean_clust_20 <- kmeans(df, 20, nstart = 25)$cluster

clust20 <- lapply(1:M,
                  function(m) k_mean_clust_20[(C_cumsum[m]+1):C_cumsum[m+1]])

clust20_r <- k_means_reorder(Y = data_Hans,
                             Z = clust20)

k_mean_clust_30 <- kmeans(df, 30, nstart = 25)$cluster

clust30 <- lapply(1:M,
                  function(m) k_mean_clust_30[(C_cumsum[m]+1):C_cumsum[m+1]])


clust30_r <- k_means_reorder(Y = data_Hans,
                             Z = clust30)


k_mean_clust_40 <- kmeans(df, 40, nstart = 25)$cluster

clust40 <- lapply(1:M,
                  function(m) k_mean_clust_40[(C_cumsum[m]+1):C_cumsum[m+1]])


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

hans_binomial <- binomial_model(data = data_Hans)

hans_binomial_reorder <- binom_cluster_reorder(Y = data_Hans,
                                               binomial_output = hans_binomial)


png(file = './plots/Hans/binomial_cluster.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_Hans,
                     Z = hans_binomial_reorder$allocation,
                     regions.name = rownames(data_Hans[[1]]),
                     mouse.index = mouse.index)

dev.off()

hans_binomial_reorder$cluster_summary %>%
  filter(significant == 'significant')

#-------------------------- Variation of information between different estimates -------------------------------------


# No1. Bayesian clustering
# No2. Binomial clustering

# No3. k-means with 20 clusters
# No4. k-means with 30 clusters
# No5. k-means with 40 clusters

clusterings <- list(mcmc_unique_hans$Z,
                    hans_binomial_reorder$allocation,
                    clust20_r,
                    clust30_r,
                    clust40_r)

mx <- matrix(0, nrow = 5, ncol = 5)

colnames(mx) <- c('bayesian_motif',
                  'binomial_motif',
                  'k_means_20',
                  'k_means_30',
                  'k_means_40')

rownames(mx) <- c('bayesian_motif',
                  'binomial_motif',
                  'k_means_20',
                  'k_means_30',
                  'k_means_40')


for(i in 1:5){
  
  mx[i,] <- sapply(1:5,
                   function(v) variation_info(unlist(clusterings[[i]]),
                                              unlist(clusterings[[v]])))
}

mx <- mx/log(base = 2, x = sum(sapply(1:length(data_Hans),
                                      function(m) ncol(data_Hans[[m]]))))

mx <- round(mx, 2)



ds <-  pivot_longer(as.data.frame(mx),
                    cols = 1:5,
                    names_to = "class.B", 
                    values_to = "vi"
)

# plot
ds$class.A <- rep(rownames(mx), each = 5)


png(file = './plots/Hans/comparison_vi.png',
    width = 500,
    height = 150)

ggplot(ds, aes(x = class.A, y = class.B, fill = vi)) +
  geom_tile() +
  labs(x = "Method A",
       y = "Method B") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "white", low="blue", midpoint = 0)

dev.off()


#--------------------------- Uncertainty in the neuron allocation ------------------------

evi_hans <- evi.contribution(Zmat = Zmat,
                             Zhat = unlist(mcmc_unique_hans$Z))

projection_by_evi(Y = data_Hans,
                  evi = evi_hans,
                  Z = mcmc_unique_hans$Z, 
                  region_name = rownames(data_Hans[[1]]))


#------------------------ Add a 10 percent noise ----------------------------------

hans_added_noise <- add_noise(mcmc_run_all_output = mcmc_all_hans,
                              Y = data_Hans,
                              regions.name = rownames(data_Hans[[1]]))

hans_added_noise_gelplot <- gel_plot(hans_added_noise$noisy_data)

hans_replicated_gelplot <- gel_plot(hans_added_noise$replicated_data)


png(file = './plots/Hans/gelplot_added_noise.png',
    width = 2000,
    height = 500)

ggarrange(hans_added_noise_gelplot[[1]],
          hans_added_noise_gelplot[[2]],
          hans_added_noise_gelplot[[3]],
          hans_added_noise_gelplot[[4]],
          nrow = 1,
          widths = c(1,1,1,1.3))

dev.off()

png(file = './plots/Hans/gelplot_replicated.png',
    width = 2000,
    height = 500)

ggarrange(hans_replicated_gelplot[[1]],
          hans_replicated_gelplot[[2]],
          hans_replicated_gelplot[[3]],
          hans_replicated_gelplot[[4]],
          nrow = 1,
          widths = c(1,1,1, 1.3))

dev.off()



#----------------------------- MCMC of the noisy data ----------------------------------------

M <- 4

df <- t(do.call(cbind, hans_added_noise$noisy_data))
df = t(apply(df, 1, function(x){return(x/sum(x))}))


C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(hans_added_noise$noisy_data[[m]]))))

# k-means
k_mean_clust_20_added_noise <- kmeans(df, 20, iter.max = 100, nstart = 25)$cluster

clust20_added_noise <- lapply(1:M,
                  function(m) k_mean_clust_20_added_noise[(C_cumsum[m]+1):C_cumsum[m+1]])


mcmc_all_hans_added_noise <- mcmc_run_all(Y = hans_added_noise$noisy_data,
                              J = 40,
                              number_iter = 20000,
                              thinning = 5,
                              burn_in = 5000,
                              adaptive_prop = 0.0001,
                              print_Z = TRUE,
                              a_gamma = 30,
                              b_gamma = 1,
                              a_alpha = 1/5,
                              b_alpha = 1/2,
                              Z.init = clust20_added_noise)


#---------------------------- MCMC of the replicated data ---------------------------------------

M <- 4

df <- t(do.call(cbind, hans_added_noise$replicated_data))
df = t(apply(df, 1, function(x){return(x/sum(x))}))


C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(hans_added_noise$replicated_data[[m]]))))

# k-means
k_mean_clust_20_replicated <- kmeans(df, 20, iter.max = 100, nstart = 25)$cluster

clust20_replicated <- lapply(1:M,
                              function(m) k_mean_clust_20_replicated[(C_cumsum[m]+1):C_cumsum[m+1]])


mcmc_all_hans_replicated <- mcmc_run_all(Y = hans_added_noise$replicated_data,
                                          J = 40,
                                          number_iter = 20000,
                                          thinning = 5,
                                          burn_in = 5000,
                                          adaptive_prop = 0.0001,
                                          print_Z = TRUE,
                                          a_gamma = 30,
                                          b_gamma = 1,
                                          a_alpha = 1/5,
                                          b_alpha = 1/2,
                                          Z.init = clust20_replicated)



#-------------------------------------------------------------------------------------------------



Zmat_added_noise = matrix(unlist(mcmc_all_hans_added_noise$Z_output), length(mcmc_all_hans_added_noise$Z_output), sum(C),byrow = TRUE)
Zmat_replicated = matrix(unlist(mcmc_all_hans_replicated$Z_output), length(mcmc_all_hans_replicated$Z_output), sum(C),byrow = TRUE)

# Number of occupied components
k_added_noise = apply(Zmat_added_noise,1,function(x){length(unique(x))})
k_replicated = apply(Zmat_replicated,1,function(x){length(unique(x))})

# Posterior similarity matrix
psm_hans_added_noise = similarity_matrix(mcmc_run_all_output = mcmc_all_hans_added_noise)
psm_hans_replicated = similarity_matrix(mcmc_run_all_output = mcmc_all_hans_replicated)

# Reordered posterior samples of z
hans_z_reordered_added_noise <- z_trace_updated(mcmc_run_all_output = mcmc_all_hans_added_noise)
hans_z_reordered_replicated <- z_trace_updated(mcmc_run_all_output = mcmc_all_hans_replicated)

# optimal clustering
hans_Z_added_noise <- opt.clustering.comb(z_trace = hans_z_reordered_added_noise,
                              post_similarity = psm_hans_added_noise,
                              max.k = max(k_added_noise))

hans_Z_replicated <- opt.clustering.comb(z_trace = hans_z_reordered_replicated,
                                          post_similarity = psm_hans_replicated,
                                          max.k = max(k_replicated))



#-- Convert to a list
C_cumsum <- c(0, cumsum(C))

hans_Z_added_noise <- lapply(1:M,
                 function(m) hans_Z_added_noise[(C_cumsum[m]+1):C_cumsum[m+1]])

hans_Z_replicated <- lapply(1:M,
                            function(m) hans_Z_replicated[(C_cumsum[m]+1):C_cumsum[m+1]])

# MCMC unique
mcmc_unique_hans_added_noise <- mcmc_run_post(mcmc_run_all_output = mcmc_all_hans_added_noise,
                                  Z = hans_Z_added_noise,
                                  thinning = 5,
                                  burn_in = 2000,
                                  number_iter = 12000,
                                  Y = hans_added_noise$noisy_data,
                                  a_gamma = 30,
                                  b_gamma = 1,
                                  regions.name = rownames(hans_added_noise$noisy_data[[1]]))

mcmc_unique_hans_replicated <- mcmc_run_post(mcmc_run_all_output = mcmc_all_hans_replicated,
                                              Z = hans_Z_replicated,
                                              thinning = 5,
                                              burn_in = 2000,
                                              number_iter = 12000,
                                              Y = hans_added_noise$replicated_data,
                                              a_gamma = 30,
                                              b_gamma = 1,
                                              regions.name = rownames(hans_added_noise$replicated_data[[1]]))


# Heatmap of projection stength of neurons in each cluster
png(file = './plots/Hans/heatmap_neuron_added_noise.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = hans_added_noise$noisy_data,
                     Z = mcmc_unique_hans_added_noise$Z,
                     regions.name = rownames(hans_added_noise$noisy_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/Hans/heatmap_neuron_replicated.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = hans_added_noise$replicated_data,
                     Z = mcmc_unique_hans_replicated$Z,
                     regions.name = rownames(hans_added_noise$replicated_data[[1]]),
                     mouse.index = mouse.index)

dev.off()


#--------------------------- K-means of the noisy data ------------------------------------------


df <- t(do.call(cbind, hans_added_noise$noisy_data))
df = t(apply(df, 1, function(x){return(x/sum(x))}))

# K = 20
k_mean_clust_20_added_noise <- kmeans(df, 20, nstart = 25)$cluster

clust20_added_noise <- lapply(1:M,
                              function(m) k_mean_clust_20_added_noise[(C_cumsum[m]+1):C_cumsum[m+1]])

clust20_r_added_noise <- k_means_reorder(Y = hans_added_noise$noisy_data,
                                         Z = clust20_added_noise)

# K = 30
k_mean_clust_30_added_noise <- kmeans(df, 30, nstart = 25)$cluster

clust30_added_noise <- lapply(1:M,
                              function(m) k_mean_clust_30_added_noise[(C_cumsum[m]+1):C_cumsum[m+1]])


clust30_r_added_noise <- k_means_reorder(Y = hans_added_noise$noisy_data,
                                         Z = clust30_added_noise)

# K = 40
k_mean_clust_40_added_noise <- kmeans(df, 40, nstart = 25)$cluster

clust40_added_noise <- lapply(1:M,
                              function(m) k_mean_clust_40_added_noise[(C_cumsum[m]+1):C_cumsum[m+1]])


clust40_r_added_noise <- k_means_reorder(Y = hans_added_noise$noisy_data,
                                         Z = clust40_added_noise)


#-------------------------- K-means of the replicated data ------------------------------------------



df <- t(do.call(cbind, hans_added_noise$replicated_data))
df = t(apply(df, 1, function(x){return(x/sum(x))}))


# K = 20
k_mean_clust_20_replicated <- kmeans(df, 20, nstart = 25)$cluster

clust20_replicated <- lapply(1:M,
                             function(m) k_mean_clust_20_replicated[(C_cumsum[m]+1):C_cumsum[m+1]])

clust20_r_replicated <- k_means_reorder(Y = hans_added_noise$replicated_data,
                                        Z = clust20_replicated)


# K = 30
k_mean_clust_30_replicated <- kmeans(df, 30, nstart = 25)$cluster

clust30_replicated <- lapply(1:M,
                             function(m) k_mean_clust_30_replicated[(C_cumsum[m]+1):C_cumsum[m+1]])

clust30_r_replicated <- k_means_reorder(Y = hans_added_noise$replicated_data,
                                        Z = clust30_replicated)

# K = 40
k_mean_clust_40_replicated <- kmeans(df, 40, nstart = 25)$cluster

clust40_replicated <- lapply(1:M,
                             function(m) k_mean_clust_40_replicated[(C_cumsum[m]+1):C_cumsum[m+1]])

clust40_r_replicated <- k_means_reorder(Y = hans_added_noise$replicated_data,
                                        Z = clust40_replicated)

#----------------------------------------------------------------------------------------------------


# Heatmap of projection strength of neurons in each cluster
png(file = './plots/Hans/k_means_20_added_noise.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = hans_added_noise$noisy_data,
                     Z = clust20_r_added_noise,
                     regions.name = rownames(hans_added_noise$noisy_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/Hans/k_means_30_added_noise.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = hans_added_noise$noisy_data,
                     Z = clust30_r_added_noise,
                     regions.name = rownames(hans_added_noise$noisy_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/Hans/k_means_40_added_noise.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = hans_added_noise$noisy_data,
                     Z = clust40_r_added_noise,
                     regions.name = rownames(hans_added_noise$noisy_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/Hans/k_means_20_replicated.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = hans_added_noise$replicated_data,
                     Z = clust20_r_replicated,
                     regions.name = rownames(hans_added_noise$replicated_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/Hans/k_means_30_replicated.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = hans_added_noise$replicated_data,
                     Z = clust30_r_replicated,
                     regions.name = rownames(hans_added_noise$replicated_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/Hans/k_means_40_replicated.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = hans_added_noise$replicated_data,
                     Z = clust40_r_replicated,
                     regions.name = rownames(hans_added_noise$replicated_data[[1]]),
                     mouse.index = mouse.index)

dev.off()


#--------------------------- Binomial clustering ----------------------------------

hans_binomial_added_noise <- binomial_model(data = hans_added_noise$noisy_data)

hans_binomial_reorder_added_noise <- binom_cluster_reorder(Y = hans_added_noise$noisy_data,
                                                           binomial_output = hans_binomial_added_noise)


png(file = './plots/Hans/binomial_cluster_added_noise.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = hans_added_noise$noisy_data,
                     Z = hans_binomial_reorder_added_noise$allocation,
                     regions.name = rownames(hans_added_noise$noisy_data[[1]]),
                     mouse.index = mouse.index)

dev.off()


hans_binomial_replicated <- binomial_model(data = hans_added_noise$replicated_data)

hans_binomial_reorder_replicated <- binom_cluster_reorder(Y = hans_added_noise$replicated_data,
                                                           binomial_output = hans_binomial_replicated)


png(file = './plots/Hans/binomial_cluster_replicated.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = hans_added_noise$replicated_data,
                     Z = hans_binomial_reorder_replicated$allocation,
                     regions.name = rownames(hans_added_noise$replicated_data[[1]]),
                     mouse.index = mouse.index)

dev.off()




# Compare variation of information
added_noise_vi <- NULL
added_noise_vi$bayesian <- variation_info(unlist(mcmc_unique_hans_replicated$Z),
                                          unlist(mcmc_unique_hans_added_noise$Z))/log(base = 2, x = sum(sapply(1:length(data_Hans),
                                                                                                               function(m) ncol(data_Hans[[m]]))))

added_noise_vi$k_mean_20 <- variation_info(unlist(clust20_r_replicated),
                                           unlist(clust20_r_added_noise))/log(base = 2, x = sum(sapply(1:length(data_Hans),
                                                                                                                function(m) ncol(data_Hans[[m]]))))


added_noise_vi$k_mean_30 <- variation_info(unlist(clust30_r_replicated),
                                           unlist(clust30_r_added_noise))/log(base = 2, x = sum(sapply(1:length(data_Hans),
                                                                                                       function(m) ncol(data_Hans[[m]]))))


added_noise_vi$k_mean_40 <- variation_info(unlist(clust40_r_replicated),
                                           unlist(clust40_r_added_noise))/log(base = 2, x = sum(sapply(1:length(data_Hans),
                                                                                                       function(m) ncol(data_Hans[[m]]))))


added_noise_vi$binomial <- variation_info(unlist(hans_binomial_reorder_replicated$allocation),
                                           unlist(hans_binomial_reorder_added_noise$allocation))/log(base = 2, x = sum(sapply(1:length(data_Hans),
                                                                                                       function(m) ncol(data_Hans[[m]]))))
