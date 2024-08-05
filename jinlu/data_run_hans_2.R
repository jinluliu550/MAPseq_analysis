

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


# Compare different methods





#------------------------------ Binomial ---------------------------------

binomial_difference <- variation_info(unlist(hans_binomial_reorder5$allocation),
                                      unlist(hans_binomial_reorder1$allocation))/
  
  log(base = 2, x = sum(C))


binomial_difference_pe <- variation_info(unlist(hans_binomial_reorder5_pe$allocation),
                                         unlist(hans_binomial_reorder1_pe$allocation))/
  
  log(base = 2, x = sum(C))
