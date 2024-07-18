# Running on the barSeq data


# Load data
data1 <- read.csv('./data/barseq/BRAIN_C9.csv')
data2 <- read.csv('./data/barseq/BRAIN_C14.csv')
data3 <- read.csv('./data/barseq/BRAIN_C28.csv')


# Round to integers
data1 <- round(data1, 0)
data2 <- round(data2, 0)
data3 <- round(data3, 0)

data_barseq <- list(t(data1),
                    t(data2),
                    t(data3))

M <- length(data_barseq)
C <- sapply(1:M, function(m) ncol(data_barseq[[m]]))


mouse.index <- c(rep(1, C[1]),
                 rep(2, C[2]),
                 rep(3, C[3]))


#Initialize z
data_barseq_cbind <- do.call(cbind, data_barseq)

df <- t(data_barseq_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))


C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(data_barseq[[m]]))))

# k-means
k_mean_clust_60 <- kmeans(df, 60, iter.max = 100, nstart = 25)$cluster

clust60 <- lapply(1:M,
                  function(m) k_mean_clust_60[(C_cumsum[m]+1):C_cumsum[m+1]])



# MCMC run
mcmc_all_barseq <- mcmc_run_all(Y = data_barseq,
                                J = 80,
                                number_iter = 20000,
                                thinning = 5,
                                burn_in = 5000,
                                adaptive_prop = 0.0001,
                                print_Z = TRUE,
                                a_gamma = 30,
                                b_gamma = 1,
                                a_alpha = 1/5,
                                b_alpha = 1/2,
                                Z.init = clust60)


Zmat = matrix(unlist(mcmc_all_barseq$Z_output), length(mcmc_all_barseq$Z_output), sum(C),byrow = TRUE)

# Number of occupied components
k = apply(Zmat,1,function(x){length(unique(x))})
plot(k, type = 'l')


# Posterior similarity matrix
psm_barseq = similarity_matrix(mcmc_run_all_output = mcmc_all_barseq)


# Reordered posterior samples of z
barseq_z_reordered <- z_trace_updated(mcmc_run_all_output = mcmc_all_barseq)


# optimal clustering
barseq_Z <- opt.clustering.comb(z_trace = barseq_z_reordered,
                                post_similarity = psm_barseq,
                                max.k = max(k))


#-- Convert to a list
C_cumsum <- c(0, cumsum(C))
barseq_Z <- lapply(1:M,
                 function(m) barseq_Z[(C_cumsum[m]+1):C_cumsum[m+1]])


# Plot of posterior similarity matrix
psm_barseq_plot <- plotpsm(psm.ind = psm_barseq$psm.within,
                           psm.tot = psm_barseq$psm.combined)

png(file = './plots/barseq/heatmap_psm_1.png',
    width = 500,
    height = 400)

psm_barseq_plot$plot.ind

dev.off()


png(file = './plots/barseq/heatmap_psm_2.png',
    width = 500,
    height = 400)

psm_barseq_plot$plot.tot

dev.off()

# MCMC unique
mcmc_unique_barseq <- mcmc_run_post(mcmc_run_all_output = mcmc_all_barseq,
                                  Z = barseq_Z,
                                  thinning = 5,
                                  burn_in = 2000,
                                  number_iter = 12000,
                                  Y = data_barseq,
                                  a_gamma = 30,
                                  b_gamma = 1,
                                  regions.name = rownames(data_barseq[[1]]))

barseq_Z_reordered <- mcmc_unique_barseq$Z


# Number of neurons by cluster and mosue
png(file = './plots/barseq/number_of_neuron_by_m.png',
    width = 1200,
    height = 500)

opt.clustering.frequency(clustering = mcmc_unique_barseq$Z)

dev.off()


# Plot of estimated projection strength
png(file = './plots/barseq/estimated_pp.png',
    width = 1500,
    height = 800)

mcmc_unique_barseq$estimated.pp.plot

dev.off()


# q tilde
png(file = './plots/barseq/q_tilde.png',
    width = 1200,
    height = 300)

mcmc_unique_barseq$q_tilde_plot

dev.off()


# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_barseq,
                                   mcmc_run_post_output = mcmc_unique_barseq,
                                   thinning = 5,
                                   burn_in = 2000,
                                   number_iter = 12000)



png(file = './plots/barseq/w_jm.png',
    width = 1800,
    height = 900)

omega_JM_mcmc$omega_JM_plot

dev.off()

# Examining the difference of omega_jm between any 2 mice:
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc)


png(file = './plots/barseq/w_jm_difference.png',
    width = 1200,
    height = 300)

difference_omega_JM$probability_plot

dev.off()


png(file = './plots/barseq/w_jm_difference_eg.png',
    width = 800,
    height = 400)

difference_in_omega_jm_plot(difference_in_omega_jm_output = difference_omega_JM,
                            mcmc_run_omega_output = omega_JM_mcmc,
                            N = 6)

dev.off()


# Heatmap of projection stength of neurons in each cluster
png(file = './plots/barseq/heatmap_neuron.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_barseq,
                     Z = mcmc_unique_barseq$Z,
                     regions.name = rownames(data_barseq[[1]]),
                     mouse.index = mouse.index)

dev.off()


# Distribution of N_{i,m} for neurons within the same cluster
data_barseq_N <- lapply(1:length(data_barseq),
                      function(m) colSums(data_barseq[[m]]))

df <- data.frame(N = unlist(data_barseq_N),
                 motif = unlist(mcmc_unique_barseq$Z))

png(file = './plots/barseq/N_sum_by_cluster.png',
    width = 1000,
    height = 400)

ggplot(df)+
  geom_boxplot(mapping = aes(x = factor(motif, levels = 1:max(unlist(mcmc_unique_barseq$Z))),
                             y = N))+
  theme_bw()+
  xlab('Cluster')

dev.off()


# Posterior predictive check with multiple replicates
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_barseq,
                      Y = data_barseq,
                      N = 3,
                      regions.name = rownames(data_barseq[[1]]))


png(file = './plots/barseq/ppc_zero.png',
    width = 600,
    height = 250)

ppc_multiple$zero.plot

dev.off()


png(file = './plots/barseq/ppc_nonzero.png',
    width = 600,
    height = 250)

ppc_multiple$non.zero.plot

dev.off()


# Posterior predictive check with single replicated data
ppc_single <- ppc_single_f(mcmc_run_all_output = mcmc_all_barseq,
                           Y = data_barseq,
                           regions.name = rownames(data_barseq[[1]]))

for(m in 1:M){
  
  png(file = paste0('./plots/barseq/ppc_single_mouse_', m, '.png'),
      width = 1200,
      height = 400)
  
  print(ppc_single[[m]])
  
  dev.off()
}


#Initialize z
data_Hans_cbind <- do.call(cbind, data_barseq)

df <- t(data_Hans_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))

k_mean_clust_40 <- kmeans(df, 40, nstart = 25)$cluster

clust40 <- lapply(1:M,
                  function(m) k_mean_clust_40[(C_cumsum[m]+1):C_cumsum[m+1]])

clust40_r <- k_means_reorder(Y = data_barseq,
                             Z = clust40)

k_mean_clust_60 <- kmeans(df, 60, nstart = 25)$cluster

clust60 <- lapply(1:M,
                  function(m) k_mean_clust_60[(C_cumsum[m]+1):C_cumsum[m+1]])


clust60_r <- k_means_reorder(Y = data_barseq,
                             Z = clust60)


k_mean_clust_80 <- kmeans(df, 80, nstart = 25)$cluster

clust80 <- lapply(1:M,
                  function(m) k_mean_clust_80[(C_cumsum[m]+1):C_cumsum[m+1]])


clust80_r <- k_means_reorder(Y = data_barseq,
                             Z = clust80)


# Heatmap of projection strength of neurons in each cluster
png(file = './plots/barseq/k_means_40.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_barseq,
                     Z = clust40_r,
                     regions.name = rownames(data_barseq[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/barseq/k_means_60.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_barseq,
                     Z = clust60_r,
                     regions.name = rownames(data_barseq[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/barseq/k_means_80.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_barseq,
                     Z = clust80_r,
                     regions.name = rownames(data_barseq[[1]]),
                     mouse.index = mouse.index)

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


#--------------------------- Binomial clustering ----------------------------------

barseq_binomial <- binomial_model(data = data_barseq)

barseq_binomial_reorder <- binom_cluster_reorder(Y = data_barseq,
                                                 binomial_output = barseq_binomial)


png(file = './plots/barseq/binomial_cluster.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = data_barseq,
                     Z = barseq_binomial_reorder$allocation,
                     regions.name = rownames(data_barseq[[1]]),
                     mouse.index = mouse.index)

dev.off()

#-------------------------- Variation of information between different estimates -------------------------------------


# No1. Bayesian clustering
# No2. Binomial clustering

# No3. k-means with 40 clusters
# No4. k-means with 60 clusters
# No5. k-means with 80 clusters

clusterings <- list(mcmc_unique_barseq$Z,
                    barseq_binomial_reorder$allocation,
                    clust40_r,
                    clust60_r,
                    clust80_r)

mx <- matrix(0, nrow = 5, ncol = 5)

colnames(mx) <- c('bayesian_motif',
                  'binomial_motif',
                  'k_means_40',
                  'k_means_60',
                  'k_means_80')

rownames(mx) <- c('bayesian_motif',
                  'binomial_motif',
                  'k_means_40',
                  'k_means_60',
                  'k_means_80')


for(i in 1:5){
  
  mx[i,] <- sapply(1:5,
                   function(v) variation_info(unlist(clusterings[[i]]),
                                              unlist(clusterings[[v]])))
}

mx <- mx/log(base = 2, x = sum(sapply(1:length(data_barseq),
                                      function(m) ncol(data_barseq[[m]]))))

mx <- round(mx, 2)



ds <-  pivot_longer(as.data.frame(mx),
                    cols = 1:5,
                    names_to = "class.B", 
                    values_to = "vi"
)

# plot
ds$class.A <- rep(rownames(mx), each = 5)


png(file = './plots/barseq/comparison_vi.png',
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

evi_barseq <- evi.contribution(Zmat = Zmat,
                               Zhat = unlist(mcmc_unique_barseq$Z))


png(file = './plots/barseq/expected_vi.png',
    width = 3000,
    height = 1300)

projection_by_evi(Y = data_barseq,
                  evi = evi_barseq,
                  Z = mcmc_unique_barseq$Z, 
                  region_name = rownames(data_barseq[[1]]))

dev.off()


#------------------------ Add a 10 percent noise ----------------------------------

barseq_added_noise <- add_noise(mcmc_run_all_output = mcmc_all_barseq,
                              Y = data_barseq,
                              regions.name = rownames(data_barseq[[1]]))

barseq_added_noise_gelplot <- gel_plot(barseq_added_noise$noisy_data)
barseq_replicated_gelplot <- gel_plot(barseq_added_noise$replicated_data)


png(file = './plots/barseq/gelplot_added_noise.png',
    width = 1500,
    height = 500)

ggarrange(barseq_added_noise_gelplot[[1]],
          barseq_added_noise_gelplot[[2]],
          barseq_added_noise_gelplot[[3]],
          nrow = 1,
          widths = c(1,1,1.3))

dev.off()

png(file = './plots/barseq/gelplot_replicated.png',
    width = 1500,
    height = 500)

ggarrange(barseq_replicated_gelplot[[1]],
          barseq_replicated_gelplot[[2]],
          barseq_replicated_gelplot[[3]],
          nrow = 1,
          widths = c(1,1,1.3))

dev.off()

# Binomial clustering

#Initialize z

#---------------------------------- K-means noisy data ------------------------------------------------


df <- t(do.call(cbind, barseq_added_noise$noisy_data))
df = t(apply(df, 1, function(x){return(x/sum(x))}))


C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(barseq_added_noise$noisy_data[[m]]))))


# K = 40
k_mean_clust_40_added_noise <- kmeans(df, 40, nstart = 25)$cluster

clust40_added_noise <- lapply(1:M,
                              function(m) k_mean_clust_40_added_noise[(C_cumsum[m]+1):C_cumsum[m+1]])

clust40_r_added_noise <- k_means_reorder(Y = barseq_added_noise$noisy_data,
                                         Z = clust40_added_noise)


# K = 60
k_mean_clust_60_added_noise <- kmeans(df, 60, nstart = 25)$cluster

clust60_added_noise <- lapply(1:M,
                              function(m) k_mean_clust_60_added_noise[(C_cumsum[m]+1):C_cumsum[m+1]])


clust60_r_added_noise <- k_means_reorder(Y = barseq_added_noise$noisy_data,
                                         Z = clust60_added_noise)


# K = 80
k_mean_clust_80_added_noise <- kmeans(df, 80, nstart = 25)$cluster

clust80_added_noise <- lapply(1:M,
                              function(m) k_mean_clust_80_added_noise[(C_cumsum[m]+1):C_cumsum[m+1]])


clust80_r_added_noise <- k_means_reorder(Y = barseq_added_noise$noisy_data,
                                         Z = clust80_added_noise)


#--------------------------------- K-means replicated data -----------------------------------------


df <- t(do.call(cbind, barseq_added_noise$replicated_data))
df = t(apply(df, 1, function(x){return(x/sum(x))}))


# K = 40
k_mean_clust_40_replicated <- kmeans(df, 40, nstart = 25)$cluster

clust40_replicated <- lapply(1:M,
                              function(m) k_mean_clust_40_replicated[(C_cumsum[m]+1):C_cumsum[m+1]])

clust40_r_replicated <- k_means_reorder(Y = barseq_added_noise$replicated_data,
                                         Z = clust40_replicated)


# K = 60
k_mean_clust_60_replicated <- kmeans(df, 60, nstart = 25)$cluster

clust60_replicated <- lapply(1:M,
                             function(m) k_mean_clust_60_replicated[(C_cumsum[m]+1):C_cumsum[m+1]])

clust60_r_replicated <- k_means_reorder(Y = barseq_added_noise$replicated_data,
                                        Z = clust60_replicated)


# K = 80
k_mean_clust_80_replicated <- kmeans(df, 80, nstart = 25)$cluster

clust80_replicated <- lapply(1:M,
                             function(m) k_mean_clust_80_replicated[(C_cumsum[m]+1):C_cumsum[m+1]])

clust80_r_replicated <- k_means_reorder(Y = barseq_added_noise$replicated_data,
                                        Z = clust80_replicated)








# Heatmap of projection strength of neurons in each cluster
png(file = './plots/barseq/k_means_40_added_noise.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = barseq_added_noise$noisy_data,
                     Z = clust40_r_added_noise,
                     regions.name = rownames(barseq_added_noise$noisy_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/barseq/k_means_60_added_noise.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = barseq_added_noise$noisy_data,
                     Z = clust60_r_added_noise,
                     regions.name = rownames(barseq_added_noise$noisy_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/barseq/k_means_80_added_noise.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = barseq_added_noise$noisy_data,
                     Z = clust60_r_added_noise,
                     regions.name = rownames(barseq_added_noise$noisy_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/barseq/k_means_40_replicated.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = barseq_added_noise$replicated_data,
                     Z = clust40_r_replicated,
                     regions.name = rownames(barseq_added_noise$replicated_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/barseq/k_means_60_replicated.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = barseq_added_noise$replicated_data,
                     Z = clust60_r_replicated,
                     regions.name = rownames(barseq_added_noise$replicated_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/barseq/k_means_80_replicated.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = barseq_added_noise$replicated_data,
                     Z = clust60_r_replicated,
                     regions.name = rownames(barseq_added_noise$replicated_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

#--------------------------- Binomial clustering ----------------------------------

barseq_binomial_added_noise <- binomial_model(data = barseq_added_noise$noisy_data)

barseq_binomial_reorder_added_noise <- binom_cluster_reorder(Y = barseq_added_noise$noisy_data,
                                                 binomial_output = barseq_binomial_added_noise)


png(file = './plots/barseq/binomial_cluster_added_noise.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = barseq_added_noise$noisy_data,
                     Z = barseq_binomial_reorder_added_noise$allocation,
                     regions.name = rownames(barseq_added_noise$noisy_data[[1]]),
                     mouse.index = mouse.index)

dev.off()


barseq_binomial_replicated <- binomial_model(data = barseq_added_noise$replicated_data)

barseq_binomial_reorder_replicated <- binom_cluster_reorder(Y = barseq_added_noise$replicated_data,
                                                             binomial_output = barseq_binomial_replicated)

png(file = './plots/barseq/binomial_cluster_replicated.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = barseq_added_noise$replicated_data,
                     Z = barseq_binomial_reorder_replicated$allocation,
                     regions.name = rownames(barseq_added_noise$replicated_data[[1]]),
                     mouse.index = mouse.index)

dev.off()


#-------------------------------------------- MCMC of the noisy data ----------------------------------------------


mcmc_all_barseq_added_noise <- mcmc_run_all(Y = barseq_added_noise$noisy_data,
                                J = 80,
                                number_iter = 20000,
                                thinning = 5,
                                burn_in = 5000,
                                adaptive_prop = 0.0001,
                                print_Z = TRUE,
                                a_gamma = 30,
                                b_gamma = 1,
                                a_alpha = 1/5,
                                b_alpha = 1/2,
                                Z.init = clust60_r_added_noise)


#------------------------------------------ MCMC of the replicated data ------------------------------------------------


mcmc_all_barseq_replicated <- mcmc_run_all(Y = barseq_added_noise$replicated_data,
                                            J = 80,
                                            number_iter = 20000,
                                            thinning = 5,
                                            burn_in = 5000,
                                            adaptive_prop = 0.0001,
                                            print_Z = TRUE,
                                            a_gamma = 30,
                                            b_gamma = 1,
                                            a_alpha = 1/5,
                                            b_alpha = 1/2,
                                            Z.init = clust60_r_replicated)



Zmat_added_noise = matrix(unlist(mcmc_all_barseq_added_noise$Z_output), length(mcmc_all_barseq_added_noise$Z_output), sum(C),byrow = TRUE)

Zmat_replicated = matrix(unlist(mcmc_all_barseq_replicated$Z_output), length(mcmc_all_barseq_replicated$Z_output), sum(C),byrow = TRUE)


# Number of occupied components

k_added_noise = apply(Zmat_added_noise,1,function(x){length(unique(x))})
k_replicated = apply(Zmat_replicated,1,function(x){length(unique(x))})


# Posterior similarity matrix
psm_barseq_added_noise = similarity_matrix(mcmc_run_all_output = mcmc_all_barseq_added_noise)
psm_barseq_replicated = similarity_matrix(mcmc_run_all_output = mcmc_all_barseq_replicated)



# Reordered posterior samples of z
barseq_z_reordered_added_noise <- z_trace_updated(mcmc_run_all_output = mcmc_all_barseq_added_noise)
barseq_z_reordered_replicated <- z_trace_updated(mcmc_run_all_output = mcmc_all_barseq_replicated)

# optimal clustering
barseq_Z_added_noise <- opt.clustering.comb(z_trace = barseq_z_reordered_added_noise,
                                            post_similarity = psm_barseq_added_noise,
                                            max.k = max(k_added_noise))

barseq_Z_replicated <- opt.clustering.comb(z_trace = barseq_z_reordered_replicated,
                                          post_similarity = psm_barseq_replicated,
                                          max.k = max(k_replicated))


#-- Convert to a list
C_cumsum <- c(0, cumsum(C))
barseq_Z_added_noise <- lapply(1:M,
                             function(m) barseq_Z_added_noise[(C_cumsum[m]+1):C_cumsum[m+1]])


barseq_Z_replicated <- lapply(1:M,
                              function(m) barseq_Z_replicated[(C_cumsum[m]+1):C_cumsum[m+1]])

# MCMC unique
mcmc_unique_barseq_added_noise <- mcmc_run_post(mcmc_run_all_output = mcmc_all_barseq_added_noise,
                                              Z = barseq_Z_added_noise,
                                              thinning = 5,
                                              burn_in = 2000,
                                              number_iter = 12000,
                                              Y = barseq_added_noise$noisy_data,
                                              a_gamma = 30,
                                              b_gamma = 1,
                                              regions.name = rownames(barseq_added_noise$noisy_data[[1]]))

mcmc_unique_barseq_replicated <- mcmc_run_post(mcmc_run_all_output = mcmc_all_barseq_replicated,
                                                Z = barseq_Z_replicated,
                                                thinning = 5,
                                                burn_in = 2000,
                                                number_iter = 12000,
                                                Y = barseq_added_noise$replicated_data,
                                                a_gamma = 30,
                                                b_gamma = 1,
                                                regions.name = rownames(barseq_added_noise$replicated_data[[1]]))



# Heat-map of projection strength of neurons in each cluster
png(file = './plots/barseq/heatmap_neuron_added_noise.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = barseq_added_noise$noisy_data,
                     Z = mcmc_unique_barseq_added_noise$Z,
                     regions.name = rownames(barseq_added_noise$noisy_data[[1]]),
                     mouse.index = mouse.index)

dev.off()

png(file = './plots/barseq/heatmap_neuron_replicated.png',
    width = 600,
    height = 600)

pp.standard.ordering2(Y = barseq_added_noise$replicated_data,
                     Z = mcmc_unique_barseq_replicated$Z,
                     regions.name = rownames(barseq_added_noise$replicated_data[[1]]),
                     mouse.index = mouse.index)

dev.off()




#---------------------------------------- Compare variation of information ----------------------------------------------

added_noise_vi <- NULL
added_noise_vi$bayesian <- variation_info(unlist(mcmc_unique_barseq_replicated$Z),
                                          unlist(mcmc_unique_barseq_added_noise$Z))/log(base = 2, x = sum(sapply(1:length(data_barseq),
                                                                                                               function(m) ncol(data_barseq[[m]]))))

added_noise_vi$k_mean_40 <- variation_info(unlist(clust40_r_replicated),
                                           unlist(clust40_r_added_noise))/log(base = 2, x = sum(sapply(1:length(data_barseq),
                                                                                                       function(m) ncol(data_barseq[[m]]))))


added_noise_vi$k_mean_60 <- variation_info(unlist(clust60_r_replicated),
                                           unlist(clust60_r_added_noise))/log(base = 2, x = sum(sapply(1:length(data_barseq),
                                                                                                       function(m) ncol(data_barseq[[m]]))))


added_noise_vi$k_mean_80 <- variation_info(unlist(clust80_r_replicated),
                                           unlist(clust80_r_added_noise))/log(base = 2, x = sum(sapply(1:length(data_barseq),
                                                                                                       function(m) ncol(data_barseq[[m]]))))


added_noise_vi$binomial <- variation_info(unlist(barseq_binomial_reorder_replicated$allocation),
                                          unlist(barseq_binomial_reorder_added_noise$allocation))/log(base = 2, x = sum(sapply(1:length(data_barseq),
                                                                                                                             function(m) ncol(data_barseq[[m]]))))



#------------------------------------- Show how a Binomial cluster splits into multiple Bayesian cluster ----------------------------


# Combine binomial results with Bayesian results
df_compare <- lapply(1:M,
             function(m){
               
               data.frame(dataset = m,
                          neuron_index = 1:C[m],
                          binomial_allocation = barseq_binomial_reorder$allocation[[m]],
                          bayesian_allocation = mcmc_unique_barseq$Z[[m]])
             })

df_compare <- do.call(rbind, df_compare)


# Top 20 binomial motifs with the largest number of neurons
large_binomial_motifs <- df_compare %>%
  group_by(binomial_allocation) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  top_n(20) %>%
  select(binomial_allocation) %>%
  pull() %>%
  sort()


for(j in large_binomial_motifs){
  
  df_j <- df_compare %>%
    filter(binomial_allocation == j)
  
  
  df1 <- data.frame(neuron_index = rep(1:nrow(df_j), each = R),
                    bayesian_allocation = rep(df_j$bayesian_allocation, each = R),
                    projection_probability = unlist(lapply(1:nrow(df_j),
                                                           function(i) data_barseq[[df_j$dataset[i]]][,df_j$neuron_index[i]]/
                                                             sum(data_barseq[[df_j$dataset[i]]][,df_j$neuron_index[i]]))),
                    brain_region = rownames(data_barseq[[1]]))
  
  df1$bayesian_allocation <- factor(df1$bayesian_allocation,
                                    levels = unique(df1$bayesian_allocation))
  
  df1$brain_region <- factor(df1$brain_region,
                             levels = rownames(data_barseq[[1]]))
  
  png(filename = paste0('plots/barseq/binomial_motif_', j, '.png'),
      width = 800,
      height = 500)
  
  
  print(df1 %>%
          ggplot()+
          geom_line(mapping = aes(x = brain_region, y = projection_probability, color = bayesian_allocation, group = neuron_index))+
          theme_bw()+
          ylab('projection probability')+
          xlab('brain region')+
          ggtitle(paste('Cluster', j))
  )
  
  dev.off()
}

