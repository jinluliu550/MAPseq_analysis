EC8_new <- read.csv('./data/EC8_new/all_brains_setE_JL.csv')


load('./data/EC8_new/mcmc_all_EC8.RData')
load('./data/EC8_new/psm_EC8.RData')
load('./data/EC8_new/EC8_Z_minvi.RData')
load('./data/EC8_new/EC8_Z_dlso.RData')
load('./data/EC8_new/EC8_Z.RData')
load('./data/EC8_new/mcmc_unique_EC8.RData')

#---------------------------------------- Binomial results --------------------------------------------

cluster_label <- read.csv('./edward/setE_results/cluster_label.csv',
                          header = FALSE)$V1

# Allocations
load('./data/EC8_new/setE_binom_allocations.RData')


cluster_summary <- read.csv('./edward/setE_results/cluster_summary.csv')

binomial_output <- list('allocation' = allocation,
                        'cluster_label' = cluster_label,
                        'cluster_summary' = cluster_summary)

binomial_output_reorder <- binom_cluster_reorder(Y = EC8_new,
                                                 binomial_output = binomial_output)

#--------------------------------------------------------------------------------------------------------

# Convert to list
digit.list <- list(one = 1,
                   two = 2,
                   three = 3,
                   four = 4,
                   five = 5,
                   six = 6)

EC8_new_list <- lapply(1:6, 
                       function(m) EC8_new %>%
                         filter(brain == names(digit.list)[m]))

EC8_new <- lapply(1:6,
                  function(m) t(round(EC8_new_list[[m]][,1:8], 0)))

EC8_EC_label <- lapply(1:6,
                       function(m) EC8_new_list[[m]]$EC)


C <- sapply(1:6, function(m) ncol(EC8_new[[m]]))
R <- 8

mcmc_all_EC8 <- mcmc_run_all(Y = EC8_new,
                             J = 150,
                             number_iter = 8000,
                             thinning = 5,
                             burn_in = 3000,
                             adaptive_prop = 0.1,
                             print_Z = TRUE,
                             
                             
                             a_gamma = 500,
                             b_gamma = 10,
                             a_alpha = 1/5,
                             b_alpha = 1/2,
                             num.cores = 10)


psm_EC8 <- similarity_matrix(mcmc_run_all_output = mcmc_all_EC8,
                             num.cores = 10,
                             run.on.pc = FALSE)



# minvi
EC8_Z_minvi <- opt.clustering(mcmc_run_all_output = mcmc_all_EC8,
                              post_similarity = psm_EC8)

# dlso
EC8_Z_dlso <- dlso_cluster_estimate(mcmc_run_all_output = mcmc_all_EC8)

# optimal clustering - 121 clusters
EC8_Z <- opt.clustering.comb(mcmc_run_all_output = mcmc_all_EC8,
                             post_similarity = psm_EC8)


# Plot of posterior similarity matrix
plotpsm(psm.ind = psm_EC8$psm.within,
        psm.tot = psm_EC8$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)


plotpsm(psm.ind = psm_EC8$psm.within,
        psm.tot = psm_EC8$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')

# MCMC unique
mcmc_unique_EC8 <- mcmc_run_post(mcmc_run_all_output = mcmc_all_EC8,
                                 Z = EC8_Z,
                                 thinning = 5,
                                 burn_in = 1500,
                                 number_iter = 3000,
                                 Y = EC8_new,
                                 a_gamma = 500,
                                 b_gamma = 10,
                                 regions.name = rownames(EC8_new[[1]]))

# Number of neurons in each cluster
png(file = './plots/EC8_new/number_of_neuron.png',
    width = 2500,
    height = 700)

opt.clustering.frequency(clustering = mcmc_unique_EC8$Z)

dev.off()

# Number of LEC and MEC neurons in each cluster
png(file = './plots/EC8_new/number_of_neuron_by_EC.png',
    width = 2500,
    height = 700)

opt.clustering.frequency1(clustering = mcmc_unique_EC8$Z,
                          EC_label = EC8_EC_label)

dev.off()


# Proportion of LEC and MEC in each cluster
png(file = './plots/EC8_new/proportion_of_EC.png',
    width = 2500,
    height = 700)

opt.clustering.frequency2(clustering = mcmc_unique_EC8$Z,
                          EC_label = EC8_EC_label)

dev.off()

# Plot of estimated projection strength
png(file = './plots/EC8_new/estimated_pp.png',
    width = 3500,
    height = 2000)

mcmc_unique_EC8$estimated.pp.plot

dev.off()

# q tilde
png(file = './plots/EC8_new/q_tilde.png',
    width = 3000,
    height = 600)

mcmc_unique_EC8$q_tilde_plot

dev.off()

# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_EC8,
                                   mcmc_run_post_output = mcmc_unique_EC8,
                                   thinning = 5,
                                   burn_in = 1500,
                                   number_iter = 3000)

png(file = './plots/EC8_new/w_jm_EC.png',
    width = 4000,
    height = 2200)

omega_JM_mcmc$omega_JM_plot

dev.off()


# Projection strength of each neuron in each cluster, color-coded by the injection site
png(file = './plots/EC8_new/projection_by_EC.png',
    width = 4000,
    height = 2200)

projection_by_EC(Y = EC8_new,
                 EC_label = EC8_EC_label,
                 Z = mcmc_unique_EC8$Z,
                 region_name = rownames(EC8_new[[1]]))

dev.off()


# Examining the difference of omega_jm between any 2 mice:
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc)


png(file = './plots/EC8_new/w_jm_difference.png',
    width = 2000,
    height = 600)

difference_omega_JM$probability_plot

dev.off()

png(file = './plots/EC8_new/w_jm_difference_eg.png',
    width = 1500,
    height = 800)

difference_in_omega_jm_plot(difference_in_omega_jm_output = difference_omega_JM,
                            mcmc_run_omega_output = omega_JM_mcmc,
                            N = 20)

dev.off()

for(m in 1:M){
  
  png(file = paste0('./plots/EC8_new/w_jm_difference_', m, '.png'),
      width = 1500,
      height = 800)
  
  difference_omega_JM$plot.for.each.data[[m]]
  
  dev.off()
}

# Heatmap of projection stength of neurons in each cluster
png(file = './plots/EC8_new/heatmap_neuron.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = EC8_new,
                     Z = mcmc_unique_EC8$Z,
                     regions.name = rownames(EC8_new[[1]]))

dev.off()


# Distribution of N_{i,m} for neurons within the same cluster
data_EC_N <- lapply(1:length(EC8_new),
                    function(m) colSums(EC8_new[[m]]))

df <- data.frame(N = unlist(data_EC_N),
                 motif = unlist(mcmc_unique_EC8$Z))

png(file = './plots/EC8_new/N_sum_by_cluster.png',
    width = 2500,
    height = 900)

ggplot(df)+
  geom_boxplot(mapping = aes(x = factor(motif, levels = 1:max(unlist(mcmc_unique_EC8$Z))),
                             y = N))+
  theme_bw()+
  xlab('Cluster')

dev.off()

# Posterior predictive check with multiple replicates
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_EC8,
                      Y = EC8_new,
                      N = 3,
                      regions.name = rownames(EC8_new[[1]]))


png(file = './plots/EC8_new/ppc_zero.png',
    width = 1200,
    height = 700)

ppc_multiple$zero.plot

dev.off()

png(file = './plots/EC8_new/ppc_nonzero.png',
    width = 1200,
    height = 700)

ppc_multiple$non.zero.plot

dev.off()

# Posterior predictive check with single replicated data
ppc_single <- ppc_single_f(mcmc_run_all_output = mcmc_all_EC8,
                           Y = EC8_new,
                           regions.name = rownames(EC8_new[[1]]))

for(m in 1:length(EC8_new)){
  
  png(file = paste0('./plots/EC8_new/ppc_single_mouse_', m, '.png'),
      width = 1200,
      height = 700)
  
  print(ppc_single[[m]])
  
  dev.off()
}


#---------------------------------------------------------------------------------------------------

data_EC_cbind <- do.call(cbind, EC8_new)

df <- scale(t(data_EC_cbind))

C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(EC8_new[[m]]))))

# K-means
k_mean_clust_30 <- kmeans(df, 30, nstart = 25)$cluster

clust30 <- lapply(1:6,
                  function(m) k_mean_clust_30[(C_cumsum[m]+1):C_cumsum[m+1]])

clust30_r <- k_means_reorder(Y = EC8_new,
                             Z = clust30)

k_mean_clust_50 <- kmeans(df, 50, nstart = 25)$cluster

clust50 <- lapply(1:6,
                  function(m) k_mean_clust_50[(C_cumsum[m]+1):C_cumsum[m+1]])


clust50_r <- k_means_reorder(Y = EC8_new,
                             Z = clust50)


k_mean_clust_70 <- kmeans(df, 70, nstart = 25)$cluster

clust70 <- lapply(1:6,
                  function(m) k_mean_clust_70[(C_cumsum[m]+1):C_cumsum[m+1]])


clust70_r <- k_means_reorder(Y = EC8_new,
                             Z = clust70)

# Heatmap of projection strength of neurons in each cluster
png(file = './plots/EC8_new/k_means_30.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = EC8_new,
                     Z = clust30_r,
                     regions.name = rownames(EC8_new[[1]]))

dev.off()

png(file = './plots/EC8_new/k_means_50.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = EC8_new,
                     Z = clust50_r,
                     regions.name = rownames(EC8_new[[1]]))

dev.off()

png(file = './plots/EC8_new/k_means_70.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = EC8_new,
                     Z = clust70_r,
                     regions.name = rownames(EC8_new[[1]]))

dev.off()

#-------------------------------- Comparison ----------------------------------------------------

png(file = './plots/EC8_new/heatmap_binomial.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = EC8_new,
                     Z = binomial_output_reorder$allocation,
                     regions.name = rownames(EC8_new[[1]]))

dev.off()

#---------------------------------------------------------------------------------------------------



# Compare results

# No1. Bayesian clustering
# No2. Binomial clustering

# No3. k-means with 30 clusters
# No4. k-means with 50 clusters
# No5. k-means with 70 clusters

clusterings <- list(mcmc_unique_EC8$Z,
                    binomial_output_reorder$allocation,
                    clust30,
                    clust50,
                    clust70)

df1 <- sapply(1:5,
              function(i) variation_info(unlist(clusterings[[1]]),
                                         unlist(clusterings[[i]])))

df2 <- sapply(2:5,
              function(i) variation_info(unlist(clusterings[[2]]),
                                         unlist(clusterings[[i]])))

df3 <- sapply(3:5,
              function(i) variation_info(unlist(clusterings[[3]]),
                                         unlist(clusterings[[i]])))

df4 <- sapply(4:5,
              function(i) variation_info(unlist(clusterings[[4]]),
                                         unlist(clusterings[[i]])))


mx <- matrix(0, nrow = 5, ncol = 5)
mx[1,] <- df1
mx[2,] <- c(NA, df2)
mx[3,] <- c(NA, NA, df3)
mx[4,] <- c(NA, NA, NA, df4)
mx[5,] <- c(rep(NA,4), 0)

mx <- round(mx, 2)

colnames(mx) <- c('bayesian_motif',
                  'binomial_motif',
                  'k_means_30',
                  'k_means_50',
                  'k_means_70')

rownames(mx) <- c('bayesian_motif',
                  'binomial_motif',
                  'k_means_30',
                  'k_means_50',
                  'k_means_70')


mx <- mx/log(base = 2, x = sum(sapply(1:length(EC8_new),
                                      function(m) ncol(EC8_new[[m]]))))

mx <- round(mx, 2)


# Combine binomial results with Bayesian results
df <- lapply(1:M,
             function(m){
               
               data.frame(dataset = m,
                          neuron_index = 1:C[m],
                          binomial_allocation = binomial_output_reorder$allocation[[m]],
                          bayesian_allocation = mcmc_unique_EC8$Z[[m]])
             })

df <- do.call(rbind, df)

# Find the Bayesian motifs with more than 100 allocations
large_bayesian <- df %>%
  group_by(bayesian_allocation) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  filter(count > 100) %>%
  select(bayesian_allocation) %>%
  pull()

for(j in large_bayesian){
  
  df_j <- df %>%
    filter(bayesian_allocation == j)
  
  df1 <- data.frame(neuron_index = rep(1:nrow(df_j), each = R),
                    binomial_allocation = rep(df_j$binomial_allocation, each = R),
                    projection_probability = unlist(lapply(1:nrow(df_j),
                                                           function(i) EC8_new[[df_j$dataset[i]]][,df_j$neuron_index[i]]/
                                                             sum(EC8_new[[df_j$dataset[i]]][,df_j$neuron_index[i]]))),
                    brain_region = rownames(EC8_new[[1]]))
  
  df1$binomial_allocation <- factor(df1$binomial_allocation,
                                    levels = unique(df1$binomial_allocation))
  
  df1$brain_region <- factor(df1$brain_region,
                             levels = rownames(EC8_new[[1]]))
  
  png(filename = paste0('plots/EC8_new/bayesian_motif_', j, '.png'),
      width = 800,
      height = 500)
  
  
  print(df1 %>%
          ggplot()+
          geom_line(mapping = aes(x = brain_region, y = projection_probability, color = binomial_allocation, group = neuron_index))+
          theme_bw()+
          ylab('projection probability')+
          xlab('brain region')+
          ggtitle(paste('Cluster', j))
  )
  
  dev.off()
}


# Find the binomial motifs with more than 100 allocations
large_binomial <- df %>%
  group_by(binomial_allocation) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  filter(count > 100) %>%
  select(binomial_allocation) %>%
  pull()


for(j in large_binomial){
  
  df_j <- df %>%
    filter(binomial_allocation == j)
  
  df1 <- data.frame(neuron_index = rep(1:nrow(df_j), each = R),
                    bayesian_allocation = rep(df_j$bayesian_allocation, each = R),
                    projection_probability = unlist(lapply(1:nrow(df_j),
                                                           function(i) EC8_new[[df_j$dataset[i]]][,df_j$neuron_index[i]]/
                                                             sum(EC8_new[[df_j$dataset[i]]][,df_j$neuron_index[i]]))),
                    brain_region = rownames(EC8_new[[1]]))
  
  df1$bayesian_allocation <- factor(df1$bayesian_allocation,
                                    levels = unique(df1$bayesian_allocation))
  
  df1$brain_region <- factor(df1$brain_region,
                             levels = rownames(EC8_new[[1]]))
  
  png(filename = paste0('plots/EC8_new/binomial_motif_', j, '.png'),
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





