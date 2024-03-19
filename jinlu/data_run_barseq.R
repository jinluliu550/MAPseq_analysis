# Running on the barSeq data


# Load data
data1 <- read.csv('./data/Bar-seq-100010/BRAIN_C9.csv')
data2 <- read.csv('./data/Bar-seq-100010/BRAIN_C14.csv')
data3 <- read.csv('./data/Bar-seq-100010/BRAIN_C28.csv')

# Load data
load('data/Bar-seq-100010/psm_barseq.RData')
load('data/Bar-seq-100010/mcmc_all_barseq.RData')
load('data/Bar-seq-100010/opt.clust0.barseq.RData')

load('data/Bar-seq-100010/mcmc_unique_barseq.RData')


# Round to integers
data1 <- round(data1, 0)
data2 <- round(data2, 0)
data3 <- round(data3, 0)

data_barseq <- list(t(data1),
                    t(data2),
                    t(data3))

M <- length(data_barseq)
C <- sapply(1:M, function(m) ncol(data_barseq[[m]]))

# Save data
# save(data_barseq, file = './data/bar_seq.RData')

# MCMC run
mcmc_all_barseq <- mcmc_run_all(Y = data_barseq,
                                J = 150,
                                number_iter = 8000,
                                thinning = 5,
                                burn_in = 3000,
                                adaptive_prop = 0.1,
                                print_Z = TRUE,
                                
                                a_gamma = 1000,
                                b_gamma = 10,
                                a_alpha = 1/5,
                                b_alpha = 1/2,
                                num.cores = 10)



psm_barseq <- similarity_matrix(mcmc_run_all_output = mcmc_all_barseq,
                                num.cores = 10,
                                run.on.pc = FALSE)


# Plot of posterior similarity matrix
plotpsm(psm.ind = psm_barseq$psm.within,
        psm.tot = psm_barseq$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)


plotpsm(psm.ind = psm_barseq$psm.within,
        psm.tot = psm_barseq$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')


# Point estimate of z
opt.clust0.barseq <- dlso_cluster_estimate(mcmc_run_all_output = mcmc_all_barseq)
save(opt.clust0.barseq, file = 'data/Bar-seq-100010/opt.clust0.barseq.RData')


# MCMC of unique parameters
mcmc_unique_barseq <- mcmc_run_post(mcmc_run_all_output = mcmc_all_barseq,
                                    Z = opt.clust0.barseq,
                                    thinning = 5,
                                    burn_in = 1500,
                                    number_iter = 3000,
                                    Y = data_barseq,
                                    a_gamma = 1000,
                                    b_gamma = 10,
                                    regions.name = rownames(data_barseq[[1]])
                                    )

save(mcmc_unique_barseq, file = 'data/Bar-seq-100010/mcmc_unique_barseq.RData')

# Number of neurons in each cluster
png('./plots/Bar-Seq-100010/N_j_barseq.png',
    width = 2000,
    height = 600)

opt.clustering.frequency(clustering = mcmc_unique_barseq$Z)

dev.off()

# Plot of estimated projection strength
png('./plots/Bar-Seq-100010/pp_j_barseq.png',
    width = 3200,
    height = 1800)

mcmc_unique_barseq$estimated.pp.plot

dev.off()


# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc_barseq <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_barseq,
                                          mcmc_run_post_output = mcmc_unique_barseq,
                                          thinning = 5,
                                          burn_in = 1500,
                                          number_iter = 3000)

# plot                                   
png('./plots/Bar-Seq-100010/w_jm.png',
    width = 3500,
    height = 2000)

omega_JM_mcmc_barseq$omega_JM_plot

dev.off()

# Difference in omega_{j,m}
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc_barseq)

# plot
png('./plots/Bar-Seq-100010/w_jm_difference.png',
    width = 2000,
    height = 500)

difference_omega_JM$probability_plot

dev.off()

# Heat-map of projection strengths of each neuron
png('./plots/Bar-Seq-100010/projection_by_neuron.png',
    width = 937,
    height = 653)

pp.standard.ordering(Y = data_barseq,
                     Z = mcmc_unique_barseq$Z,
                     regions.name = rownames(data_barseq[[1]]))

dev.off()


# Posterior predictive checks
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_barseq,
                      Y = data_barseq,
                      N = 3,
                      regions.name = rownames(data_barseq[[1]]))

ppc_multiple$zero.plot
ppc_multiple$non.zero.plot


# Distribution of N_{i,m} for neurons within the same cluster
data_barseq_N <- lapply(1:length(data_barseq),
                        function(m) colSums(data_barseq[[m]]))

df <- data.frame(N = unlist(data_barseq_N),
                 motif = unlist(mcmc_unique_barseq$Z))

png('./plots/Bar-Seq-100010/N_distribution.png',
    width = 3000,
    height = 1200)

ggplot(df)+
  geom_boxplot(mapping = aes(x = factor(motif, levels = 1:max(unlist(mcmc_unique_barseq$Z))),
                             y = N))+
  theme_bw()+
  xlab('Cluster')

dev.off()

#-------------------------------------------------------------------------------------------------

# K-means clustering
data_barseq_cbind <- do.call(cbind, data_barseq)

df <- scale(t(data_barseq_cbind))

total.withinss <- data.frame(k = 1:50,
                             total.ss = sapply(1:50,
                                               function(i){
                                                 
                                                 
                                                 km.res <- kmeans(df, i, nstart = 25)
                                                 km.res$tot.withinss
                                               }))

ggplot(total.withinss,
       mapping = aes(x = k, y = total.ss))+
  geom_line(color = 'blue')+
  geom_point(color = 'blue')+
  xlab('Total number of clusters')+
  ylab('Total within sum of square')+
  theme_bw()


# Case with 30 clusters
k_mean_clust_30 <- kmeans(df, 30, nstart = 25)$cluster

clust30 <- list(k_mean_clust_30[1:605],
              k_mean_clust_30[606:(606+5082)],
              k_mean_clust_30[(606+5082+1):(605+5082+704)])

k_mean_reorder_cluster <- binom_cluster_reorder(Y = data_barseq,
                                                Z = clust30)



png('plots/Bar-Seq-100010/k_means_30.png',
    height = 600,
    width = 1000)

pp.standard.ordering(Y = data_barseq,
                     Z = k_mean_reorder_cluster$Z,
                     regions.name = rownames(data_barseq[[1]]))

dev.off()

# Case with 50 clusters
k_mean_clust_50 <- kmeans(df, 50, nstart = 25)$cluster

clust50 <- list(k_mean_clust_50[1:605],
              k_mean_clust_50[606:(606+5082)],
              k_mean_clust_50[(606+5082+1):(605+5082+704)])

k_mean_reorder_cluster <- binom_cluster_reorder(Y = data_barseq,
                                                Z = clust50)



png('plots/Bar-Seq-100010/k_means_50.png',
    height = 600,
    width = 1000)

pp.standard.ordering(Y = data_barseq,
                     Z = k_mean_reorder_cluster$Z,
                     regions.name = rownames(data_barseq[[1]]))

dev.off()

# Case with 70 clusters
k_mean_clust_70 <- kmeans(df, 70, nstart = 25)$cluster

clust70 <- list(k_mean_clust_70[1:605],
                k_mean_clust_70[606:(606+5082)],
                k_mean_clust_70[(606+5082+1):(605+5082+704)])

k_mean_reorder_cluster <- binom_cluster_reorder(Y = data_barseq,
                                                Z = clust70)

png('plots/Bar-Seq-100010/k_means_70.png',
    height = 600,
    width = 1000)

pp.standard.ordering(Y = data_barseq,
                     Z = k_mean_reorder_cluster$Z,
                     regions.name = rownames(data_barseq[[1]]))
dev.off()



#------------------------------------- Binomial clustering ---------------------------------------

# Reorder neurons
binomial_result <- binomial_model(data = data_barseq)

# Reorder clusters
binomial_reorder <- binom_cluster_reorder(Y = data_barseq,
                                          Z = binomial_result$allocation)

png('plots/Bar-Seq-100010/binomial_heatmap.png',
    height = 653,
    width = 937)

pp.standard.ordering(Y = data_barseq,
                     Z = binomial_reorder$Z,
                     regions.name = rownames(data_barseq[[1]]))

dev.off()

#---------------------------------------------------------------------------------------------------



# Compare results

# No1. Bayesian clustering
# No2. Binomial clustering

# No3. k-means with 30 clusters
# No4. k-means with 50 clusters
# No5. k-means with 70 clusters

clusterings <- list(mcmc_unique_barseq$Z,
                    binomial_reorder$Z,
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

mx <- mx/log(base = 2, x = sum(sapply(1:length(data_barseq),
                                      function(m) ncol(data_barseq[[m]]))))

mx <- round(mx, 2)


# Combine binomial results with Bayesian results
df <- lapply(1:M,
             function(m){
               
               data.frame(dataset = m,
                          neuron_index = 1:C[m],
                          binomial_allocation = binomial_reorder$Z[[m]],
                          bayesian_allocation = mcmc_unique_barseq$Z[[m]])
             })

df <- do.call(rbind, df)

# Find the bayesian motifs with more than 100 allocations
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
                                                           function(i) data_barseq[[df_j$dataset[i]]][,df_j$neuron_index[i]]/
                                                             sum(data_barseq[[df_j$dataset[i]]][,df_j$neuron_index[i]]))),
                    brain_region = rownames(data_barseq[[1]]))
  
  df1$binomial_allocation <- factor(df1$binomial_allocation,
                                    levels = unique(df1$binomial_allocation))
  
  
  png(filename = paste0('plots/Bar-Seq-100010/bayesian_motif_', j, '.png'),
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
                                                           function(i) data_barseq[[df_j$dataset[i]]][,df_j$neuron_index[i]]/
                                                             sum(data_barseq[[df_j$dataset[i]]][,df_j$neuron_index[i]]))),
                    brain_region = rownames(data_barseq[[1]]))
  
  df1$bayesian_allocation <- factor(df1$bayesian_allocation,
                                    levels = unique(df1$bayesian_allocation))
  
  
  png(filename = paste0('plots/Bar-Seq-100010/binomial_motif_', j, '.png'),
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
