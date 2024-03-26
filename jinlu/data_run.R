# Run the model on the MAPseq data

# pre-processing
LEC <- read.csv('./data/EC/lec_newSet_800_plusBrain.csv')
MEC <- read.csv('./data/EC/mec_newSet_800_plusBrain.csv')

R <- 9
M <- 6


# Load data
load("./data/EC/mcmc_all_sample.RData")
load("./data/EC/psm_EC.RData")
load("./data/EC/opt.clust0.EC.RData")
load("./data/EC/mcmc_unique_EC.RData")

digit.list <- list(one = 1,
                   two = 2,
                   three = 3,
                   four = 4,
                   five = 5,
                   six = 6)

# Change labels from text to numbers
LEC$brain <- sapply(1:nrow(LEC),
                    function(i) which(names(digit.list) == LEC$brain[i]))

MEC$brain <- sapply(1:nrow(MEC),
                    function(i) which(names(digit.list) == MEC$brain[i]))

# round data
data <- round(rbind(LEC, MEC), 0)


data_by_mouse <- lapply(1:M,
                        function(m) data %>%
                          filter(brain == m) %>%
                          select(-brain) %>%
                          as.matrix())

data_by_mouse <- lapply(1:M,
                        function(m) t(data_by_mouse[[m]]))




# Label of LEC and MEC of each mouse
LEC$EC_label <- 'LEC'
MEC$EC_label <- 'MEC'

data_by_mouse2 <- lapply(1:M,
                         function(m) rbind(LEC,MEC) %>%
                           filter(brain == m) %>%
                           select(EC_label) %>%
                           pull())

# A summary data frame
data3 <- data %>%
  mutate(neuron = c(paste('LEC neuron', 1:nrow(LEC)),
                    paste('MEC neuron', 1:nrow(MEC))))

data3_rowsums <- rowSums(data3[,1:R])

for(i in 1:R){
  
  data3[,i] <- data3[,i]/data3_rowsums
}


neuron_in_each_mouse <- lapply(1:M, 
                               function(m) data3 %>%
                                 filter(brain == m) %>%
                                 select(neuron) %>%
                                 pull())

C <- sapply(1:M,
            function(m) ncol(data_by_mouse[[m]]))





# MCMC run
mcmc_all_sample <- mcmc_run_all(Y = data_by_mouse,
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


psm_EC <- similarity_matrix(mcmc_run_all_output = mcmc_all_sample,
                         num.cores = 10,
                         run.on.pc = FALSE)

# Plot of posterior similarity matrix
plotpsm(psm.ind = psm_EC$psm.within,
        psm.tot = psm_EC$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)


plotpsm(psm.ind = psm_EC$psm.within,
        psm.tot = psm_EC$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')

# Optimal clustering
opt.clust0.EC <- dlso_cluster_estimate(mcmc_run_all_output = mcmc_all_sample)



# MCMC of unique parameters
mcmc_unique_EC <- mcmc_run_post(mcmc_run_all_output = mcmc_all_sample,
                             Z = opt.clust0.EC,
                             thinning = 5,
                             burn_in = 1500,
                             number_iter = 3000,
                             Y = data_by_mouse,
                             a_gamma = 500,
                             b_gamma = 10,
                             regions.name = rownames(data_by_mouse[[1]]))

# Number of neurons in each cluster - 148 clusters
png(file = './plots/EC/number_of_neuron.png',
    width = 2500,
    height = 700)

opt.clustering.frequency(clustering = mcmc_unique_EC$Z)

dev.off()

# Number of LEC and MEC neurons in each cluster
png(file = './plots/EC/number_of_neuron_by_EC.png',
    width = 2500,
    height = 700)

opt.clustering.frequency1(clustering = mcmc_unique_EC$Z,
                          EC_label = data_by_mouse2)

dev.off()

# Proportion of LEC and MEC in each cluster
png(file = './plots/EC/proportion_of_EC.png',
    width = 2500,
    height = 700)

opt.clustering.frequency2(clustering = mcmc_unique_EC$Z,
                          EC_label = data_by_mouse2)

dev.off()

# Plot of estimated projection strength
png(file = './plots/EC/estimated_pp.png',
    width = 3500,
    height = 2000)

mcmc_unique_EC$estimated.pp.plot

dev.off()

# Plot of q tilde
png(file = './plots/EC/q_tilde.png',
    width = 3000,
    height = 600)

mcmc_unique_EC$q_tilde_plot

dev.off()


# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_sample,
                                   mcmc_run_post_output = mcmc_unique_EC,
                                   thinning = 5,
                                   burn_in = 1500,
                                   number_iter = 3000)

png(file = './plots/EC/w_jm_EC.png',
    width = 4000,
    height = 2200)

omega_JM_mcmc$omega_JM_plot

dev.off()

# Projection strength of each neuron in each cluster, color-coded by the injection site

png(file = './plots/EC/projection_by_EC.png',
    width = 4000,
    height = 2200)

projection_by_EC(Y = data_by_mouse,
                 EC_label = data_by_mouse2,
                 Z = mcmc_unique_EC$Z,
                 region_name = rownames(data_by_mouse[[1]]))

dev.off()

# Examining the difference of omega_jm between any 2 mice:
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc)


png(file = './plots/EC/w_jm_difference.png',
    width = 2000,
    height = 600)

difference_omega_JM$probability_plot

dev.off()

png(file = './plots/EC/w_jm_difference_eg.png',
    width = 1500,
    height = 800)

difference_in_omega_jm_plot(difference_in_omega_jm_output = difference_omega_JM,
                            mcmc_run_omega_output = omega_JM_mcmc,
                            N = 20)

dev.off()


# Heat-map of projection strengths of each neuron

png(file = './plots/EC/heatmap_neuron.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = data_by_mouse,
                     Z = mcmc_unique_EC$Z,
                     regions.name = rownames(data_by_mouse[[1]]))

dev.off()

# Distribution of N_{i,m} for neurons within the same cluster
data_EC_N <- lapply(1:length(data_by_mouse),
                        function(m) colSums(data_by_mouse[[m]]))

df <- data.frame(N = unlist(data_EC_N),
                 motif = unlist(mcmc_unique_EC$Z))

png(file = './plots/EC/N_sum_by_cluster.png',
    width = 2500,
    height = 900)

ggplot(df)+
  geom_boxplot(mapping = aes(x = factor(motif, levels = 1:max(unlist(mcmc_unique_EC$Z))),
                             y = N))+
  theme_bw()+
  xlab('Cluster')

dev.off()

# Posterior predictive checks
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_sample,
                      Y = data_by_mouse,
                      N = 3,
                      regions.name = rownames(data_by_mouse[[1]]))


png(file = './plots/EC/ppc_zero.png',
    width = 1200,
    height = 700)

ppc_multiple$zero.plot

dev.off()

png(file = './plots/EC/ppc_nonzero.png',
    width = 1200,
    height = 700)

ppc_multiple$non.zero.plot

dev.off()

ppc_single <- ppc_single_f(mcmc_run_all_output = mcmc_all_sample,
                           Y = data_by_mouse,
                           regions.name = rownames(data_by_mouse[[1]]))

for(m in 1:M){
  
  png(file = paste0('./plots/EC/ppc_single_mouse_', m, '.png'),
      width = 1200,
      height = 700)
  
  print(ppc_single[[m]])
  
  dev.off()
}


#---------------------------------------------- Binomial results ----------------------------------------
binomial_result_EC <- binomial_model(data = data_by_mouse)

binomial_reorder_EC <- binom_cluster_reorder(Y = data_by_mouse,
                                             Z = binomial_result_EC$allocation)


png('plots/EC/binomial_heatmap.png',
    height = 900,
    width = 900)

pp.standard.ordering(Y = data_by_mouse,
                     Z = binomial_reorder_EC$Z,
                     regions.name = rownames(data_by_mouse[[1]]))

dev.off()


#---------------------------------------------------------------------------------------------------

data_EC_cbind <- do.call(cbind, data_by_mouse)

df <- scale(t(data_EC_cbind))

C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(data_by_mouse[[m]]))))

# K-means
k_mean_clust_30 <- kmeans(df, 30, nstart = 25)$cluster

clust30 <- lapply(1:6,
                  function(m) k_mean_clust_30[(C_cumsum[m]+1):C_cumsum[m+1]])

k_mean_clust_50 <- kmeans(df, 50, nstart = 25)$cluster

clust50 <- lapply(1:6,
                  function(m) k_mean_clust_50[(C_cumsum[m]+1):C_cumsum[m+1]])


k_mean_clust_70 <- kmeans(df, 70, nstart = 25)$cluster

clust70 <- lapply(1:6,
                  function(m) k_mean_clust_70[(C_cumsum[m]+1):C_cumsum[m+1]])


#---------------------------------------------------------------------------------------------------



# Compare results

# No1. Bayesian clustering
# No2. Binomial clustering

# No3. k-means with 30 clusters
# No4. k-means with 50 clusters
# No5. k-means with 70 clusters

clusterings <- list(mcmc_unique_EC$Z,
                    binomial_reorder_EC$Z,
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

mx <- mx/log(base = 2, x = sum(sapply(1:length(data_by_mouse),
                                      function(m) ncol(data_by_mouse[[m]]))))

mx <- round(mx, 2)


# Combine binomial results with Bayesian results
df <- lapply(1:M,
             function(m){
               
               data.frame(dataset = m,
                          neuron_index = 1:C[m],
                          binomial_allocation = binomial_reorder_EC$Z[[m]],
                          bayesian_allocation = mcmc_unique_EC$Z[[m]])
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
                                                           function(i) data_by_mouse[[df_j$dataset[i]]][,df_j$neuron_index[i]]/
                                                             sum(data_by_mouse[[df_j$dataset[i]]][,df_j$neuron_index[i]]))),
                    brain_region = rownames(data_by_mouse[[1]]))
  
  df1$binomial_allocation <- factor(df1$binomial_allocation,
                                    levels = unique(df1$binomial_allocation))
  
  
  png(filename = paste0('plots/EC/bayesian_motif_', j, '.png'),
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
                                                           function(i) data_by_mouse[[df_j$dataset[i]]][,df_j$neuron_index[i]]/
                                                             sum(data_by_mouse[[df_j$dataset[i]]][,df_j$neuron_index[i]]))),
                    brain_region = rownames(data_by_mouse[[1]]))
  
  df1$bayesian_allocation <- factor(df1$bayesian_allocation,
                                    levels = unique(df1$bayesian_allocation))
  
  
  png(filename = paste0('plots/EC/binomial_motif_', j, '.png'),
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





