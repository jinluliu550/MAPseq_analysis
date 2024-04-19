EC8_new <- read.csv('./data/EC8_new/all_brains_setE_JL.csv')



# Acceptance probabilities
plot(mcmc_all_EC8$acceptance_prob$omega, type = 'l')
plot(mcmc_all_EC8$acceptance_prob$alpha, type = 'l')
plot(mcmc_all_EC8$acceptance_prob$alpha_zero, type = 'l')
plot(mcmc_all_EC8$acceptance_prob$q_star, type = 'l')
plot(mcmc_all_EC8$acceptance_prob$gamma_star, type = 'l')
plot(mcmc_all_EC8$acceptance_prob$alpha_h, type = 'l')

# Trace plots
plot(mcmc_all_EC8$alpha_output, type = 'l')
plot(mcmc_all_EC8$alpha_zero_output, type = 'l')


#---------------------------------------- Binomial results --------------------------------------------

cluster_label <- read.csv('./edward/setE_results/cluster_label.csv',
                          header = FALSE)$V1

# Allocations
load('./data/EC8_new/setE_binom_allocations.RData')

load('./edward/setE_results/setE_motif_groups.RData')



cluster_summary <- read.csv('./edward/setE_results/cluster_summary.csv')


binomial_output <- list('allocation' = allocation,
                        'cluster_label' = cluster_label,
                        'cluster_summary' = cluster_summary)

binomial_output_reorder <- binom_cluster_reorder(Y = EC8_new,
                                                 binomial_output = binomial_output)

binomial_output_reorder$cluster_summary <- NULL

save(binomial_output_reorder, file = './data/EC8_new/binomial_output_reorder.RData')


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

#-------------------------------- Empirical Analysis ---------------------------------------

# Number of neurons in each mouse
png(file = './plots/EC8_new/number_of_neurons_in_each_m.png')

data.frame(mouse = paste('mouse', 1:M),
           C = C) %>%
  ggplot()+
  geom_bar(mapping = aes(x = mouse,
                         y = C),
           stat = 'identity')+
  ylab('Number of neurons in each mouse')+
  theme_bw()

dev.off()

# Proportion of LEC and MEC neurons in each mouse
proportion_of_EC <- lapply(1:6,
                           function(m) data.frame(mouse = paste('mouse', m),
                                                  EC = c('LEC','MEC'),
                                                  count = c(length(which(EC8_EC_label[[m]] == 'LEC')),
                                                            length(which(EC8_EC_label[[m]] == 'MEC')))
                                                  )
                           )

proportion_of_EC <- do.call(rbind, proportion_of_EC)

png(file = './plots/EC8_new/prop_of_EC_in_each_m.png')

ggplot(proportion_of_EC, aes(x = mouse, y = count, fill = EC))+
  geom_bar(position = 'fill', stat = 'identity')+
  ylab('proportion of LEC and MEC')+
  theme_bw()

dev.off()

#-------------------------------------------------------------------------------------------

#Initialize z
data_EC_cbind <- do.call(cbind, EC8_new)

df <- t(data_EC_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))

C_cumsum <- c(0, cumsum(sapply(1:6, function(m) ncol(EC8_new[[m]]))))

# k-means
k_mean_clust_70 <- kmeans(df, 70, iter.max = 100, nstart = 25)$cluster

clust70 <- lapply(1:6,
                  function(m) k_mean_clust_70[(C_cumsum[m]+1):C_cumsum[m+1]])



mcmc_all_EC8 <- mcmc_run_all(Y = EC8_new,
                             J = 100,
                             number_iter = 15000,
                             thinning = 5,
                             burn_in = 5000,
                             adaptive_prop = 0.0001,
                             print_Z = TRUE,
                             a_gamma = 30,
                             b_gamma = 1,
                             a_alpha = 1/5,
                             b_alpha = 1/2,
                             Z.init = clust70)


Zmat = matrix(unlist(mcmc_all_EC8$Z_output), length(mcmc_all_EC8$Z_output), sum(C),byrow = TRUE)

# Number of occupied components
k = apply(Zmat,1,function(x){length(unique(x))})

# Posterior similarity matrix
psm_EC8 = similarity_matrix(mcmc_run_all_output = mcmc_all_EC8)

# Reordered posterior samples of z
EC8_z_reordered <- z_trace_updated(mcmc_run_all_output = mcmc_all_EC8)

# optimal clustering
EC8_Z <- opt.clustering.comb(z_trace = EC8_z_reordered,
                             post_similarity = psm_EC8,
                             max.k = max(k))


#-- Convert to a list
C_cumsum <- c(0, cumsum(C))
EC8_Z <- lapply(1:6,
                           function(m) EC8_Z[(C_cumsum[m]+1):C_cumsum[m+1]])

# Plot of posterior similarity matrix
png(file = './plots/EC8_new/heatmap_psm_1_new.png',
    width = 664,
    height = 664)

plotpsm(psm.ind = psm_EC8$psm.within,
        psm.tot = psm_EC8$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)

dev.off()

png(file = './plots/EC8_new/heatmap_psm_2_new.png',
    width = 664,
    height = 664)

plotpsm(psm.ind = psm_EC8$psm.within,
        psm.tot = psm_EC8$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')

dev.off()

# MCMC unique
mcmc_unique_EC8 <- mcmc_run_post(mcmc_run_all_output = mcmc_all_EC8,
                                 Z = EC8_Z,
                                 thinning = 5,
                                 burn_in = 2000,
                                 number_iter = 12000,
                                 Y = EC8_new,
                                 a_gamma = 30,
                                 b_gamma = 1,
                                 regions.name = rownames(EC8_new[[1]]))

EC8_Z_reordered <- mcmc_unique_EC8$Z

# Number of neurons in each cluster
png(file = './plots/EC8_new/number_of_neuron_new.png',
    width = 2500,
    height = 700)

opt.clustering.frequency(clustering = mcmc_unique_EC8$Z)

dev.off()

# Number of LEC and MEC neurons in each cluster
png(file = './plots/EC8_new/number_of_neuron_by_EC_new.png',
    width = 2500,
    height = 700)

opt.clustering.frequency1(clustering = mcmc_unique_EC8$Z,
                          EC_label = EC8_EC_label)

dev.off()


# Proportion of LEC and MEC in each cluster
png(file = './plots/EC8_new/proportion_of_EC_new.png',
    width = 2500,
    height = 700)

opt.clustering.frequency2(clustering = mcmc_unique_EC8$Z,
                          EC_label = EC8_EC_label)

dev.off()

# Plot of estimated projection strength
png(file = './plots/EC8_new/estimated_pp_new.png',
    width = 3500,
    height = 2000)

mcmc_unique_EC8$estimated.pp.plot

dev.off()

# q tilde
png(file = './plots/EC8_new/q_tilde_new.png',
    width = 3000,
    height = 600)

mcmc_unique_EC8$q_tilde_plot

dev.off()

# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_EC8,
                                   mcmc_run_post_output = mcmc_unique_EC8,
                                   thinning = 5,
                                   burn_in = 2000,
                                   number_iter = 12000)

png(file = './plots/EC8_new/w_jm_EC_new.png',
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


png(file = './plots/EC8_new/w_jm_difference_new.png',
    width = 2000,
    height = 600)

difference_omega_JM$probability_plot

dev.off()

png(file = './plots/EC8_new/w_jm_difference_eg_new.png',
    width = 1500,
    height = 800)

difference_in_omega_jm_plot(difference_in_omega_jm_output = difference_omega_JM,
                            mcmc_run_omega_output = omega_JM_mcmc,
                            N = 20)

dev.off()

for(m in 1:M){
  
  png(file = paste0('./plots/EC8_new/w_jm_difference_new_', m, '.png'),
      width = 1500,
      height = 800)
  
  print(difference_omega_JM$plot.for.each.data[[m]])
  
  dev.off()
}

# Heatmap of projection stength of neurons in each cluster
png(file = './plots/EC8_new/heatmap_neuron_new.png',
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

png(file = './plots/EC8_new/N_sum_by_cluste_new.png',
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


png(file = './plots/EC8_new/ppc_zero_new.png',
    width = 1200,
    height = 700)

ppc_multiple$zero.plot

dev.off()

png(file = './plots/EC8_new/ppc_nonzero_new.png',
    width = 1200,
    height = 700)

ppc_multiple$non.zero.plot

dev.off()

# Posterior predictive check with single replicated data
ppc_single <- ppc_single_f(mcmc_run_all_output = mcmc_all_EC8,
                           Y = EC8_new,
                           regions.name = rownames(EC8_new[[1]]))

for(m in 1:length(EC8_new)){
  
  png(file = paste0('./plots/EC8_new/ppc_single_mouse_new_', m, '.png'),
      width = 1200,
      height = 700)
  
  print(ppc_single[[m]])
  
  dev.off()
}


#---------------------------------------------------------------------------------------------------



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

mx <- matrix(0, nrow = 5, ncol = 5)

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


for(i in 1:5){
  
  mx[i,] <- sapply(1:5,
                   function(v) variation_info(unlist(clusterings[[i]]),
                                              unlist(clusterings[[v]])))
}

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


df$EC_label <- unlist(EC8_EC_label)




#-----------------------------------------------------------------------------------------------

# Write a function to choose 50 Bayesian motifs with the largest posterior probability of having w_j > tau

find_large_cluster <- function(mcmc_run_omega_output,
                               tau = 0.01,
                               top.n){
  
  # Trace of omega
  omega.trace <- mcmc_run_omega_output$omega_J_M_output
  
  # Number of clusters
  J <- nrow(omega.trace[[1]])
  
  # Posterior probability of having w_j > tau
  posterior.probability <- sapply(1:J,
                                  function(j){
                                    
                                    omega.trace.j <- sapply(1:length(omega.trace),
                                                            function(t) omega.trace[[t]][j])
                                    
                                    mean(omega.trace.j > tau)
                                  })
  
  return(which(rank(desc(posterior.probability)) <= top.n))
}


# top 50 large Bayesian motifs
Bayesian_motifs_large <- find_large_cluster(mcmc_run_omega_output = omega_JM_mcmc,
                                            tau = 0.02,
                                            top.n = 30)



# For each large Bayesian motif, calculate LEC/MEC proportion
for(j in Bayesian_motifs_large){
  
  df.j <- df %>%
    filter(bayesian_allocation == j)
  
  # Bayesian
  LEC_prop_bayesian <- length(which(df.j$EC_label == 'LEC'))/nrow(df.j)
  n_bayesian <- nrow(df.j)
  
  # Binomial
  LEC_binom <- df.j %>%
    group_by(binomial_allocation) %>%
    summarise(LEC = length(which(EC_label == 'LEC'))/n(),
              count = n()) %>%
    mutate(MEC = 1-LEC) %>%
    mutate(binomial_allocation = paste0('binomial motif ', binomial_allocation)) %>%
    rename(motif = binomial_allocation)
  
  # Combine both
  df.j <- rbind(data.frame(motif = paste0('bayesian motif ', j),
                           LEC = LEC_prop_bayesian,
                           count = n_bayesian,
                           MEC = 1-LEC_prop_bayesian),
                
                LEC_binom)
  
  df.j$motif <- factor(df.j$motif, levels = unique(df.j$motif))
  
  png(filename = paste0('./plots/EC8_new/large_bayesian_motifs/large_bayesian_motifs_', j, '.png'),
      width = 50*nrow(df.j))
  
  print(df.j %>%
          pivot_longer(cols = c(2,4)) %>%
          rename(EC_group = name) %>%
          ggplot()+
          geom_bar(mapping = aes(x = motif,
                                 y = value,
                                 fill = EC_group),
                   stat = 'identity')+
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
          ggtitle(paste0('bayesian motif ', j))+
          annotate('text',
                   x = unique(df.j$motif),
                   y = 0.95,
                   label = df.j$count,
                   size = 3,
                   angle = 90)+
          geom_hline(yintercept = mean(unlist(EC8_EC_label) == 'MEC'), col = "red"))
  
  dev.off()

}


for(j in Bayesian_motifs_large){
  
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
  
  png(filename = paste0('plots/EC8_new/large_bayesian_motifs/bayesian_motif_', j, '.png'),
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

#-----------------------------------------------------------------------------------------------

# Top 30 binomial motifs with the largest number of neurons

large_binomial_motifs <- df %>%
  group_by(binomial_allocation) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  top_n(30) %>%
  select(binomial_allocation) %>%
  pull() %>%
  sort()

# For each large Binomial motif, calculate LEC/MEC proportion
for(j in large_binomial_motifs){
  
  df.j <- df %>%
    filter(binomial_allocation == j)
  
  # Binomial
  LEC_prop_binomial <- length(which(df.j$EC_label == 'LEC'))/nrow(df.j)
  n_binomial <- nrow(df.j)
  
  # Bayesian
  LEC_bayes <- df.j %>%
    group_by(bayesian_allocation) %>%
    summarise(LEC = length(which(EC_label == 'LEC'))/n(),
              count = n()) %>%
    mutate(MEC = 1-LEC) %>%
    mutate(bayesian_allocation = paste0('bayesian motif ', bayesian_allocation)) %>%
    rename(motif = bayesian_allocation)
  
  # Combine both
  df.j <- rbind(data.frame(motif = paste0('binomial motif ', j),
                           LEC = LEC_prop_binomial,
                           count = n_binomial,
                           MEC = 1-LEC_prop_binomial),
                
                LEC_bayes)
  
  df.j$motif <- factor(df.j$motif, levels = unique(df.j$motif))
  
  png(filename = paste0('./plots/EC8_new/large_binomial_motifs/large_binomial_motifs_', j, '.png'))
  
  print(df.j %>%
          pivot_longer(cols = c(2,4)) %>%
          rename(EC_group = name) %>%
          ggplot()+
          geom_bar(mapping = aes(x = motif,
                                 y = value,
                                 fill = EC_group),
                   stat = 'identity')+
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
          ggtitle(paste0('binomial motif ', j))+
          annotate('text',
                   x = unique(df.j$motif),
                   y = 0.95,
                   label = df.j$count,
                   size = 3,
                   angle = 90)+
          geom_hline(yintercept = mean(unlist(EC8_EC_label) == 'MEC'), col = "red"))
  
  dev.off()
  
}



for(j in large_binomial_motifs){
  
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
  
  png(filename = paste0('plots/EC8_new/large_binomial_motifs/binomial_motif_', j, '.png'),
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


save.image("~/OneDrive - University of Edinburgh/BrainConnectivity/data/EC8/mcmc_EC8_longeriter.RData")

#-------------------------------- Create Excel sheet to summarize ------------------------------------


for(j in Bayesian_motifs_large){
  
  #-- Plot 1
  
  # Plot of estimated projection strength
  projection.strength.j <- data.frame(med = mcmc_unique_EC8$proj_prob_med[j,],
                                      lower_bound = mcmc_unique_EC8$proj_prob_lower[j,],
                                      upper_bound = mcmc_unique_EC8$proj_prob_upper[j,],
                                      region.name = factor(rownames(EC8_new[[1]]),
                                                           levels = rownames(EC8_new[[1]])))
  
  # Projection region in the Bayesian cluster
  projecting.region.j <- rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[j,] > 0.5]
  
  plot1 <- projection.strength.j %>%
    ggplot(mapping = aes(x = region.name, y = med))+
    geom_line(color = 'black', group = 1)+
    geom_point()+
    geom_errorbar(aes(ymin = lower_bound,
                      ymax = upper_bound),
                  width = 0.1)+
    ylim(c(0,1))+
    theme_bw()+
    ylab('Estimated projection probability')+
    ggtitle(paste('Cluster', j))+
    annotate('text',
             x = 'ACA',
             y = 0.8,
             label = paste0('[', knitr::combine_words(projecting.region.j, and = ""), ']'),
             size = 10)
  
  #-- Plot 2
  
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
  
  
  plot2 <- df1 %>%
    ggplot()+
    geom_line(mapping = aes(x = brain_region, y = projection_probability, color = binomial_allocation, group = neuron_index))+
    theme_bw()+
    ylab('projection probability')+
    xlab('brain region')+
    ggtitle(paste('Cluster', j))
  
  #-- Plot 3
  
  df2.0 <- t(data_EC_cbind)
  df2.0 = t(apply(df2.0, 1, function(x){return(x/sum(x))}))
  
  df2 <- data.frame(projection.probability = as.vector(t(df2.0[which(unlist(EC8_Z_reordered)==j),])),
                    region.name = factor(rownames(EC8_new[[1]]),
                                         levels = rownames(EC8_new[[1]])),
                    EC_group = rep(unlist(EC8_EC_label)[which(unlist(EC8_Z_reordered)==j)],
                                   each = R),
                    neuron = factor(rep(1:length(which(unlist(EC8_Z_reordered)==j)), 
                                        each = nrow(EC8_new[[1]])),
                                    levels = rep(1:length(which(unlist(EC8_Z_reordered)==j)))))
  
  
  # Line graph
  plot3 <- ggplot(df2)+
    geom_line(mapping = aes(x = region.name,
                            y = projection.probability,
                            colour = EC_group,
                            group = interaction(neuron, EC_group)))+
    theme_bw()+
    xlab('region')+
    ylab('projection strengths')+
    ggtitle(paste('cluster', j))
  
  #-- Plot 4
  
  
  LEC_prop_bayesian <- length(which(df_j$EC_label == 'LEC'))/nrow(df_j)
  n_bayesian <- nrow(df_j)
  
  # Binomial
  LEC_binom <- df_j %>%
    group_by(binomial_allocation) %>%
    summarise(LEC = length(which(EC_label == 'LEC'))/n(),
              count = n()) %>%
    mutate(MEC = 1-LEC) %>%
    mutate(binomial_allocation = paste0('binomial motif ', binomial_allocation)) %>%
    rename(motif = binomial_allocation)
  
  # Combine both
  df_j <- rbind(data.frame(motif = paste0('bayesian motif ', j),
                           LEC = LEC_prop_bayesian,
                           count = n_bayesian,
                           MEC = 1-LEC_prop_bayesian),
                
                LEC_binom)
  
  df_j$motif <- factor(df_j$motif, levels = unique(df_j$motif))
  
  plot4 <- df_j %>%
    pivot_longer(cols = c(2,4)) %>%
    rename(EC_group = name) %>%
    ggplot()+
    geom_bar(mapping = aes(x = motif,
                           y = value,
                           fill = EC_group),
             stat = 'identity')+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0('bayesian motif ', j))+
    annotate('text',
             x = unique(df_j$motif),
             y = 0.95,
             label = df_j$count,
             size = 6,
             angle = 90)+
    geom_hline(yintercept = mean(unlist(EC8_EC_label) == 'MEC'), col = "red")
  
  
  # Binomial clusters
  binomial_cluster <- parse_number(LEC_binom$motif)
  table00 <- data.frame(binomial_cluster = binomial_cluster,
                        projecting_region = binomial_output_reorder$cluster_label[binomial_cluster])
  
  
  png(filename = paste0('plots/EC8_new/large_bayesian_motifs/binomial_motif_', j, '.png'),
      width = 1500,
      height = 2000)
  
  
  grid.arrange(
    tableGrob(table00),
    plot1,
    plot2,
    plot3,
    plot4,
    ncol = 2,
    widths = c(1.5, 1),
    clip = FALSE
  )
  
  dev.off()
  
  
  
  
}

