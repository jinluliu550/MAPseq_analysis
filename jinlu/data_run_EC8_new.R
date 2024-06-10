EC8_new <- read.csv('./data/EC8_new/all_brains_setE_new.csv')


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
M <- 6


#---------------------------------------- Binomial results --------------------------------------------

binomial_output <- binomial_model(data = EC8_new)


binomial_output_reorder <- binom_cluster_reorder(Y = EC8_new,
                                                 binomial_output = binomial_output)



#-------------------------------- Empirical Analysis ---------------------------------------

# Number of neurons in each mouse
png(file = './plots/EC8_new2/number_of_neurons_in_each_m.png',
    width = 300,
    height = 200)

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

png(file = './plots/EC8_new2/prop_of_EC_in_each_m.png',
    width = 400,
    height = 200)

ggplot(proportion_of_EC, aes(x = mouse, y = count, fill = EC))+
  geom_bar(position = 'fill', stat = 'identity')+
  ylab('proportion of LEC and MEC')+
  theme_bw()+
  scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))

dev.off()

#-------------------------------------------------------------------------------------------

#Initialize z
data_EC_cbind <- do.call(cbind, EC8_new)

df <- t(data_EC_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))

C_cumsum <- c(0, cumsum(sapply(1:6, function(m) ncol(EC8_new[[m]]))))

# k-means
k_mean_clust_60 <- kmeans(df, 60, iter.max = 100, nstart = 25)$cluster

clust60 <- lapply(1:6,
                  function(m) k_mean_clust_60[(C_cumsum[m]+1):C_cumsum[m+1]])



mcmc_all_EC8 <- mcmc_run_all(Y = EC8_new,
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


Zmat = matrix(unlist(mcmc_all_EC8$Z_output), length(mcmc_all_EC8$Z_output), sum(C),byrow = TRUE)

# Number of occupied components
k = apply(Zmat,1,function(x){length(unique(x))})
plot(k, type = 'l')

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

length(unique(unlist(EC8_Z)))

# Plot of posterior similarity matrix
png(file = './plots/EC8_new2/heatmap_psm_1.png',
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

png(file = './plots/EC8_new2/heatmap_psm_2.png',
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


# Number of LEC and MEC neurons in each cluster
png(file = './plots/EC8_new2/number_of_neuron_by_EC.png',
    width = 2500,
    height = 700)

opt.clustering.frequency1(clustering = mcmc_unique_EC8$Z,
                          EC_label = EC8_EC_label)

dev.off()



# Proportion of LEC and MEC in each cluster
png(file = './plots/EC8_new2/proportion_of_EC.png',
    width = 2500,
    height = 700)

opt.clustering.frequency2(clustering = mcmc_unique_EC8$Z,
                          EC_label = EC8_EC_label)

dev.off()

# Number of neurons by cluster and mosue
png(file = './plots/EC8_new2/number_of_neuron_by_m.png',
    width = 2500,
    height = 700)

opt.clustering.frequency(clustering = mcmc_unique_EC8$Z)


dev.off()

opt.clustering.frequency.large(clustering = mcmc_unique_EC8$Z)

# Plot of estimated projection strength
png(file = './plots/EC8_new2/estimated_pp.png',
    width = 3500,
    height = 2000)

mcmc_unique_EC8$estimated.pp.plot

dev.off()

# q tilde
png(file = './plots/EC8_new2/q_tilde.png',
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

png(file = './plots/EC8_new2/w_jm_EC.png',
    width = 4000,
    height = 2200)

omega_JM_mcmc$omega_JM_plot

dev.off()


# Projection strength of each neuron in each cluster, color-coded by the injection site
png(file = './plots/EC8_new2/projection_by_EC.png',
    width = 4000,
    height = 2200)

projection_by_EC(Y = EC8_new,
                 EC_label = EC8_EC_label,
                 Z = mcmc_unique_EC8$Z,
                 region_name = rownames(EC8_new[[1]]))

dev.off()


# Examining the difference of omega_jm between any 2 mice:
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc)

for(m in 1:M){
  
  print(length(which(difference_omega_JM$significant_obs$data1 == m | difference_omega_JM$significant_obs$data2 == m)))
}


png(file = './plots/EC8_new2/w_jm_difference.png',
    width = 2000,
    height = 600)

difference_omega_JM$probability_plot

dev.off()

png(file = './plots/EC8_new2/w_jm_difference_eg.png',
    width = 1500,
    height = 800)

difference_in_omega_jm_plot(difference_in_omega_jm_output = difference_omega_JM,
                            mcmc_run_omega_output = omega_JM_mcmc,
                            N = 20)

dev.off()

for(m in 1:M){
  
  png(file = paste0('./plots/EC8_new2/w_jm_difference_', m, '.png'),
      width = 1500,
      height = 800)
  
  print(difference_omega_JM$plot.for.each.data[[m]])
  
  dev.off()
}

# Heatmap of projection strength of neurons in each cluster
png(file = './plots/EC8_new2/heatmap_neuron.png',
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

png(file = './plots/EC8_new2/N_sum_by_cluster.png',
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


png(file = './plots/EC8_new2/ppc_zero.png',
    width = 1200,
    height = 700)

ppc_multiple$zero.plot

dev.off()

png(file = './plots/EC8_new2/ppc_nonzero.png',
    width = 1200,
    height = 700)

ppc_multiple$non.zero.plot

dev.off()

# Posterior predictive check with single replicated data
ppc_single <- ppc_single_f(mcmc_run_all_output = mcmc_all_EC8,
                           Y = EC8_new,
                           regions.name = rownames(EC8_new[[1]]))

for(m in 1:length(EC8_new)){
  
  png(file = paste0('./plots/EC8_new2/ppc_single_mouse_', m, '.png'),
      width = 1200,
      height = 700)
  
  print(ppc_single[[m]])
  
  dev.off()
}


#---------------------------------------------------------------------------------------------------

df <- t(data_EC_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))


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
png(file = './plots/EC8_new2/k_means_30.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = EC8_new,
                     Z = clust30_r,
                     regions.name = rownames(EC8_new[[1]]))

dev.off()

png(file = './plots/EC8_new2/k_means_50.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = EC8_new,
                     Z = clust50_r,
                     regions.name = rownames(EC8_new[[1]]))

dev.off()

png(file = './plots/EC8_new2/k_means_70.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = EC8_new,
                     Z = clust70_r,
                     regions.name = rownames(EC8_new[[1]]))

dev.off()

#-------------------------------- Comparison ----------------------------------------------------

png(file = './plots/EC8_new2/heatmap_binomial.png',
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




#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------

# top 20 large Bayesian motifs
Bayesian_motifs_large <- find_large_cluster(mcmc_run_omega_output = omega_JM_mcmc,
                                            tau = 0.01,
                                            top.n = 20)


# Top 20 binomial motifs with the largest number of neurons
large_binomial_motifs <- df %>%
  group_by(binomial_allocation) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  top_n(20) %>%
  select(binomial_allocation) %>%
  pull() %>%
  sort()

# For each large Binomial motif, calculate LEC/MEC proportion
for(j in large_binomial_motifs[1:20]){
  
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
  
  png(filename = paste0('./plots/EC8_new/large_binomial_motifs_new80/large_binomial_motifs_', j, '.png'))
  
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
  
  png(filename = paste0('plots/EC8_new/large_binomial_motifs_new80/binomial_motif_', j, '.png'),
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
  
  group.colors <- c('LEC'= "#F8766D",'MEC'= "#00BFC4")
  
  # Line graph
  plot3 <- ggplot(df2)+
    geom_line(mapping = aes(x = region.name,
                            y = projection.probability,
                            colour = EC_group,
                            group = interaction(neuron, EC_group)))+
    theme_bw()+
    xlab('region')+
    ylab('projection strengths')+
    ggtitle(paste('cluster', j))+
    scale_color_manual(values=group.colors)
  
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
             size = 4,
             angle = 90)+
    geom_hline(yintercept = mean(unlist(EC8_EC_label) == 'MEC'), col = "red")+
    scale_fill_manual(values = group.colors)
  
  
  # Binomial clusters
  binomial_cluster <- parse_number(LEC_binom$motif)
  
  if(length(binomial_cluster) <= 10){
    
    
    table00 <- data.frame(binomial_cluster = binomial_cluster,
                          projecting_region = binomial_output_reorder$cluster_label[binomial_cluster])
    
    
    
  }else{
    table00 <- data.frame(binomial_cluster = binomial_cluster[1:round(length(binomial_cluster)/2)],
                          projecting_region = binomial_output_reorder$cluster_label[binomial_cluster][1:round(length(binomial_cluster)/2)])
    
    table01 <- data.frame(binomial_cluster = binomial_cluster[(round(length(binomial_cluster)/2)+1):length(binomial_cluster)],
                          projecting_region = binomial_output_reorder$cluster_label[binomial_cluster][(round(length(binomial_cluster)/2)+1):length(binomial_cluster)])
  }
  
  
  png(filename = paste0('plots/EC8_new/large_bayesian_motifs_new80/bayesian_motif_', j, '.png'),
      width = 1500,
      height = 2000)
  
  if(length(binomial_cluster) <= 10){
    
    ggempty <- ggplot() + theme_void()
    
    grid.arrange(
      tableGrob(table00),
      ggempty,
      plot1,
      plot2,
      plot3,
      plot4,
      ncol = 2,
      widths = c(1, 1),
      clip = FALSE
    )
    
    
  }else{
    grid.arrange(
      tableGrob(table00),
      tableGrob(table01),
      plot1,
      plot2,
      plot3,
      plot4,
      ncol = 2,
      widths = c(1, 1),
      clip = FALSE
    )
  }
  
  dev.off()
  
}


for(j in 1:length(unique(unlist(mcmc_unique_EC8$Z)))){
  
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
  
  group.colors <- c('LEC'= "#F8766D",'MEC'= "#00BFC4")
  
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
    ggtitle(paste('cluster', j))+
    scale_color_manual(values=group.colors)
  
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
             size = 4,
             angle = 90)+
    geom_hline(yintercept = mean(unlist(EC8_EC_label) == 'MEC'), col = "red")+
    scale_fill_manual(values = group.colors)
  
  
  # Binomial clusters
  binomial_cluster <- parse_number(LEC_binom$motif)
  
  if(length(binomial_cluster) <= 10){
    
    
    table00 <- data.frame(binomial_cluster = binomial_cluster,
                          projecting_region = binomial_output_reorder$cluster_label[binomial_cluster])
    
    
    
  }else{
    table00 <- data.frame(binomial_cluster = binomial_cluster[1:round(length(binomial_cluster)/2)],
                          projecting_region = binomial_output_reorder$cluster_label[binomial_cluster][1:round(length(binomial_cluster)/2)])
    
    table01 <- data.frame(binomial_cluster = binomial_cluster[(round(length(binomial_cluster)/2)+1):length(binomial_cluster)],
                          projecting_region = binomial_output_reorder$cluster_label[binomial_cluster][(round(length(binomial_cluster)/2)+1):length(binomial_cluster)])
  }
  
  
  png(filename = paste0('plots/EC8_new/all_bayesian_motifs_new80/bayesian_motif_', j, '.png'),
      width = 1500,
      height = 2000)
  
  if(length(binomial_cluster) <= 10){
    
    ggempty <- ggplot() + theme_void()
    
    grid.arrange(
      tableGrob(table00),
      ggempty,
      plot1,
      plot2,
      plot3,
      plot4,
      ncol = 2,
      widths = c(1, 1),
      clip = FALSE
    )
    
    
  }else{
    grid.arrange(
      tableGrob(table00),
      tableGrob(table01),
      plot1,
      plot2,
      plot3,
      plot4,
      ncol = 2,
      widths = c(1, 1),
      clip = FALSE
    )
  }
  
  dev.off()
  
  
  
}


# Number of neurons in each cluster for the large Bayesian clusters

png(filename = paste0('plots/EC8_new/large_cluster_J80/number_of_neurons_by_EC_new.png'),
    width = 1500,
    height = 700)

opt.clustering.frequency1_large(clustering = mcmc_unique_EC8$Z,
                                EC_label = EC8_EC_label,
                                large_cluster_index = Bayesian_motifs_large)

dev.off()


png(filename = paste0('plots/EC8_new/large_cluster_J80/number_of_neurons_new.png'),
    width = 1500,
    height = 700)

opt.clustering.frequency_large(clustering = mcmc_unique_EC8$Z,
                               large_cluster_index = Bayesian_motifs_large)

dev.off()

png(filename = paste0('plots/EC8_new/large_cluster_J80/proportion_of_EC_new.png'),
    width = 1500,
    height = 700)

opt.clustering.frequency2_large(clustering = mcmc_unique_EC8$Z,
                                EC_label = EC8_EC_label,
                                large_cluster_index = Bayesian_motifs_large)

dev.off()


png(filename = paste0('plots/EC8_new/large_cluster_J80/estimated_pp_new.png'),
    width = 1500,
    height = 1000)

estimated_projection_strength_large(mcmc_run_unique_output = mcmc_unique_EC8,
                                    large_cluster_index = Bayesian_motifs_large,
                                    region.name = rownames(EC8_new[[1]]))

dev.off()


png(filename = paste0('plots/EC8_new/large_cluster_J80/q_tilde_new.png'),
    width = 1500,
    height = 700)

q_tilde_plot_large(mcmc_run_unique_output = mcmc_unique_EC8,
                   large_cluster_index = Bayesian_motifs_large,
                   region.name = rownames(EC8_new[[1]]))

dev.off()


# Gel plot

gel_plot_EC <- gel_plot(Y = EC8_new)

ggarrange(gel_plot_EC[[1]], 
          gel_plot_EC[[2]],
          gel_plot_EC[[3]],
          gel_plot_EC[[4]],
          gel_plot_EC[[5]],
          gel_plot_EC[[6]])


# Scatter plot


# VIS and SS
region1 <- 'VIS'
region2 <- 'SS'

list_region1 <- lapply(1:length(Y),
                       function(m){
                         
                         
                         Y_m <- Y[[m]]
                         
                         Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
                         
                         data.frame(mouse = m,
                                    region = region1,
                                    projection_strength = as.vector(Y_m_prop[region1,]),
                                    neuron = factor(1:ncol(Y_m)))
                         
                       })

df_region1 <- do.call(rbind, list_region1)

# region 2
list_region2 <- lapply(1:length(Y),
                       function(m){
                         
                         
                         Y_m <- Y[[m]]
                         
                         Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
                         
                         data.frame(mouse = m,
                                    region = region2,
                                    projection_strength = as.vector(Y_m_prop[region2,]),
                                    neuron = factor(1:ncol(Y_m)))
                         
                       })

df_region2 <- do.call(rbind, list_region2)


df <- rbind(df_region1,
            df_region2)

df$mouse <- factor(paste('mouse', df$mouse),
                   levels = paste('mouse', 1:length(Y)))

# Histogram
histogram <- ggplot(df, aes(x=projection_strength))+
  geom_histogram(binwidth = 0.05)+
  facet_grid(vars(region), vars(mouse))+
  xlab('y/N')+
  theme_bw()

# Scatter plot
scatter_plot <- df %>%
  pivot_wider(names_from = region,
              values_from = projection_strength) %>%
  ggplot(mapping = aes(x = VIS, y = SS))+
  geom_point(mapping = aes(color = mouse))+
  theme_bw()+
  xlim(c(0,1))+
  ylim(c(0,1))

png(filename = 'plots/EC8_new/VIS_SS_histogram.png',
    width = 1200,
    height = 400)

histogram

dev.off()


png(filename = 'plots/EC8_new/VIS_SS_scatterplot.png',
    width = 500,
    height = 500)

scatter_plot

dev.off()

# ACA and PTL
region1 <- 'ACA'
region2 <- 'PTL'

list_region1 <- lapply(1:length(Y),
                       function(m){
                         
                         
                         Y_m <- Y[[m]]
                         
                         Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
                         
                         data.frame(mouse = m,
                                    region = region1,
                                    projection_strength = as.vector(Y_m_prop[region1,]),
                                    neuron = factor(1:ncol(Y_m)))
                         
                       })

df_region1 <- do.call(rbind, list_region1)

# region 2
list_region2 <- lapply(1:length(Y),
                       function(m){
                         
                         
                         Y_m <- Y[[m]]
                         
                         Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
                         
                         data.frame(mouse = m,
                                    region = region2,
                                    projection_strength = as.vector(Y_m_prop[region2,]),
                                    neuron = factor(1:ncol(Y_m)))
                         
                       })

df_region2 <- do.call(rbind, list_region2)


df <- rbind(df_region1,
            df_region2)

df$mouse <- factor(paste('mouse', df$mouse),
                   levels = paste('mouse', 1:length(Y)))

# Histogram
histogram <- ggplot(df, aes(x=projection_strength))+
  geom_histogram(binwidth = 0.05)+
  facet_grid(vars(region), vars(mouse))+
  xlab('y/N')+
  theme_bw()

# Scatter plot
scatter_plot <- df %>%
  pivot_wider(names_from = region,
              values_from = projection_strength) %>%
  ggplot(mapping = aes(x = ACA, y = PTL))+
  geom_point(mapping = aes(color = mouse))+
  theme_bw()+
  xlim(c(0,1))+
  ylim(c(0,1))

png(filename = 'plots/EC8_new/ACA_PTL_histogram.png',
    width = 1200,
    height = 400)

histogram

dev.off()


png(filename = 'plots/EC8_new/ACA_PTL_scatterplot.png',
    width = 500,
    height = 500)

scatter_plot

dev.off()


# PFC and MO
region1 <- 'PFC'
region2 <- 'MO'

list_region1 <- lapply(1:length(Y),
                       function(m){
                         
                         
                         Y_m <- Y[[m]]
                         
                         Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
                         
                         data.frame(mouse = m,
                                    region = region1,
                                    projection_strength = as.vector(Y_m_prop[region1,]),
                                    neuron = factor(1:ncol(Y_m)))
                         
                       })

df_region1 <- do.call(rbind, list_region1)

# region 2
list_region2 <- lapply(1:length(Y),
                       function(m){
                         
                         
                         Y_m <- Y[[m]]
                         
                         Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
                         
                         data.frame(mouse = m,
                                    region = region2,
                                    projection_strength = as.vector(Y_m_prop[region2,]),
                                    neuron = factor(1:ncol(Y_m)))
                         
                       })

df_region2 <- do.call(rbind, list_region2)


df <- rbind(df_region1,
            df_region2)

df$mouse <- factor(paste('mouse', df$mouse),
                   levels = paste('mouse', 1:length(Y)))

# Histogram
histogram <- ggplot(df, aes(x=projection_strength))+
  geom_histogram(binwidth = 0.05)+
  facet_grid(vars(region), vars(mouse))+
  xlab('y/N')+
  theme_bw()

# Scatter plot
scatter_plot <- df %>%
  pivot_wider(names_from = region,
              values_from = projection_strength) %>%
  ggplot(mapping = aes(x = PFC, y = MO))+
  geom_point(mapping = aes(color = mouse))+
  theme_bw()+
  xlim(c(0,1))+
  ylim(c(0,1))

png(filename = 'plots/EC8_new/PFC_MO_histogram.png',
    width = 1200,
    height = 400)

histogram

dev.off()


png(filename = 'plots/EC8_new/PFC_MO_scatterplot.png',
    width = 500,
    height = 500)

scatter_plot

dev.off()


# ORB and RSC
region1 <- 'ORB'
region2 <- 'RSC'

list_region1 <- lapply(1:length(Y),
                       function(m){
                         
                         
                         Y_m <- Y[[m]]
                         
                         Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
                         
                         data.frame(mouse = m,
                                    region = region1,
                                    projection_strength = as.vector(Y_m_prop[region1,]),
                                    neuron = factor(1:ncol(Y_m)))
                         
                       })

df_region1 <- do.call(rbind, list_region1)

# region 2
list_region2 <- lapply(1:length(Y),
                       function(m){
                         
                         
                         Y_m <- Y[[m]]
                         
                         Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
                         
                         data.frame(mouse = m,
                                    region = region2,
                                    projection_strength = as.vector(Y_m_prop[region2,]),
                                    neuron = factor(1:ncol(Y_m)))
                         
                       })

df_region2 <- do.call(rbind, list_region2)


df <- rbind(df_region1,
            df_region2)

df$mouse <- factor(paste('mouse', df$mouse),
                   levels = paste('mouse', 1:length(Y)))

# Histogram
histogram <- ggplot(df, aes(x=projection_strength))+
  geom_histogram(binwidth = 0.05)+
  facet_grid(vars(region), vars(mouse))+
  xlab('y/N')+
  theme_bw()

# Scatter plot
scatter_plot <- df %>%
  pivot_wider(names_from = region,
              values_from = projection_strength) %>%
  ggplot(mapping = aes(x = ORB, y = RSC))+
  geom_point(mapping = aes(color = mouse))+
  theme_bw()+
  xlim(c(0,1))+
  ylim(c(0,1))

png(filename = 'plots/EC8_new/ORB_RSC_histogram.png',
    width = 1200,
    height = 400)

histogram

dev.off()


png(filename = 'plots/EC8_new/ORB_RSC_scatterplot.png',
    width = 500,
    height = 500)

scatter_plot

dev.off()
