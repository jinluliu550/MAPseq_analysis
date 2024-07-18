
df <- t(data_Hans_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))

# large Bayesian motifs
bayesian_motif_large <- c(7,20,22)

# significant Binomial motifs
hans_binomial_reorder$cluster_summary$cluster_index <- 1:nrow(hans_binomial_reorder$cluster_summary)

binomial_motif_large <- hans_binomial_reorder$cluster_summary %>%
  filter(significant == 'significant') %>%
  select(cluster_index) %>%
  pull()


# Large Bayesian motifs 
bayesian_motif_pp <- lapply(bayesian_motif_large,
                            function(j){
                              
                              data.frame(cluster = j,
                                         pp = as.vector(t(df[which(unlist(mcmc_unique_hans$Z) == j),])),
                                         region.name = rownames(data_Hans[[1]]),
                                         cell.index = rep(which(unlist(mcmc_unique_hans$Z) == j), each = 6),
                                         mouse.index = rep(mouse.index[which(unlist(mcmc_unique_hans$Z) == j)], each = 6))
                            })

bayesian_motif_pp <- do.call(rbind, bayesian_motif_pp)


bayesian_motif_pp %>%
  mutate(cell.index = as.factor(cell.index),
         mouse.index = as.factor(mouse.index),
         cluster = as.factor(cluster),
         region.name = factor(region.name, levels = rownames(data_Hans[[1]]))) %>%
  ggplot(mapping = aes(x = region.name,
                       y = pp,
                       color = mouse.index,
                       group = cell.index))+
  geom_line()+
  facet_wrap(~cluster, ncol = 3)+
  theme_bw()+
  xlab('region')+
  ylab('projection strength')+
  guides(color = guide_legend(title="mouse"))

# Significant Binomial motifs
binomial_motif_pp <- lapply(binomial_motif_large,
                            function(j){
                              
                              data.frame(cluster = j,
                                         pp = as.vector(t(df[which(unlist(hans_binomial_reorder$allocation) == j),])),
                                         region.name = rownames(data_Hans[[1]]),
                                         cell.index = rep(which(unlist(hans_binomial_reorder$allocation) == j), each = 6),
                                         mouse.index = rep(mouse.index[which(unlist(hans_binomial_reorder$allocation) == j)], each = 6))
                            })

binomial_motif_pp <- do.call(rbind, binomial_motif_pp)


binomial_motif_pp %>%
  mutate(cell.index = as.factor(cell.index),
         mouse.index = as.factor(mouse.index),
         cluster = as.factor(cluster),
         region.name = factor(region.name, levels = rownames(data_Hans[[1]]))) %>%
  ggplot(mapping = aes(x = region.name,
                       y = pp,
                       color = mouse.index,
                       group = cell.index))+
  geom_line()+
  facet_wrap(~cluster, ncol = 3)+
  theme_bw()+
  xlab('region')+
  ylab('projection strength')+
  guides(color = guide_legend(title="mouse"))

# Large motifs with k-means
large_k_mean_cluster <- which(as.vector(table(unlist(clust30_r))) >= 20)


# Large K-means
k_means_motif_pp <- lapply(large_k_mean_cluster,
                            function(j){
                              
                              data.frame(cluster = j,
                                         pp = as.vector(t(df[which(unlist(clust30_r) == j),])),
                                         region.name = rownames(data_Hans[[1]]),
                                         cell.index = rep(which(unlist(clust30_r) == j), each = 6),
                                         mouse.index = rep(mouse.index[which(unlist(clust30_r) == j)], each = 6))
                            })

k_means_motif_pp <- do.call(rbind, k_means_motif_pp)


k_means_motif_pp %>%
  mutate(cell.index = as.factor(cell.index),
         mouse.index = as.factor(mouse.index),
         cluster = as.factor(cluster),
         region.name = factor(region.name, levels = rownames(data_Hans[[1]]))) %>%
  ggplot(mapping = aes(x = region.name,
                       y = pp,
                       color = mouse.index,
                       group = cell.index))+
  geom_line()+
  facet_wrap(~cluster, ncol = 3)+
  theme_bw()+
  xlab('region')+
  ylab('projection strength')+
  guides(color = guide_legend(title="mouse"))