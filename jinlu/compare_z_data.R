# Code for comparing allocation obtained from different sources

z_binomial_LEC <- read.csv('./data/motif_cell_assignments_lec.csv')
z_binomial_MEC <- read.csv('./data/motif_cell_assignments_mec.csv')


# Total number of motifs
motifs_union <- union(z_binomial_LEC$Motif,
                      z_binomial_MEC$Motif)

# Name of brain regions
region_names <- rownames(data_by_mouse[[1]])

#----------------------------------- A summary of motifs from the binomial model ---------------------------------------

motifs_union_summary <- lapply(1:length(motifs_union),
                               function(i){
                                 
                                 data.frame(motif_label = i,
                                            region_name = region_names,
                                            projecting = sapply(region_names,
                                                                function(region_name) str_detect(motifs_union[i], region_name))
                                 )
                               })

motifs_union_summary <- do.call(rbind, motifs_union_summary)

motifs_union_summary$projecting <- ifelse(motifs_union_summary$projecting == TRUE, 1, 0)

motifs_union_summary$motif_label <- factor(motifs_union_summary$motif_label,
                                           levels = 1:length(motifs_union))

motifs_union_summary$region_name <- factor(motifs_union_summary$region_name,
                                           levels = region_names)


# Index of motif in LEC
motifs_LEC_index <- sapply(1:length(z_binomial_LEC$Motif), 
                           function(i) which(motifs_union == z_binomial_LEC$Motif[i]))

# Index of motif in MEC
motifs_MEC_index <- sapply(1:length(z_binomial_MEC$Motif), 
                           function(i) which(motifs_union == z_binomial_MEC$Motif[i]))

# Allocation of LEC
LEC_allocation <- data.frame(neuron = paste('LEC neuron', 1:nrow(LEC)),
                             allocation_binom = 0)

for(i in 1:length(z_binomial_LEC$Motif)){
  
  # Convert the output form to a number string - easier to work with
  label_string <- strsplit(as.character(z_binomial_LEC$cell_ids[i]), split = " ")[[1]]
  
  label_string <- as.numeric(gsub(".*?([0-9]+).*", "\\1", label_string))
  
  LEC_allocation$allocation_binom[label_string+1] <- motifs_LEC_index[i]
}

# Allocation of MEC
MEC_allocation <- data.frame(neuron = paste('MEC neuron', 1:nrow(MEC)),
                             allocation_binom = 0)

for(i in 1:length(z_binomial_MEC$Motif)){
  
  # Convert the output form to a number string - easier to work with
  label_string <- strsplit(as.character(z_binomial_MEC$cell_ids[i]), split = " ")[[1]]
  
  label_string <- as.numeric(gsub(".*?([0-9]+).*", "\\1", label_string))
  
  MEC_allocation$allocation_binom[label_string+1] <- motifs_MEC_index[i]
}

# Binomial Motifs
binomial_allocation <- rbind(LEC_allocation, MEC_allocation)

# Reordered Binomial Motifs
df <- binom_cluster_reorder(Y = data,
                            Z = binomial_allocation$allocation_binom)

binomial_allocation$allocation_binom <- df$Z

# Old cluster index
old_cluster_index <- df$old_ordering


# A combined data frame
df_combined <- merge(allocation_bayesian,
                     binomial_allocation,
                     by = 'neuron')

df_combined <- df_combined[sapply(1:nrow(df_combined), 
                                  function(i) which(df_combined$neuron == binomial_allocation$neuron[i])),]


df_combined_summary <- df_combined %>%
  group_by(allocation_bayes,
           allocation_binom) %>%
  summarise(count = n())



df_combined_summary$allocation_bayes <- factor(df_combined_summary$allocation_bayes,
                                               levels = 1:max(df_combined_summary$allocation_bayes))

df_combined_summary$allocation_binom <- factor(df_combined_summary$allocation_binom,
                                               levels = 1:max(df_combined_summary$allocation_binom))


df_combined_summary %>%
  ggplot(mapping = aes(x = allocation_binom,
                       y = allocation_bayes,
                       fill = count))+
  geom_tile(color = "gray")+
  theme_bw()+
  scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                       values = scales::rescale(c(0, 3, 330)),
                       limits = c(0,330))+
  xlab('binomial motif')+
  ylab('bayesian motif')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# Heat-map of reordered clusters in the binomial model
list00 <- lapply(1:length(motifs_union),
                 function(i) motifs_union_summary[((i-1)*8 + 1):(i*8),])

list0 <- lapply(old_cluster_index,
                function(i) list00[[i]])

df0 <- do.call(rbind, list0)
df0$motif_label <- rep(1:length(motifs_union), 
                       each = nrow(data_by_mouse[[1]]))

df0 %>%
  ggplot(mapping = aes(x = motif_label,
                       y = region_name,
                       fill = projecting))+
  geom_tile(color = "gray")+
  theme_bw()+
  scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                       values = scales::rescale(c(0, 0.25, 1)),
                       limits = c(0,1))+
  xlab('cluster')+
  ylab('regions')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")



#---------------------- For larger clusters found in the Bayesian model, color-code by the binomial allocation

for(j in c(1,6,21,53,55,63,97)){
  
  
  png(filename = paste0('bayesian_motif_', j, '.png'),
      width = 800,
      height = 500)
  
  print(df_combined %>%
          filter(allocation_bayes == j) %>%
          left_join(data3, by = 'neuron') %>%
          pivot_longer(4:11, names_to = 'brain_region') %>%
          mutate(brain_region = factor(brain_region,
                                       levels = rownames(data_by_mouse[[1]])),
                 allocation_binom = factor(allocation_binom)) %>%
          ggplot()+
          geom_line(mapping = aes(x = brain_region, y = value, color = allocation_binom, group = neuron))+
          theme_bw()+
          ylab('projection probability')+
          xlab('brain region')+
          ggtitle(paste('Cluster', j)))
  
  dev.off()
}


#-------------------------- Large motifs in the binomial model ------------------

large_binomial_motif <- df_combined %>%
  group_by(allocation_binom) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  filter(count > 100)

for(j in large_binomial_motif$allocation_binom){
  
  png(filename = paste0('binom_motif_', j, '.png'),
      width = 800,
      height = 500)
  print(df_combined %>%
          filter(allocation_binom == j) %>%
          left_join(data3, by = 'neuron') %>%
          pivot_longer(4:11, names_to = 'brain_region') %>%
          mutate(brain_region = factor(brain_region,
                                       levels = rownames(data_by_mouse[[1]])),
                 allocation_bayes = factor(allocation_bayes)) %>%
          ggplot()+
          geom_line(mapping = aes(x = brain_region, y = value, color = allocation_bayes, group = neuron))+
          theme_bw()+
          ylab('projection probability')+
          xlab('brain region')+
          ggtitle(paste('cluster', j)))
  dev.off()
}


