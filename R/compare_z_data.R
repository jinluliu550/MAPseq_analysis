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

motifs_union_summary %>%
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

binomial_allocation <- rbind(LEC_allocation, MEC_allocation)


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



# To show percentages
df_combined_summary %>%
  group_by(allocation_bayes) %>%
  summarize(percentage = count/sum(count),
            allocation_binom = allocation_binom) %>%
  
  ggplot(mapping = aes(x = allocation_binom,
                       y = allocation_bayes,
                       fill = percentage))+
  geom_tile(color = "gray")+
  theme_bw()+
  scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                       values = scales::rescale(c(0, 0.05, 1)),
                       limits = c(0,1))+
  xlab('binomial motif')+
  ylab('bayesian motif')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
