







# Bootstrapping
B <- 1000


figure_1f_df_bootstrap <- lapply(1:B,
                                 function(b){
                                   
                                   df1 <- lapply(2:5,
                                                 function(m){
                                                   
                                                   # sample index
                                                   sample.index <- sample(1:C[m], size = C[m], replace = TRUE)
                                                   
                                                   # Bootstrapping
                                                   EC8_new_m_bootstrap <- EC8_new[[m]][,sample.index]
                                                   EC8_EC_m_bootstrap <- EC8_EC_label[[m]][sample.index]
                                                   
                                                   # Overall MEC proportion in mouse m
                                                   MEC_proportion_m <- mean(EC8_EC_m_bootstrap == 'MEC')
                                                   
                                                   df0 <- lapply(1:R,
                                                                 function(r){
                                                                   
                                                                   # index of neuron projecting to region r
                                                                   neuron.index.project.to.r <- which(EC8_new_m_bootstrap[r,] != 0)
                                                                   
                                                                   # Proportion of MEC neurons
                                                                   MEC.proportion <- length(which(EC8_EC_m_bootstrap[neuron.index.project.to.r] == 'MEC'))/
                                                                     length(neuron.index.project.to.r)
                                                                   
                                                                   
                                                                   data.frame(bootstrap.num = b,
                                                                              mouse = m,
                                                                              region = regions.name[r],
                                                                              conditional_probability_diff = MEC.proportion - MEC_proportion_m)
                                                                 })
                                                   
                                                   do.call(rbind, df0)
                                                   
                                                 })
                                   
                                   do.call(rbind, df1)
                                   
                                 })

figure_1f_df_bootstrap <- do.call(rbind, figure_1f_df_bootstrap)


figure_1f_df_bootstrap_summary <- figure_1f_df_bootstrap %>%
  group_by(mouse, region) %>%
  summarise(median = mean(conditional_probability_diff),
            diff_lb = mean(conditional_probability_diff)  - qt(0.9875,df=100-1)*sqrt(var(conditional_probability_diff))/sqrt(100),
            diff_ub = mean(conditional_probability_diff)  + qt(0.9875,df=100-1)*sqrt(var(conditional_probability_diff))/sqrt(100))

figure_1f_df_bootstrap_summary$mouse <- factor(figure_1f_df_bootstrap_summary$mouse, levels = 2:5)
figure_1f_df_bootstrap_summary$region <- factor(figure_1f_df_bootstrap_summary$region, levels = rownames(EC8_new[[1]]))

ggplot(figure_1f_df_bootstrap_summary,
       aes(x = region,
           y = median,
           fill = mouse))+
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=diff_lb,
                    ymax = diff_ub), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  ylab('P(from MEC|project to r, mouse m) - P(from MEC|mouse m)')


# One sided test
PFC_ORB_p_value <- figure_1f_df_bootstrap %>%
  filter(region %in% c('PFC', 'ORB')) %>%
  mutate(region = factor(region, levels = rownames(EC8_new[[1]]))) %>%
  group_by(mouse, region) %>%
  summarise(p_value_one_sided = mean(conditional_probability_diff > 0)) %>%
  mutate(significant = ifelse(p_value_one_sided < 0.05/4, 'significant', 'insignificant'))

all_other_region_p_value <- figure_1f_df_bootstrap %>%
  filter(region %notin% c('PFC', 'ORB')) %>%
  mutate(region = factor(region, levels = rownames(EC8_new[[1]]))) %>%
  group_by(mouse, region) %>%
  summarise(p_value_one_sided = mean(conditional_probability_diff < 0)) %>%
  mutate(significant = ifelse(p_value_one_sided < 0.05/4, 'significant', 'insignificant'))


p_value <- rbind(PFC_ORB_p_value,
                 all_other_region_p_value) %>%
  arrange(mouse, region)


ggplot(figure_1f_df_bootstrap_summary)+
  
  geom_bar(stat="identity", color="black", 
           position=position_dodge(),
           
           aes(x = region,
               y = median,
               fill = mouse)) +
  
  geom_errorbar(aes(x = region,
                    y = median,
                    group = mouse,
                    
                    ymin=diff_lb,
                    ymax = diff_ub), width=.2,
                position=position_dodge(.9))+
  
  theme_bw()+
  
  ylab('P(from MEC|project to r, mouse m) - P(from MEC|mouse m)')+
  
  geom_point(data = p_value,
             aes(x = region,
                 y = 0.3,
                 group = mouse,
                 color = significant),
             stat = 'identity',
             position = position_dodge(width = 0.9),
             color = ifelse(p_value$significant == 'significant', 'black', 'white'),
             shape = 8,
             show.legend = FALSE)+
  
  guides(shape = "none")+
  
  ylim(c(-0.3, 0.35))


# Figure 3b


# Cluster 1, 6, 14, 5
figure_3b_cluster <- c(1,6,14,5)

projecting.region.name <- sapply(1:4, 
                                 function(i) rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[figure_3b_cluster[i],] > 0.5])

df_text <- data.frame(label  = sapply(1:4,
                                      function(i){
                                        
                                        paste0(knitr::combine_words(projecting.region.name[[i]], and = ""))
                                      }),
                      
                      cluster = paste('cluster', figure_3b_cluster))

df2 <- lapply(figure_3b_cluster,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_med[j,],
                           lower_projection_strength = mcmc_unique_EC8$proj_prob_lower[j,],
                           upper_projection_strength = mcmc_unique_EC8$proj_prob_upper[j,])
              })

df2 <- do.call(rbind, df2)
df2$region <- factor(df2$region,
                     levels = rownames(EC8_new[[1]]))



LEC_proportion_overall <- table(unlist(EC8_EC_label))[1]/length(unlist(EC8_EC_label))
MEC_proportion_overall <- 1-LEC_proportion_overall

# Proportion of LEC and MEC
figure_3b_proportion <- lapply(figure_3b_cluster,
                               function(j){
                                 
                                 LEC_proportion <- table(unlist(EC8_EC_label)[which(unlist(mcmc_unique_EC8$Z) == j)])[1]/length(which(unlist(mcmc_unique_EC8$Z) == j))
                                 MEC_proportion <- 1-LEC_proportion
                                 
                                 data.frame(cluster = j,
                                            EC = c('LEC','MEC'),
                                            proportion = c(LEC_proportion, MEC_proportion))
                               })

figure_3b_proportion <- do.call(rbind, figure_3b_proportion)



# For each cluster
plot1 <- lapply(figure_3b_cluster,
                function(j){
                  
                  df2 %>%
                    filter(cluster == paste('cluster', j)) %>%
                    ggplot(mapping = aes(x = region,
                                         y = mean_projection_strength,
                                         group = 1))+
                    geom_line()+
                    geom_errorbar(mapping = aes(ymin = lower_projection_strength, 
                                                ymax = upper_projection_strength),
                                  width = 0.1)+
                    geom_point()+
                    theme_bw()+
                    ylab('projection probability')+
                    # geom_text(data = df_text[which(df_text$cluster == paste('cluster', j)),],
                    #           mapping = aes(x = 'VIS', y = 0.6, label = paste(strsplit(label, ', ')[[1]], collapse = '\n')),
                    #           hjust   = -0.1)+
                    theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          axis.text = element_text(size = 10,
                                                   color = 'black')
                          )+
                    scale_y_continuous(breaks = c(0,1), limits = c(0,1))
                })


plot2 <- lapply(figure_3b_cluster,
                function(j){
                  
                  figure_3b_proportion %>%
                    filter(cluster == j) %>%
                    ggplot(aes(x = factor(cluster, levels = figure_3b_cluster),
                               y = proportion,
                               fill = EC))+
                    geom_bar(stat = 'identity')+
                    theme_bw()+
                    xlab('')+
                    ylab('')+
                    scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                    geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                               color = "red",
                               linewidth = 1.5)+
                    coord_flip()+
                    ggtitle(paste0((length(which(unlist(mcmc_unique_EC8$Z)==j))), ' neurons, ',
                                   length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'LEC']==j)), ' LEC, ',
                                   length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'MEC']==j)), ' MEC'))+
                    theme(legend.title=element_blank())+
                    theme(legend.position="none")+
                    xlab('')+
                    scale_y_continuous(breaks = c(0,1))+
                    theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          plot.title = element_text(size=12))+
                    theme(axis.text = element_text(size = 12)) 
                  
                })

plot_all <- lapply(1:4, function(i) ggarrange(plot2[[i]], plot1[[i]], nrow = 2, heights = c(0.45,1)))


png(file = './plots/EC8_new2/figure_3b.png',
    width = 900,
    height = 900/4)

ggarrange(plot_all[[1]],
          plot_all[[2]],
          plot_all[[3]],
          plot_all[[4]],
          nrow = 1)

dev.off()





# Figure 3c
figure_3c_cluster <- c(35,25,41,53)



projecting.region.name <- sapply(1:4, 
                                 function(i) rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[figure_3c_cluster[i],] > 0.5])

df_text <- data.frame(label  = sapply(1:4,
                                      function(i){
                                        
                                        paste0(knitr::combine_words(projecting.region.name[[i]], and = ""))
                                      }),
                      
                      cluster = paste('cluster', figure_3c_cluster))

df2 <- lapply(figure_3c_cluster,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_med[j,],
                           lower_projection_strength = mcmc_unique_EC8$proj_prob_lower[j,],
                           upper_projection_strength = mcmc_unique_EC8$proj_prob_upper[j,])
              })

df2 <- do.call(rbind, df2)
df2$region <- factor(df2$region,
                     levels = rownames(EC8_new[[1]]))




# Proportion of LEC and MEC
figure_3c_proportion <- lapply(figure_3c_cluster,
                               function(j){
                                 
                                 LEC_proportion <- length(which(unlist(EC8_EC_label)[which(unlist(mcmc_unique_EC8$Z) == j)] == 'LEC'))/
                                   length(which(unlist(mcmc_unique_EC8$Z) == j))
                                 
                                 MEC_proportion <- 1-LEC_proportion
                                 
                                 data.frame(cluster = j,
                                            EC = c('LEC','MEC'),
                                            proportion = c(LEC_proportion, MEC_proportion))
                               })

figure_3c_proportion <- do.call(rbind, figure_3c_proportion)



# For each cluster
plot1 <- lapply(figure_3c_cluster,
                function(j){
                  
                  df2 %>%
                    filter(cluster == paste('cluster', j)) %>%
                    ggplot(mapping = aes(x = region,
                                         y = mean_projection_strength,
                                         group = 1))+
                    geom_line()+
                    geom_errorbar(mapping = aes(ymin = lower_projection_strength, 
                                                ymax = upper_projection_strength),
                                  width = 0.1)+
                    geom_point()+
                    theme_bw()+
                    ylab('projection probability')+
                    # geom_text(data = df_text[which(df_text$cluster == paste('cluster', j)),],
                    #           mapping = aes(x = 'VIS', y = 0.5, label = paste(strsplit(label, ', ')[[1]], collapse = '\n')),
                    #           hjust   = -0.1)+
                    theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          axis.text = element_text(size = 10,
                                                   color = 'black'))+
                    scale_y_continuous(breaks = c(0,1), limits = c(0,1))
                })

plot2 <- lapply(figure_3c_cluster,
                function(j){
                  
                  figure_3c_proportion %>%
                    filter(cluster == j) %>%
                    ggplot(aes(x = factor(cluster, levels = figure_3c_cluster),
                               y = proportion,
                               fill = EC))+
                    geom_bar(stat = 'identity')+
                    theme_bw()+
                    xlab('')+
                    ylab('')+
                    scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                    geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                               color = "red",
                               linewidth = 1.5)+
                    coord_flip()+
                    ggtitle(paste0((length(which(unlist(mcmc_unique_EC8$Z)==j))), ' neurons, ',
                                   length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'LEC']==j)), ' LEC, ',
                                   length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'MEC']==j)), ' MEC'))+
                    theme(legend.title=element_blank())+
                    theme(legend.position="none")+
                    xlab('')+
                    scale_y_continuous(breaks = c(0,1))+
                    theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          plot.title = element_text(size=12))+
                    theme(axis.text = element_text(size = 12)) 
                  
                })

plot_all <- lapply(1:4, function(i) ggarrange(plot2[[i]], plot1[[i]], nrow = 2, heights = c(0.45,1)))


png(file = './plots/EC8_new2/figure_3c.png',
    width = 900,
    height = 900/4)

ggarrange(plot_all[[1]],
          plot_all[[2]],
          plot_all[[3]],
          plot_all[[4]],
          nrow = 1)

dev.off()


# Figure 3d
figure_3d_cluster <- c(48,54)



projecting.region.name <- sapply(1:2, 
                                 function(i) rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[figure_3d_cluster[i],] > 0.5])

df_text <- data.frame(label  = sapply(1:2,
                                      function(i){
                                        
                                        paste0(knitr::combine_words(projecting.region.name[[i]], and = ""))
                                      }),
                      
                      cluster = paste('cluster', figure_3d_cluster))

df2 <- lapply(figure_3d_cluster,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_med[j,],
                           lower_projection_strength = mcmc_unique_EC8$proj_prob_lower[j,],
                           upper_projection_strength = mcmc_unique_EC8$proj_prob_upper[j,])
              })

df2 <- do.call(rbind, df2)
df2$region <- factor(df2$region,
                     levels = rownames(EC8_new[[1]]))




# Proportion of LEC and MEC
figure_3d_proportion <- lapply(figure_3d_cluster,
                               function(j){
                                 
                                 LEC_proportion <- length(which(unlist(EC8_EC_label)[which(unlist(mcmc_unique_EC8$Z) == j)] == 'LEC'))/
                                   length(which(unlist(mcmc_unique_EC8$Z) == j))
                                 
                                 MEC_proportion <- 1-LEC_proportion
                                 
                                 data.frame(cluster = j,
                                            EC = c('LEC','MEC'),
                                            proportion = c(LEC_proportion, MEC_proportion))
                               })

figure_3d_proportion <- do.call(rbind, figure_3d_proportion)



# For each cluster
plot1 <- lapply(figure_3d_cluster,
                function(j){
                  
                  df2 %>%
                    filter(cluster == paste('cluster', j)) %>%
                    ggplot(mapping = aes(x = region,
                                         y = mean_projection_strength,
                                         group = 1))+
                    geom_line()+
                    geom_errorbar(mapping = aes(ymin = lower_projection_strength, 
                                                ymax = upper_projection_strength),
                                  width = 0.1)+
                    geom_point()+
                    theme_bw()+
                    ylab('projection probability')+
                    # geom_text(data = df_text[which(df_text$cluster == paste('cluster', j)),],
                    #           mapping = aes(x = 'VIS', y = 0.5, label = paste(strsplit(label, ', ')[[1]], collapse = '\n')),
                    #           hjust   = -0.1)+
                    theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          axis.text = element_text(size = 10,
                                                   color = 'black'))+
                    scale_y_continuous(breaks = c(0,1), limits = c(0,1))
                })



plot2 <- lapply(figure_3d_cluster,
                function(j){
                  
                  
                  figure_3d_proportion %>%
                    filter(cluster == j) %>%
                    ggplot(aes(x = factor(cluster, levels = figure_3d_cluster),
                               y = proportion,
                               fill = EC))+
                    geom_bar(stat = 'identity')+
                    theme_bw()+
                    xlab('')+
                    ylab('')+
                    scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                    geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                               color = "red",
                               linewidth = 1.5)+
                    coord_flip()+
                    ggtitle(paste0((length(which(unlist(mcmc_unique_EC8$Z)==j))), ' neurons, ',
                                   length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'LEC']==j)), ' LEC, ',
                                   length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'MEC']==j)), ' MEC'))+
                    theme(legend.title=element_blank())+
                    theme(legend.position="none")+
                    xlab('')+
                    scale_y_continuous(breaks = c(0,1))+
                    theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          plot.title = element_text(size=12))+
                    theme(axis.text = element_text(size = 12)) 
                  
                })

plot_all <- lapply(1:2, function(i) ggarrange(plot2[[i]], plot1[[i]], nrow = 2, heights = c(0.45,1)))


png(file = './plots/EC8_new2/figure_3d.png',
    width = 450,
    height = 450/2)

ggarrange(plot_all[[1]],
          plot_all[[2]],
          nrow = 1)

dev.off()


# Supplementary Figure 3 - plot B

# Neurons which are mostly LEC and with at least 30 neurons
LEC_MEC_proportion_large_j <- lapply(1:J,
                                     function(j){
                                       
                                       if(length(which(unlist(mcmc_unique_EC8$Z) == j)) >= 30){
                                         LEC_proportion <- length(which(unlist(EC8_EC_label)[which(unlist(mcmc_unique_EC8$Z) == j)] == 'LEC'))/
                                           length(which(unlist(mcmc_unique_EC8$Z) == j))
                                         
                                         MEC_proportion <- 1-LEC_proportion
                                         
                                         data.frame(cluster = j,
                                                    EC = c('LEC','MEC'),
                                                    proportion = c(LEC_proportion, MEC_proportion))
                                       }
                                       
                                     })

LEC_MEC_proportion_large_j <- do.call(rbind, LEC_MEC_proportion_large_j)

large_LEC_dom_cluster <- LEC_MEC_proportion_large_j %>%
  group_by(cluster) %>%
  filter(proportion[1] > LEC_proportion_overall + 0.1) %>%
  ungroup() %>%
  select(cluster) %>%
  pull() %>%
  unique()

figure_3bs_LEC <- large_LEC_dom_cluster[which(large_LEC_dom_cluster %notin% figure_3b_cluster)]

# Reorder plots
figure_3bs_LEC_Nj <- sapply(figure_3bs_LEC,
                            function(j) length(which(unlist(mcmc_unique_EC8$Z)==j)))

figure_3bs_LEC <- figure_3bs_LEC[order(desc(figure_3bs_LEC_Nj))]


projecting.region.name <- sapply(figure_3bs_LEC, 
                                 function(j) rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[j,] > 0.5])

df_text <- data.frame(label  = sapply(1:length(figure_3bs_LEC),
                                      function(i){
                                        
                                        paste0(knitr::combine_words(projecting.region.name[[i]], and = ""))
                                      }),
                      
                      cluster = paste('cluster', figure_3bs_LEC))

df2 <- lapply(figure_3bs_LEC,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_med[j,],
                           lower_projection_strength = mcmc_unique_EC8$proj_prob_lower[j,],
                           upper_projection_strength = mcmc_unique_EC8$proj_prob_upper[j,])
              })

df2 <- do.call(rbind, df2)
df2$region <- factor(df2$region,
                     levels = rownames(EC8_new[[1]]))




# Proportion of LEC and MEC
figure_3bs_proportion <- lapply(figure_3bs_LEC,
                               function(j){
                                 
                                 LEC_proportion <-  length(which(unlist(EC8_EC_label)[which(unlist(mcmc_unique_EC8$Z) == j)] == 'LEC'))/
                                   length(which(unlist(mcmc_unique_EC8$Z) == j))
                                 MEC_proportion <- 1-LEC_proportion
                                 
                                 data.frame(cluster = j,
                                            EC = c('LEC','MEC'),
                                            proportion = c(LEC_proportion, MEC_proportion))
                               })

figure_3bs_proportion <- do.call(rbind, figure_3bs_proportion)


# For each cluster
plot1 <- lapply(figure_3bs_LEC,
                function(j){
                  
                  df2 %>%
                    filter(cluster == paste('cluster', j)) %>%
                    ggplot(mapping = aes(x = region,
                                         y = mean_projection_strength,
                                         group = 1))+
                    geom_line()+
                    geom_errorbar(mapping = aes(ymin = lower_projection_strength, 
                                                ymax = upper_projection_strength),
                                  width = 0.1)+
                    geom_point()+
                    theme_bw()+
                    ylab('projection probability')+
                    # geom_text(data = df_text[which(df_text$cluster == paste('cluster', j)),],
                    #           mapping = aes(x = 'VIS', y = 0.5, label = paste(strsplit(label, ', ')[[1]], collapse = '\n')),
                    #           hjust   = -0.1)+
                    theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          axis.text = element_text(size = 10,
                                                   color = 'black'))+
                    scale_y_continuous(breaks = c(0,1), limits = c(0,1))
                })


plot2 <- lapply(figure_3bs_LEC,
                function(j){
                  
                  
                  figure_3bs_proportion %>%
                    filter(cluster == j) %>%
                    ggplot(aes(x = factor(cluster, levels = figure_3bs_LEC),
                               y = proportion,
                               fill = EC))+
                    geom_bar(stat = 'identity')+
                    theme_bw()+
                    xlab('')+
                    ylab('')+
                    scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                    geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                               color = "red",
                               linewidth = 1.5)+
                    coord_flip()+
                    ggtitle(paste0((length(which(unlist(mcmc_unique_EC8$Z)==j))), ' neurons, ',
                                   length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'LEC']==j)), ' LEC, ',
                                   length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'MEC']==j)), ' MEC'))+
                    theme(legend.title=element_blank())+
                    theme(legend.position="none")+
                    xlab('')+
                    scale_y_continuous(breaks = c(0,1))+
                    theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          plot.title = element_text(size = 12))+
                    theme(axis.text = element_text(size = 12))
                  
                  
                })

plot_all <- lapply(1:length(figure_3bs_LEC), function(i) ggarrange(plot2[[i]], plot1[[i]], nrow = 2, heights = c(0.45,1)))


png(file = './plots/EC8_new2/figure_3bs.png',
    width = 900,
    height = 900/4)

ggarrange(plot_all[[1]],
          plot_all[[2]],
          plot_all[[3]],
          plot_all[[4]],
          nrow = 1)

dev.off()







# Supplementary Figure 3 - C
large_MEC_dom_cluster <- LEC_MEC_proportion_large_j %>%
  group_by(cluster) %>%
  filter(proportion[2] > MEC_proportion_overall + 0.1) %>%
  ungroup() %>%
  select(cluster) %>%
  pull() %>%
  unique()




figure_3cs_cluster <- large_MEC_dom_cluster[which(large_MEC_dom_cluster %notin% figure_3c_cluster)]


# Reorder plots
figure_3cs_Nj <- sapply(figure_3cs_cluster,
                        function(j) length(which(unlist(mcmc_unique_EC8$Z)==j)))

figure_3cs_cluster <- figure_3cs_cluster[order(desc(figure_3cs_Nj))]


projecting.region.name <- sapply(figure_3cs_cluster, 
                                 function(j) rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[j,] > 0.5])

df_text <- data.frame(label  = sapply(1:length(figure_3cs_cluster),
                                      function(i){
                                        
                                        paste0(knitr::combine_words(projecting.region.name[[i]], and = ""))
                                      }),
                      
                      cluster = paste('cluster', figure_3cs_cluster))

df2 <- lapply(figure_3cs_cluster,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_med[j,],
                           lower_projection_strength = mcmc_unique_EC8$proj_prob_lower[j,],
                           upper_projection_strength = mcmc_unique_EC8$proj_prob_upper[j,])
              })

df2 <- do.call(rbind, df2)
df2$region <- factor(df2$region,
                     levels = rownames(EC8_new[[1]]))





# Proportion of LEC and MEC
figure_3cs_proportion <- lapply(figure_3cs_cluster,
                                function(j){
                                  
                                  LEC_proportion <- length(which(unlist(EC8_EC_label)[which(unlist(mcmc_unique_EC8$Z) == j)] == 'LEC'))/
                                    length(which(unlist(mcmc_unique_EC8$Z) == j))
                                  MEC_proportion <- 1-LEC_proportion
                                  
                                  data.frame(cluster = j,
                                             EC = c('LEC','MEC'),
                                             proportion = c(LEC_proportion, MEC_proportion))
                                })

figure_3cs_proportion <- do.call(rbind, figure_3cs_proportion)


# For each cluster
plot1 <- lapply(figure_3cs_cluster,
                function(j){
                  
                  df2 %>%
                    filter(cluster == paste('cluster', j)) %>%
                    ggplot(mapping = aes(x = region,
                                         y = mean_projection_strength,
                                         group = 1))+
                    geom_line()+
                    geom_errorbar(mapping = aes(ymin = lower_projection_strength, 
                                                ymax = upper_projection_strength),
                                  width = 0.1)+
                    geom_point()+
                    theme_bw()+
                    ylab('projection probability')+
                    # geom_text(data = df_text[which(df_text$cluster == paste('cluster', j)),],
                    #           mapping = aes(x = 'VIS', y = 0.5, label = paste(strsplit(label, ', ')[[1]], collapse = '\n')),
                    #           hjust   = -0.1)+
                    theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          axis.text = element_text(size = 10,
                                                   color = 'black'))+
                    scale_y_continuous(breaks = c(0,1), limits = c(0,1))
                })

plot2 <- lapply(figure_3cs_cluster,
                function(j){
                  
                  figure_3cs_proportion %>%
                    filter(cluster == j) %>%
                    ggplot(aes(x = factor(cluster, levels = figure_3cs_cluster),
                               y = proportion,
                               fill = EC))+
                    geom_bar(stat = 'identity')+
                    theme_bw()+
                    xlab('')+
                    ylab('')+
                    scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                    geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                               color = "red",
                               linewidth = 1.5)+
                    coord_flip()+
                    ggtitle(paste0((length(which(unlist(mcmc_unique_EC8$Z)==j))), ' neurons, ',
                                   length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'LEC']==j)), ' LEC, ',
                                   length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'MEC']==j)), ' MEC'))+
                    theme(legend.title=element_blank())+
                    theme(legend.position="none")+
                    xlab('')+
                    scale_y_continuous(breaks = c(0,1))+
                    theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          plot.title = element_text(size = 12))+
                    theme(axis.text = element_text(size = 12))
                  
                  
                })

plot_all <- lapply(1:length(figure_3cs_cluster), function(i) ggarrange(plot2[[i]], plot1[[i]], nrow = 2, heights = c(0.45,1)))


png(file = './plots/EC8_new2/figure_3cd_supp_plot1.png',
    width = 1125,
    height = 1125/5)

ggarrange(plot_all[[1]],
          plot_all[[2]],
          plot_all[[3]],
          plot_all[[4]],
          plot_all[[5]],
          nrow = 1)

dev.off()


png(file = './plots/EC8_new2/figure_3cd_supp_plot2.png',
    width = 1125,
    height = 1125/5)

ggarrange(plot_all[[6]],
          plot_all[[7]],
          plot_all[[8]],
          plot_all[[9]],
          plot_all[[10]],
          nrow = 1)

dev.off()


# Figure 3D - supplementary
figure_3ds_cluster <- 34

projecting.region.name <- sapply(figure_3ds_cluster, 
                                 function(j) rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[j,] > 0.5])

df_text <- data.frame(label  = sapply(1:length(figure_3ds_cluster),
                                      function(i){
                                        
                                        paste0(knitr::combine_words(projecting.region.name[[i]], and = ""))
                                      }),
                      
                      cluster = paste('cluster', figure_3ds_cluster))

df2 <- lapply(figure_3ds_cluster,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_med[j,],
                           lower_projection_strength = mcmc_unique_EC8$proj_prob_lower[j,],
                           upper_projection_strength = mcmc_unique_EC8$proj_prob_upper[j,])
              })

df2 <- do.call(rbind, df2)
df2$region <- factor(df2$region,
                     levels = rownames(EC8_new[[1]]))



# Proportion of LEC and MEC
figure_3ds_proportion <- lapply(figure_3ds_cluster,
                                function(j){
                                  
                                  LEC_proportion <- length(which(unlist(EC8_EC_label)[which(unlist(mcmc_unique_EC8$Z) == j)] == 'LEC'))/
                                    length(which(unlist(mcmc_unique_EC8$Z) == j))
                                  MEC_proportion <- 1-LEC_proportion
                                  
                                  data.frame(cluster = j,
                                             EC = c('LEC','MEC'),
                                             proportion = c(LEC_proportion, MEC_proportion))
                                })

figure_3ds_proportion <- do.call(rbind, figure_3ds_proportion)


plot1 <- 
  df2 %>%
  filter(cluster == paste('cluster', 34)) %>%
  ggplot(mapping = aes(x = region,
                       y = mean_projection_strength,
                       group = 1))+
  geom_line()+
  geom_errorbar(mapping = aes(ymin = lower_projection_strength, 
                              ymax = upper_projection_strength),
                width = 0.1)+
  geom_point()+
  theme_bw()+
  ylab('projection probability')+
  # geom_text(data = df_text[which(df_text$cluster == paste('cluster', 34)),],
  #           mapping = aes(x = 'VIS', y = 0.5, label = paste(strsplit(label, ', ')[[1]], collapse = '\n')),
  #           hjust   = -0.1)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 10,
                                 color = 'black'))+
  scale_y_continuous(breaks = c(0,1), limits = c(0,1))

plot2 <- 
  figure_3ds_proportion %>%
  filter(cluster == 34) %>%
  ggplot(aes(x = factor(cluster, levels = figure_3ds_cluster),
             y = proportion,
             fill = EC))+
  geom_bar(stat = 'identity')+
  theme_bw()+
  xlab('')+
  ylab('')+
  scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
  geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
             color = "red",
             linewidth = 1.5)+
  coord_flip()+
  ggtitle(paste0((length(which(unlist(mcmc_unique_EC8$Z)==j))), ' neurons, ',
                 length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'LEC']==j)), ' LEC, ',
                 length(which(unlist(mcmc_unique_EC8$Z)[unlist(EC8_EC_label) == 'MEC']==j)), ' MEC'))+
  theme(legend.title=element_blank())+
  theme(legend.position="none")+
  xlab('')+
  scale_y_continuous(breaks = c(0,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 12))+
  theme(axis.text = element_text(size = 12))




png(file = './plots/EC8_new2/figure_3d_supp.png',
    width = 225,
    height = 225)

ggarrange(plot2, plot1, nrow = 2, heights = c(0.45,1))

dev.off()

















# Finding the projecting regions of each cluster
projecting.region.name.all <- sapply(1:J, 
                                 function(j){
                                   
                                   name <- rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[j,] > 0.5]
                                   paste0('[', knitr::combine_words(name, and = ""), ']')
                                 })

table(projecting.region.name.all)


#Initialize z
data_EC_cbind <- do.call(cbind, EC8_new)

df <- t(data_EC_cbind)
df = apply(df, 1, function(x){return(x/sum(x))})


# Projection to VIS and RSC
VIS_RSC_cluster <- which(projecting.region.name.all == '[VIS, RSC]')

EC8_VIS_RSC <- df[,which(unlist(mcmc_unique_EC8$Z) %in% VIS_RSC_cluster)]
EC8_allocation_VIS_RSC <- unlist(mcmc_unique_EC8$Z)[which(unlist(mcmc_unique_EC8$Z) %in% VIS_RSC_cluster)]

df_VIS_RSC <- data.frame(projection_strength = as.vector(EC8_VIS_RSC),
                         cluster = rep(EC8_allocation_VIS_RSC, each = 8),
                         mouse_index = rep(1:length(EC8_allocation_VIS_RSC), each = 8),
                         brain_region = rownames(EC8_new[[1]]))

df_VIS_RSC$cluster <- factor(df_VIS_RSC$cluster,
                             levels = VIS_RSC_cluster)
df_VIS_RSC$brain_region <- factor(df_VIS_RSC$brain_region,
                                  levels = rownames(EC8_new[[1]]))
df_VIS_RSC$mouse_index <- factor(df_VIS_RSC$mouse_index,
                                 levels = 1:max(df_VIS_RSC$mouse_index))

# Proportion of LEC and MEC in each cluster
VIS_RSC_EC_prop <- rbind(data.frame(cluster = 25,
                                    EC = c('LEC', 'MEC'),
                                    count = c(length(which(unlist(EC8_EC_label)[unlist(mcmc_unique_EC8$Z)==25]=='LEC')),
                                              length(which(unlist(EC8_EC_label)[unlist(mcmc_unique_EC8$Z)==25]=='MEC')))),
                         
                         data.frame(cluster = 37,
                                    EC = c('LEC', 'MEC'),
                                    count = c(length(which(unlist(EC8_EC_label)[unlist(mcmc_unique_EC8$Z)==37]=='LEC')),
                                              length(which(unlist(EC8_EC_label)[unlist(mcmc_unique_EC8$Z)==37]=='MEC')))))

VIS_RSC_EC_prop$cluster <- factor(VIS_RSC_EC_prop$cluster,
                                  levels = c(25, 37))

VIS_RSC_EC_prop <- VIS_RSC_EC_prop %>%
  group_by(cluster) %>%
  mutate(proportion = c(count[1]/sum(count),
                        count[2]/sum(count))) %>%
  ungroup()



plot11 <- VIS_RSC_EC_prop %>%
  ggplot(aes(x = cluster,
             y = proportion,
             fill = EC))+
  geom_bar(stat = 'identity')+
  theme_bw()+
  xlab('cluster')+
  ylab('neuron proportion')+
  scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
  geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
             color = "red",
             linewidth = 1.5)+
  coord_flip()+
  xlab('')+
  scale_y_continuous(breaks = c(0,1))+
  ylab('')+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 10))+
  theme(legend.position="none")

plot12 <- df_VIS_RSC %>%
  ggplot()+
  geom_line(mapping = aes(x = brain_region,
                          y = projection_strength,
                          color = cluster,
                          group = mouse_index))+
  theme_bw()+
  xlab('region')+
  ylab('projection probability')+
  ggtitle('[VIS, RSC]')+
  scale_y_continuous(breaks = c(0,1), limits = c(0,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 10,
                                 color = 'black'),
        legend.position = c(0.2, 0.8))


plot1 <- ggarrange(plot11, plot12, nrow = 2, heights = c(0.3,1))




# Projection to SS and VIS
SS_VIS_cluster <- which(projecting.region.name.all == '[SS, VIS]')

EC8_SS_VIS <- df[,which(unlist(mcmc_unique_EC8$Z) %in% SS_VIS_cluster)]
EC8_allocation_SS_VIS <- unlist(mcmc_unique_EC8$Z)[which(unlist(mcmc_unique_EC8$Z) %in% SS_VIS_cluster)]


df_SS_VIS <- data.frame(projection_strength = as.vector(EC8_SS_VIS),
                         cluster = rep(EC8_allocation_SS_VIS, each = 8),
                         mouse_index = rep(1:length(EC8_allocation_SS_VIS), each = 8),
                         brain_region = rownames(EC8_new[[1]]))

df_SS_VIS$cluster <- factor(df_SS_VIS$cluster,
                             levels = SS_VIS_cluster)
df_SS_VIS$brain_region <- factor(df_SS_VIS$brain_region,
                                  levels = rownames(EC8_new[[1]]))
df_SS_VIS$mouse_index <- factor(df_SS_VIS$mouse_index,
                                 levels = 1:max(df_SS_VIS$mouse_index))


# Proportion of LEC and MEC in each cluster
SS_VIS_EC_prop <- rbind(data.frame(cluster = 16,
                                    EC = c('LEC', 'MEC'),
                                    count = c(length(which(unlist(EC8_EC_label)[unlist(mcmc_unique_EC8$Z)==16]=='LEC')),
                                              length(which(unlist(EC8_EC_label)[unlist(mcmc_unique_EC8$Z)==16]=='MEC')))),
                         
                         data.frame(cluster = 34,
                                    EC = c('LEC', 'MEC'),
                                    count = c(length(which(unlist(EC8_EC_label)[unlist(mcmc_unique_EC8$Z)==34]=='LEC')),
                                              length(which(unlist(EC8_EC_label)[unlist(mcmc_unique_EC8$Z)==34]=='MEC')))))

SS_VIS_EC_prop$cluster <- factor(SS_VIS_EC_prop$cluster,
                                  levels = c(16, 34))

SS_VIS_EC_prop <- SS_VIS_EC_prop %>%
  group_by(cluster) %>%
  mutate(proportion = c(count[1]/sum(count),
                        count[2]/sum(count))) %>%
  ungroup()



plot21 <- SS_VIS_EC_prop %>%
  ggplot(aes(x = cluster,
             y = proportion,
             fill = EC))+
  geom_bar(stat = 'identity')+
  theme_bw()+
  xlab('cluster')+
  ylab('neuron proportion')+
  scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
  geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
             color = "red",
             linewidth = 1.5)+
  coord_flip()+
  xlab('')+
  scale_y_continuous(breaks = c(0,1))+
  ylab('')+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 10))+
  theme(legend.position="none")


plot22 <- df_SS_VIS %>%
  ggplot()+
  geom_line(mapping = aes(x = brain_region,
                          y = projection_strength,
                          color = cluster,
                          group = mouse_index))+
  theme_bw()+
  xlab('region')+
  ylab('projection probability')+
  ggtitle('[SS, VIS]')+
  scale_y_continuous(breaks = c(0,1), limits = c(0,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 10,
                                 color = 'black'),
        legend.position = c(0.2, 0.8))


plot2 <- ggarrange(plot21, plot22, nrow = 2, heights = c(0.3,1))

ggarrange(plot1, plot2)

