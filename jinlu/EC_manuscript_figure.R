
# Figure 3b


# Cluster 1, 6, 14, 5
figure_3b_cluster <- c(1,6,14,5)

projecting.region.name <- sapply(1:4, 
                                 function(i) rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[figure_3b_cluster[i],] > 0.5])

df_text <- data.frame(label  = sapply(1:4,
                                      function(i){
                                        
                                        paste0('[', knitr::combine_words(projecting.region.name[[i]], and = ""), ']')
                                      }),
                      
                      cluster = paste('cluster', figure_3b_cluster))

df2 <- lapply(figure_3b_cluster,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_mean[j,])
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
                    geom_point()+
                    theme_bw()+
                    ylab('projection strength')+
                    geom_text(data = df_text[which(df_text$cluster == paste('cluster', j)),],
                              mapping = aes(x = -Inf, y = 0.9, label = label),
                              hjust   = -0.1)+
                    ylim(c(0,1))
                })

plot2 <- lapply(figure_3b_cluster,
                function(j){
                  
                  if(j == 5){
                    
                    figure_3b_proportion %>%
                      filter(cluster == j) %>%
                      ggplot(aes(x = factor(cluster, levels = figure_3b_cluster),
                                 y = proportion,
                                 fill = EC))+
                      geom_bar(stat = 'identity')+
                      theme_bw()+
                      xlab('cluster')+
                      ylab('neuron proportion')+
                      scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                      geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                                 color = "red")+
                      coord_flip()+
                      ggtitle(paste0('cluster ', j, ' (', length(which(unlist(mcmc_unique_EC8$Z)==j)), ' neurons)'))+
                      theme(legend.title=element_blank())+
                      xlab('')+
                      scale_y_continuous(breaks = c(0,1))
                    
                  }else{
                    
                    figure_3b_proportion %>%
                      filter(cluster == j) %>%
                      ggplot(aes(x = factor(cluster, levels = figure_3b_cluster),
                                 y = proportion,
                                 fill = EC))+
                      geom_bar(stat = 'identity')+
                      theme_bw()+
                      xlab('cluster')+
                      ylab('neuron proportion')+
                      scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                      geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                                 color = "red")+
                      coord_flip()+
                      ggtitle(paste0('cluster ', j, ' (', length(which(unlist(mcmc_unique_EC8$Z)==j)), ' neurons)'))+
                      theme(legend.title=element_blank())+
                      theme(legend.position="none")+
                      xlab('')+
                      scale_y_continuous(breaks = c(0,1))
                  }
                  
                })

plot_all <- lapply(1:4, function(i) ggarrange(plot2[[i]], plot1[[i]], nrow = 2, heights = c(0.3,1)))


png(file = './plots/EC8_new2/figure_3b.png',
    width = 1200,
    height = 300)

ggarrange(plot_all[[1]],
          plot_all[[2]],
          plot_all[[3]],
          plot_all[[4]],
          nrow = 1)

dev.off()





# Figure 3c
figure_3c_cluster <- c(25,41,37,55)



projecting.region.name <- sapply(1:4, 
                                 function(i) rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[figure_3c_cluster[i],] > 0.5])

df_text <- data.frame(label  = sapply(1:4,
                                      function(i){
                                        
                                        paste0('[', knitr::combine_words(projecting.region.name[[i]], and = ""), ']')
                                      }),
                      
                      cluster = paste('cluster', figure_3c_cluster))

df2 <- lapply(figure_3c_cluster,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_mean[j,])
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
                    geom_point()+
                    theme_bw()+
                    ylab('projection strength')+
                    geom_text(data = df_text[which(df_text$cluster == paste('cluster', j)),],
                              mapping = aes(x = -Inf, y = 0.9, label = label),
                              hjust   = -0.1)+
                    ylim(c(0,1))
                })

plot2 <- lapply(figure_3c_cluster,
                function(j){
                  
                  if(j == 55){
                    
                    figure_3c_proportion %>%
                      filter(cluster == j) %>%
                      ggplot(aes(x = factor(cluster, levels = figure_3c_cluster),
                                 y = proportion,
                                 fill = EC))+
                      geom_bar(stat = 'identity')+
                      theme_bw()+
                      xlab('cluster')+
                      ylab('neuron proportion')+
                      scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                      geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                                 color = "red")+
                      coord_flip()+
                      ggtitle(paste0('cluster ', j, ' (', length(which(unlist(mcmc_unique_EC8$Z)==j)), ' neurons)'))+
                      theme(legend.title=element_blank())+
                      xlab('')+
                      scale_y_continuous(breaks = c(0,1))
                    
                  }else{
                    
                    figure_3c_proportion %>%
                      filter(cluster == j) %>%
                      ggplot(aes(x = factor(cluster, levels = figure_3c_cluster),
                                 y = proportion,
                                 fill = EC))+
                      geom_bar(stat = 'identity')+
                      theme_bw()+
                      xlab('cluster')+
                      ylab('neuron proportion')+
                      scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                      geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                                 color = "red")+
                      coord_flip()+
                      ggtitle(paste0('cluster ', j, ' (', length(which(unlist(mcmc_unique_EC8$Z)==j)), ' neurons)'))+
                      theme(legend.title=element_blank())+
                      theme(legend.position="none")+
                      xlab('')+
                      scale_y_continuous(breaks = c(0,1))
                  }
                  
                })

plot_all <- lapply(1:4, function(i) ggarrange(plot2[[i]], plot1[[i]], nrow = 2, heights = c(0.3,1)))


png(file = './plots/EC8_new2/figure_3c.png',
    width = 1200,
    height = 300)

ggarrange(plot_all[[1]],
          plot_all[[2]],
          plot_all[[3]],
          plot_all[[4]],
          nrow = 1)

dev.off()


# Figure 3d
figure_3d_cluster <- c(48,54,34)



projecting.region.name <- sapply(1:3, 
                                 function(i) rownames(EC8_new[[1]])[mcmc_unique_EC8$q_tilde_001[figure_3d_cluster[i],] > 0.5])

df_text <- data.frame(label  = sapply(1:3,
                                      function(i){
                                        
                                        paste0('[', knitr::combine_words(projecting.region.name[[i]], and = ""), ']')
                                      }),
                      
                      cluster = paste('cluster', figure_3d_cluster))

df2 <- lapply(figure_3d_cluster,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_mean[j,])
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
                    geom_point()+
                    theme_bw()+
                    ylab('projection strength')+
                    geom_text(data = df_text[which(df_text$cluster == paste('cluster', j)),],
                              mapping = aes(x = -Inf, y = 0.9, label = label),
                              hjust   = -0.1)+
                    ylim(c(0,1))
                })



plot2 <- lapply(figure_3d_cluster,
                function(j){
                  
                  
                  if(j == 34){
                    
                    figure_3d_proportion %>%
                      filter(cluster == j) %>%
                      ggplot(aes(x = factor(cluster, levels = figure_3d_cluster),
                                 y = proportion,
                                 fill = EC))+
                      geom_bar(stat = 'identity')+
                      theme_bw()+
                      xlab('cluster')+
                      ylab('neuron proportion')+
                      scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                      geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                                 color = "red")+
                      coord_flip()+
                      ggtitle(paste0('cluster ', j, ' (', length(which(unlist(mcmc_unique_EC8$Z)==j)), ' neurons)'))+
                      theme(legend.title=element_blank())+
                      xlab('')+
                      scale_y_continuous(breaks = c(0,1))
                    
                  }else{
                    
                    figure_3d_proportion %>%
                      filter(cluster == j) %>%
                      ggplot(aes(x = factor(cluster, levels = figure_3d_cluster),
                                 y = proportion,
                                 fill = EC))+
                      geom_bar(stat = 'identity')+
                      theme_bw()+
                      xlab('cluster')+
                      ylab('neuron proportion')+
                      scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                      geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                                 color = "red")+
                      coord_flip()+
                      ggtitle(paste0('cluster ', j, ' (', length(which(unlist(mcmc_unique_EC8$Z)==j)), ' neurons)'))+
                      theme(legend.title=element_blank())+
                      theme(legend.position="none")+
                      xlab('')+
                      scale_y_continuous(breaks = c(0,1))
                  }
                  
                })

plot_all <- lapply(1:3, function(i) ggarrange(plot2[[i]], plot1[[i]], nrow = 2, heights = c(0.3,1)))


png(file = './plots/EC8_new2/figure_3d.png',
    width = 900,
    height = 300)

ggarrange(plot_all[[1]],
          plot_all[[2]],
          plot_all[[3]],
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
                                        
                                        paste0('[', knitr::combine_words(projecting.region.name[[i]], and = ""), ']')
                                      }),
                      
                      cluster = paste('cluster', figure_3bs_LEC))

df2 <- lapply(figure_3bs_LEC,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_mean[j,])
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
                    geom_point()+
                    theme_bw()+
                    ylab('projection strength')+
                    geom_text(data = df_text[which(df_text$cluster == paste('cluster', j)),],
                              mapping = aes(x = -Inf, y = 0.9, label = label),
                              hjust   = -0.1)+
                    ylim(c(0,1))
                })


plot2 <- lapply(figure_3bs_LEC,
                function(j){
                  
                  
                  if(j == figure_3bs_LEC[length(figure_3bs_LEC)]){
                    
                    figure_3bs_proportion %>%
                      filter(cluster == j) %>%
                      ggplot(aes(x = factor(cluster, levels = figure_3bs_LEC),
                                 y = proportion,
                                 fill = EC))+
                      geom_bar(stat = 'identity')+
                      theme_bw()+
                      xlab('cluster')+
                      ylab('neuron proportion')+
                      scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                      geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                                 color = "red")+
                      coord_flip()+
                      ggtitle(paste0('cluster ', j, ' (', length(which(unlist(mcmc_unique_EC8$Z)==j)), ' neurons)'))+
                      theme(legend.title=element_blank())+
                      xlab('')+
                      scale_y_continuous(breaks = c(0,1))
                    
                  }else{
                    
                    figure_3bs_proportion %>%
                      filter(cluster == j) %>%
                      ggplot(aes(x = factor(cluster, levels = figure_3bs_LEC),
                                 y = proportion,
                                 fill = EC))+
                      geom_bar(stat = 'identity')+
                      theme_bw()+
                      xlab('cluster')+
                      ylab('neuron proportion')+
                      scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                      geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                                 color = "red")+
                      coord_flip()+
                      ggtitle(paste0('cluster ', j, ' (', length(which(unlist(mcmc_unique_EC8$Z)==j)), ' neurons)'))+
                      theme(legend.title=element_blank())+
                      theme(legend.position="none")+
                      xlab('')+
                      scale_y_continuous(breaks = c(0,1))
                  }
                })

plot_all <- lapply(1:length(figure_3bs_LEC), function(i) ggarrange(plot2[[i]], plot1[[i]], nrow = 2, heights = c(0.3,1)))


png(file = './plots/EC8_new2/figure_3bs.png',
    width = 1200,
    height = 300)

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
                                        
                                        paste0('[', knitr::combine_words(projecting.region.name[[i]], and = ""), ']')
                                      }),
                      
                      cluster = paste('cluster', figure_3cs_cluster))

df2 <- lapply(figure_3cs_cluster,
              function(j){
                
                data.frame(cluster = paste('cluster', j),
                           region = rownames(EC8_new[[1]]),
                           mean_projection_strength = mcmc_unique_EC8$proj_prob_mean[j,])
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
                    geom_point()+
                    theme_bw()+
                    ylab('projection strength')+
                    geom_text(data = df_text[which(df_text$cluster == paste('cluster', j)),],
                              mapping = aes(x = -Inf, y = 0.9, label = label),
                              hjust   = -0.1)+
                    ylim(c(0,1))
                })

plot2 <- lapply(figure_3cs_cluster,
                function(j){
                  
                  
                  if(j %in% c(64, 30)){
                    
                    figure_3cs_proportion %>%
                      filter(cluster == j) %>%
                      ggplot(aes(x = factor(cluster, levels = figure_3cs_cluster),
                                 y = proportion,
                                 fill = EC))+
                      geom_bar(stat = 'identity')+
                      theme_bw()+
                      xlab('cluster')+
                      ylab('neuron proportion')+
                      scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                      geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                                 color = "red")+
                      coord_flip()+
                      ggtitle(paste0('cluster ', j, ' (', length(which(unlist(mcmc_unique_EC8$Z)==j)), ' neurons)'))+
                      theme(legend.title=element_blank())+
                      xlab('')+
                      scale_y_continuous(breaks = c(0,1))
                    
                  }else{
                    
                    figure_3cs_proportion %>%
                      filter(cluster == j) %>%
                      ggplot(aes(x = factor(cluster, levels = figure_3cs_cluster),
                                 y = proportion,
                                 fill = EC))+
                      geom_bar(stat = 'identity')+
                      theme_bw()+
                      xlab('cluster')+
                      ylab('neuron proportion')+
                      scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
                      geom_hline(yintercept=MEC_proportion_overall, linetype="solid",
                                 color = "red")+
                      coord_flip()+
                      ggtitle(paste0('cluster ', j, ' (', length(which(unlist(mcmc_unique_EC8$Z)==j)), ' neurons)'))+
                      theme(legend.title=element_blank())+
                      theme(legend.position="none")+
                      xlab('')+
                      scale_y_continuous(breaks = c(0,1))
                  }
                })

plot_all <- lapply(1:length(figure_3cs_cluster), function(i) ggarrange(plot2[[i]], plot1[[i]], nrow = 2, heights = c(0.3,1)))


png(file = './plots/EC8_new2/figure_3cd_supp_plot1.png',
    width = 1500,
    height = 300)

ggarrange(plot_all[[1]],
          plot_all[[2]],
          plot_all[[3]],
          plot_all[[4]],
          plot_all[[5]],
          nrow = 1)

dev.off()


png(file = './plots/EC8_new2/figure_3cd_supp_plot2.png',
    width = 1500,
    height = 300)

ggarrange(plot_all[[6]],
          plot_all[[7]],
          plot_all[[8]],
          plot_all[[9]],
          plot_all[[10]],
          nrow = 1)

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


png(file = './plots/EC8_new2/VIS_RSC.png',
    width = 400,
    height = 300)

df_VIS_RSC %>%
  ggplot()+
  geom_line(mapping = aes(x = brain_region,
                          y = projection_strength,
                          color = cluster,
                          group = mouse_index))+
  theme_bw()+
  xlab('region')+
  ylab('projection strength')+
  ggtitle('[VIS, RSC]')

dev.off()


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


png(file = './plots/EC8_new2/SS_VIS.png',
    width = 400,
    height = 300)

df_SS_VIS %>%
  ggplot()+
  geom_line(mapping = aes(x = brain_region,
                          y = projection_strength,
                          color = cluster,
                          group = mouse_index))+
  theme_bw()+
  xlab('region')+
  ylab('projection strength')+
  ggtitle('[SS, VIS]')

dev.off()

