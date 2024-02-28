# In code below, we summarize the projecting location of each cluster using a heat-map

cluster_by_projecting_location <- function(mcmc_unique_mcmc_output){
  
  df <- mcmc_unique_mcmc_output$cluster.label.summary
  
  
  df1 <- data.frame(cluster = paste('cluster', rep(df$cluster, each = 7)),
                    projecting_location = colnames(df)[5:11],
                    fill = unlist(lapply(1:nrow(df),
                                         function(i) ifelse(df[i,5:11] == 'Yes', 1, 0))))
  
  
  # Convert to factors
  df1$cluster <- factor(df1$cluster,
                        levels = paste('cluster', 1:nrow(df)))
  
  df1$projecting_location <- factor(df1$projecting_location,
                                    levels = c('project.to.front',
                                               'project.to.front.middle',
                                               'project.to.middle',
                                               'project.to.back.middle',
                                               'project.to.back',
                                               'project.to.front.back',
                                               'project.to.all.area'))
  
  df1 %>%
    ggplot(mapping = aes(x = cluster,
                         y = projecting_location,
                         fill = fill))+
    geom_tile()+
    theme_bw()+
    scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                         values = scales::rescale(c(0, 0.25, 1)),
                         limits = c(0,1))+
    xlab('cluster')+
    ylab('regions')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")
  
}