
#' Plot to show the number of neurons in each cluster, and the contribution of each mouse in each cluster
#'

opt.clustering.frequency <- function(clustering,
                                     main = ''){

  cluster.names <- paste('cluster', 1:length(unique(unlist(clustering))))

  loop.result <- lapply(1:length(clustering), function(m){

    data.frame(cluster = factor(cluster.names,
                                levels = cluster.names),
               mouse = paste('mouse', m),
               counts = as.vector(table(factor(clustering[[m]],
                                               levels = 1:length(cluster.names)))))
  })

  # Frequency table
  z.frequency <- do.call(rbind, loop.result)

  # Bar chart
  z.frequency %>%
    ggplot(mapping = aes(x = cluster, y = counts, fill = mouse))+
    geom_bar(stat="identity")+
    theme_bw()+
    ylab('Number of neurons')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(main)
}
