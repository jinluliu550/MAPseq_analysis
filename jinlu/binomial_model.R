
binomial_model <- function(data){
  
  M <- length(data)
  
  C <- sapply(1:length(data),
              function(m) ncol(data[[m]]))
  
  data <- do.call(cbind, data)
  
  # Find the set of unique clusters
  cluster_neuron <- sapply(1:ncol(data),
                           function(c) paste(rownames(data)[data[,c] > 0],
                                             collapse = " "))
  
  # Unique set of cluster
  cluster_label <- unique(cluster_neuron)
  
  # Region names
  region_names <- rownames(data)
  
  # Neuron allocation
  allocation <- sapply(1:ncol(data),
                       function(c) which(cluster_label == cluster_neuron[c]))
  
  # Cluster summary
  projecting_region <- lapply(1:length(cluster_label),
                              function(j){
                                
                                decision <- ifelse(region_names %in% str_split(cluster_label[j], " ")[[1]], 1, 0)
                              })
  
  # Separate allocations
  C_cumsum <- c(0, cumsum(C))
  
  allocation <- lapply(1:M,
                       function(m) allocation[(C_cumsum[m]+1):(C_cumsum[m+1])])
  
  return(list('allocation' = allocation,
              'cluster_label' = cluster_label))
              # 'cluster_summary' = cluster_summary_df))
}

