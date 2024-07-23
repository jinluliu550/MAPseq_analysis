


binomial_model2 <- function(data){
  
  
  # If the data is supplied in a list
  if(is.list(data)){
    
    M <- length(data)
    C <- sapply(1:length(data), function(m) ncol(data[[m]]))
    data <- do.call(cbind, data)
    
  }else{
    
    M <- 1
    C <- ncol(data)
  }
  
  R <- nrow(data)
  
  
  
  # Calculate Nf, N_0 and p_ea
  Nf <- ncol(data)
  
  si <- sapply(1:nrow(data),
               function(r) length(which(data[r,] >= 5)))
  
  N_0 <- round(fzero(function(x) Nf - x*(1-prod(1-si/x)), Nf)$x, 0)
  
  
  p_i <- si/N_0
  
  
  
  # Find the set of unique clusters
  cluster_neuron <- sapply(1:ncol(data),
                           function(c) paste(rownames(data)[data[,c] >= 5],
                                             collapse = " "))
  
  # Unique set of cluster
  cluster_label <- unique(cluster_neuron)
  
  # Region names
  region_names <- rownames(data)
  
  # Neuron allocation
  allocation <- sapply(1:ncol(data),
                       function(c) which(cluster_label == cluster_neuron[c]))
  
  
  # Separate allocations
  C_cumsum <- c(0, cumsum(C))
  
  allocation <- lapply(1:M,
                       function(m) allocation[(C_cumsum[m]+1):(C_cumsum[m+1])])
  
  
  cluster_summary <- data.frame(cluster_name = cluster_label,
                                N_projecting_regions = sapply(1:length(cluster_label),
                                                              function(i) length(str_split(cluster_label[[i]], pattern = ' ')[[1]]))
  )
  
  
  
  # Expected number of projections
  cluster_summary$size_expect <- sapply(1:nrow(cluster_summary),
                           function(j){
                             
                             prod(c(p_i[which(rownames(data) %in%  str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])],
                                    1-p_i[which(rownames(data) %notin%  str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])]))*N_0
                           }
  )
  
  
  # Observed number of projections
  cluster_summary$size_observed <- sapply(1:nrow(cluster_summary),
                                          function(j) length(which(cluster_neuron == cluster_label[j])))
  
  cluster_summary <- cluster_summary %>%
    mutate(cluster.type = ifelse(size_expect < size_observed, 'over-represented', 'under-represented'))
  
  
  # p-value
  cluster_summary$p_value <- sapply(1:nrow(cluster_summary),
                                    function(j) binom.test(x = cluster_summary$size_observed[j], 
                                                           n = N_0, 
                                                           p = prod(c(p_i[which(rownames(data) %in%  str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])],
                                                                      1-p_i[which(rownames(data) %notin%  str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])])), 
                                                           alternative = 'two.sided')$p.value)
  
  
  threshold_p_value <- -log(0.05/(M*length(which(cluster_summary$N_projecting_regions > 1))), base = 10)
  
  plot_output <- cluster_summary %>%
    filter(N_projecting_regions > 1) %>%
    ggplot(mapping = aes(x = log(size_observed/size_expect, base = 2),
                         y = -log(p_value, base = 10)))+
    geom_point()+
    geom_hline(yintercept = threshold_p_value, linetype = 'dashed', color = 'grey')+
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey')+
    theme_bw()
  
  
  # Significance
  cluster_summary$significant <- ifelse(-log(cluster_summary$p_value, base = 10) < threshold_p_value,
                                        'insignificant',
                                        'significant')
  
  cluster_summary %>%
    filter(significant == 'significant')
  
  return(list('allocation' = allocation,
              'cluster_label' = cluster_label,
              'cluster_summary' = cluster_summary,
              'plot_output' = plot_output))
}


