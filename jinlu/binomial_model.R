


binomial_model <- function(data){
  
  
  # If the data is supplied in a list
  if(is.list(data)){
    
    C <- sapply(1:length(data), function(m) ncol(data[[m]]))
    data <- do.call(cbind, data)
    
  }else{
    
    C <- ncol(data)
  }
  
  R <- nrow(data)
  
  
  
  # Calculate Nf, N_0 and p_e
  Nf <- ncol(data)
  
  si <- sapply(1:nrow(data),
               function(r) length(which(data[r,] != 0)))
  
  N_0 <- round(fzero(function(x) Nf - x*(1-prod(1-si/x)), Nf)$x, 0)
  
  p_e <- fzero(function(x) Nf - N_0*(1-(1-x)^nrow(data)), 0.5)$x
  
  # Compute the expected count for different number of projecting motifs
  expected_number <- sapply(1:nrow(data),
                            function(r){
                              
                              N_0*p_e^(r)*(1-p_e)^(R-r)
                            })
  
  
  
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
  
  
  # Separate allocations
  C_cumsum <- c(0, cumsum(C))
  
  allocation <- lapply(1:M,
                       function(m) allocation[(C_cumsum[m]+1):(C_cumsum[m+1])])
  
  
  cluster_summary <- data.frame(cluster_name = cluster_label,
                                N_projecting_regions = sapply(1:length(cluster_label),
                                                                      function(i) length(str_split(cluster_label[[i]], pattern = ' ')[[1]]))
                                )
  
  # Expected number of projections
  cluster_summary$size_expect <- expected_number[cluster_summary$N_projecting_regions]
  
  # Observed number of projections
  cluster_summary$size_observed <- sapply(1:nrow(cluster_summary),
                                          function(j) length(which(cluster_neuron == cluster_label[j])))
  
  cluster_summary <- cluster_summary %>%
    mutate(cluster.type = ifelse(size_expect < size_observed, 'over-represented', 'under-represented'))
  
  # p-value
  cluster_summary$p_value <- sapply(1:nrow(cluster_summary),
                                    function(j) binom.test(x = cluster_summary$size_observed[j], 
                                                           n = ncol(data), 
                                                           p = (p_e)^cluster_summary$N_projecting_regions[j]*(1-p_e)^(R-cluster_summary$N_projecting_regions[j]), 
                                                           alternative = 'two.sided')$p.value)
  
  
  # Threshold
  threshold_p_value <- -log(0.05/nrow(cluster_summary), base = 10)
  
  cluster_summary$significant <- ifelse(-log(cluster_summary$p_value, base = 10) < threshold_p_value,
                                        'insignificant',
                                        'significant')
  
  
  
  return(list('allocation' = allocation,
              'cluster_label' = cluster_label,
              'cluster_summary' = cluster_summary))
}


