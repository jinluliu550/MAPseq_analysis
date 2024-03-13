
# Function to reorder clusters
binom_cluster_reorder <- function(Y,
                                  Z){
  
  # Calculate proportions
  Y_proportion <- Y
  N_data <- rowSums(Y)
  
  for(i in 1:ncol(Y)){
    
    Y_proportion[,i] <- Y_proportion[,i]/N_data
  }
  
  # Calculate the average projection strength of neurons in each cluster
  average_q <- lapply(sort(unique(Z)),
                      function(j){
                        
                        data_j <- matrix(Y_proportion[which(Z == j),],
                                            
                                         nrow = length(which(Z==j)))
                        
                        matrix(colMeans(data_j),
                               nrow = 1)
                      })
  
  average_q <- do.call(rbind, average_q)
  
  # Reorder clusters
  strong.proj <- apply(average_q,
                       1,
                       function(x) which(x == max(x)))
  
  
  df0 <- lapply(1:ncol(average_q),
                function(i){
                  
                  # Clusters with the strongest projection to region i
                  strong.proj.i <- which(strong.proj==i)
                  
                  if(length(strong.proj.i) != 0){
                    
                    qj_estimate_i <- matrix(average_q[strong.proj.i,],
                                            nrow = length(strong.proj.i))
                    
                    output <- strong.proj.i[order(-qj_estimate_i[,i])]
                  }else{
                    
                    output <- NULL
                  }
                  
                }
  )
  
  old_ordering <- unlist(df0)
  
  Z_new <- sapply(1:length(Z),
                  function(i) which(old_ordering == Z[i]))
  
  return(list('Z' = Z_new,
              'old_ordering' = old_ordering))
}

