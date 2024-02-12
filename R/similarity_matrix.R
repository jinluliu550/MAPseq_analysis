
#' Posterior similarity matrix
#'

similarity_matrix <- function(mcmc_run_all_output,
                              num.cores,
                              run.on.pc = TRUE){
  
  # Input from MCMC
  M <- mcmc_run_all_output$M
  C <- mcmc_run_all_output$C
  
  # Allocations
  allocation_list <- mcmc_run_all_output$Z_output
  
  
  
  ##-- Cumulative sums of cells
  C_cum <- c(0, cumsum(C))
  
  ##-------------------- Similarity between cells from the same dataset and across datasets --------------
  
  psm.within <- NULL
  for(m in 1:M) psm.within[[m]] <- NULL
  
  ##-- For all possible combinations of datasets
  for(m1 in rev(1:M)){
    for(m2 in 1:m1){
      
      
      
      # Allocation of mouse m1
      allocation_df_m1 <- lapply(1:length(allocation_list),
                                 function(t){
                                   
                                   allocation_list[[t]][[m1]]
                                 })
      
      allocation_df_m1 <- do.call(rbind, allocation_df_m1)
      
      # Allocation of mouse m2
      allocation_df_m2 <- lapply(1:length(allocation_list),
                                 function(t){
                                   
                                   allocation_list[[t]][[m2]]
                                 })
      
      allocation_df_m2 <- do.call(rbind, allocation_df_m2)
      
      # Combine the matrix
      allocation_df_m1_m2 <- cbind(allocation_df_m1,
                                   allocation_df_m2)
      
      # # Register Cores
      if(run.on.pc == FALSE){
        
        cl <- makeCluster(num.cores,
                          type = "FORK")
      }else{
        
        cl <- makeCluster(num.cores)
      }
      
      # Calculate similarity in each iteration
      similarity <- t(pbapply(allocation_df_m1_m2,
                              1,
                              function(x) ifelse(rep(x[1:C[m1]], each = C[m2]) == rep(x[(C[m1]+1):(C[m1]+C[m2])], C[m1]),
                                                 1,
                                                 0),
                              cl = cl)
      )
      
      stopCluster(cl)
      
      # Combine everything
      psm.within[[m1]][[m2]] <- matrix(colSums(similarity)/length(allocation_list),
                                       nrow = C[m1],
                                       ncol = C[m2],
                                       byrow = TRUE)
      
      
    }
  }
  
  ##-------------------------------- Computing a combined matrix -------------------------
  
  combined_matrix <- matrix(0, nrow = max(C_cum), ncol = max(C_cum))
  
  for(m1 in 1:M){
    for(m2 in 1:m1){
      
      combined_matrix[(C_cum[m1]+1):C_cum[m1+1], (C_cum[m2]+1):C_cum[m2+1]] <- psm.within[[m1]][[m2]]
    }
    
    if(m1 != M){
      for(m2 in (m1+1):M){
        
        combined_matrix[(C_cum[m1]+1):C_cum[m1+1], (C_cum[m2]+1):C_cum[m2+1]] <- t(psm.within[[m2]][[m1]])
      }
    }
  }
  
  ##-------------------------------- Return both -----------------------------------------
  
  return(list('psm.within' = psm.within,
              'psm.combined' = combined_matrix))
  
}
