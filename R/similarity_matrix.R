
#' Posterior similarity matrix

similarity_matrix <- function(mcmc_run_all_output,
                              num.cores,
                              run.on.pc = TRUE){
  
  
  
  # Input from MCMC
  M <- mcmc_run_all_output$M
  C <- mcmc_run_all_output$C
  Z_trace <- mcmc_run_all_output$Z_output
  
  
  ##-- Cumulative sums of cells
  C_cum <- c(0, cumsum(C))
  
  ##-------------------- Similarity between cells from the same dataset and across datasets --------------
  
  psm.within <- NULL
  for(m in 1:M) psm.within[[m]] <- NULL
  
  ##-- For all possible combinations of datasets
  for(m1 in rev(1:M)){
    for(m2 in 1:m1){
      
      # Register Cores
      if(run.on.pc == FALSE){
        
        cl <- makeCluster(num.cores,
                          type = "FORK")
      }else{
        
        cl <- makeCluster(num.cores)
      }
      
      loop.result <- pblapply(1:length(Z_trace),
                              cl = cl,
                              FUN = function(iter){
                                
                                Z_trace_m1_iter <- Z_trace[[iter]][[m1]]
                                Z_trace_m2_iter <- Z_trace[[iter]][[m2]]
                                
                                psm.empty <- matrix(0, nrow = C[m1], ncol = C[m2])
                                for(i in 1:C[m1]){
                                  psm.empty[i,] <- psm.empty[i,] + ifelse(Z_trace_m2_iter == Z_trace_m1_iter[i],
                                                                          1,
                                                                          0)
                                }
                                
                                psm.empty
                                
                              })
      
      # Stop parallel computing
      stopCluster(cl)
      
      
      # Sum
      loop.result <- Reduce('+', loop.result)
      psm.within[[m1]][[m2]] <- loop.result/length(Z_trace)
      
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

