

dlso_cluster_estimate <- function(mcmc_run_all_output){
  
  # trace of z
  z_trace <- mcmc_run_all_output$Z_output
  
  # trace of z in a matrix
  z_trace_mx <- lapply(1:length(z_trace),
                       function(t) matrix(unlist(z_trace[[t]]),
                                          nrow = 1))
  
  
  z_trace_mx <- do.call(rbind, z_trace_mx)
  
  z_point_estimate <- salso::dlso(truth = z_trace_mx)
  
  # Cut into data
  M <- mcmc_run_all_output$M
  C <- mcmc_run_all_output$C
  
  C_cumsum <- c(0, cumsum(C))
  
  z_point_estimate <- lapply(1:M,
                             function(m) z_point_estimate[(C_cumsum[m]+1):C_cumsum[m+1]])
  
  # Reorder
  z_point_estimate_r <- lapply(1:M,
                             function(m){
                               
                               sapply(1:C[m],
                                      function(c){
                                        
                                        which(sort(unique(unlist(z_point_estimate))) == z_point_estimate[[m]][c])
                                      })
                             })
  
  return(z_point_estimate_r)
}