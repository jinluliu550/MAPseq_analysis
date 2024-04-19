
find_large_cluster <- function(mcmc_run_omega_output,
                               tau = 0.01,
                               top.n){
  
  # Trace of omega
  omega.trace <- mcmc_run_omega_output$omega_J_M_output
  
  # Number of clusters
  J <- nrow(omega.trace[[1]])
  
  # Posterior probability of having w_j > tau
  posterior.probability <- sapply(1:J,
                                  function(j){
                                    
                                    omega.trace.j <- sapply(1:length(omega.trace),
                                                            function(t) omega.trace[[t]][j])
                                    
                                    mean(omega.trace.j > tau)
                                  })
  
  return(which(rank(desc(posterior.probability)) <= top.n))
}
