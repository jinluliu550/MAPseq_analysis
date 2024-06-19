#' Posterior predictive check with single replicated data
#'

add_noise <- function(mcmc_run_all_output,
                      Y,
                      regions.name = NULL){
  
  
  
  length_per_chain <- length(mcmc_run_all_output$Z_output)
  
  C <- mcmc_run_all_output$C
  R <- mcmc_run_all_output$R
  M <- mcmc_run_all_output$M
  
  if(is.null(regions.name)){
    
    regions.name <- paste('region', 1:R)
  }
  
  
  # Sum of counts for each neuron cell
  N_CM <- lapply(1:M,
                 function(m) colSums(Y[[m]]))
  
  # theta
  target.index <- sample(length_per_chain,
                         size = 1,
                         replace = FALSE)
  
  
  theta <- list(q_j_star = mcmc_run_all_output$q_star_1_J_output[[target.index]],
                gamma_j_star = mcmc_run_all_output$gamma_star_1_J_output[[target.index]],
                Z = mcmc_run_all_output$Z_output[[target.index]])
  
  # Simulate replicated data Y
  
  Y_rep <- lapply(1:M,
                  function(m){
                    
                    Y_rep_d <- do.call(cbind, lapply(1:C[m],
                                                     function(c){
                                                       
                                                       prob.c <- rdirichlet(n = 1,
                                                                            theta$q_j_star[theta$Z[[m]][c],]*theta$gamma_j_star[theta$Z[[m]][c]])
                                                       
                                                       rmultinom(n = 1, size = N_CM[[m]][c], prob = prob.c)
                                                     }))
                    
                    rownames(Y_rep_d) <- regions.name
                    Y_rep_d
                  })
  
  # Add some noise to the Y_rep
  noisy_data <- lapply(1:M,
                       function(m){
                         
                         noisy_data_m <- lapply(1:R,
                                                function(r){
                                                  
                                                  noisy_data_m0 <- lapply(1:ncol(Y[[m]]),
                                                                          function(i){
                                                                            
                                                                            prob_mi <- rep(0, R)
                                                                            prob_mi[r] <- 0.9
                                                                            prob_mi[-r] <- 0.1/(R-1)
                                                                            
                                                                            matrix(rmultinom(n = 1, size = Y[[m]][r,i], prob = prob_mi),
                                                                                   ncol = 1)
                                                                          })
                                                  
                                                  do.call(cbind, noisy_data_m0)
                                                })
                         
                         # Add up and round to integer
                         noisy_data_m <- round(Reduce('+', noisy_data_m)/R, 0)
                         rownames(noisy_data_m) <- regions.name
                         
                         noisy_data_m
                       })
  
  
  return(list('noisy_data' = noisy_data,
              'replicated_data' = Y_rep))
  
  
  
}
