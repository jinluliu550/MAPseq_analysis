
##-------------------------- q_star ------------------------------


q_star_logprob <- function(Y,
                           Z,
                           q_j_star,
                           gamma_j_star,
                           alpha_h,
                           j){
  
  
  # Bind Y by columns
  Y_bind <- do.call(cbind,Y)
  
  
  # Find the cell index of the ones with allocation equal to j
  cell_index_j <- which(unlist(Z)==j)
  
  # Sum of neuron counts across all regions for these cells
  N_CM_j <- colSums(as.matrix(Y_bind[,cell_index_j]))
  
  # Number of regions
  R <- nrow(Y[[1]])
  
  # Initialize the log probability - Dirichlet prior
  log_prob <- sum((alpha_h-1)*log(q_j_star))
  
  # log probability - distribution of Y
  log_prob <- log_prob + sum(colSums(lgamma(as.matrix(Y_bind[,cell_index_j]) + matrix(rep(q_j_star*gamma_j_star,
                                                                                          length(cell_index_j)),
                                                                                      nrow = R)
  ))) - sum(lgamma(q_j_star*gamma_j_star))*length(cell_index_j) -
    
    sum(colSums(lgamma(as.matrix(Y_bind[,cell_index_j] + 1))))
  
  
  return(log_prob)
  
}

q_star_mcmc_f <- function(Y,
                        Z,
                        q_star_1_J,
                        gamma_1_J_star,
                        alpha_h,
                        covariance,
                        mean_x,
                        tilde_s,
                        iter_num,
                        adaptive_prop){
  
  J <- nrow(q_star_1_J)
  R <- ncol(q_star_1_J)
  n <- iter_num
  
  # Define the inputs
  q_star_1_J_old <- q_star_1_J
  covariance_old <- covariance
  mean_X_old <- mean_x
  tilde_s_old <- tilde_s
  X_old <- matrix(0, nrow = J, ncol = R-1) # Dimension = J x R-1
  for(j in 1:J){
    
    X_old[j,] <- log(q_star_1_J_old[j,1:(R-1)]/q_star_1_J_old[j,R])
  }
  
  
  # Register cores
  cl <- makeCluster(num.cores)
  registerDoParallel(cl)
  
  # For each component
  loop.result <- foreach(j = 1:J,
                         .packages = c('mvtnorm',
                                       'extraDistr',
                                       'stats'),
                         .export = c('q_star_logprob')) %dopar% {
                           
                           
                           # If component j is non-empty
                           
                           if(length(which(unlist(Z)==j)) != 0){
                             
                             if(n < 100){
                               
                               # Length of R - 1
                               X_j_new <- as.vector(rmvnorm(n = 1,
                                                            mean = X_old[j,],
                                                            sigma = diag(0.01,nrow = R-1, ncol = R-1)))
                               
                             }else{
                               
                               # Length of R - 1
                               X_j_new <- as.vector(rmvnorm(n = 1,
                                                            mean = X_old[j,],
                                                            sigma = (2.4^2)/(R-1)*(covariance_old[[j]] + adaptive_prop*diag(1,nrow = R-1, ncol = R-1))))
                               
                             }
                             
                             # q_star_new
                             q_star_j_new <- c(exp(X_j_new)/(1+sum(exp(X_j_new))),1/(1+sum(exp(X_j_new))))
                             
                             # Acceptance probability
                             log.accept.prob <- q_star_logprob(Y = Y,
                                                               Z = Z,
                                                               q_j_star = q_star_j_new,
                                                               gamma_j_star = gamma_1_J_star[j],
                                                               alpha_h = alpha_h,
                                                               j = j)-
                               
                               q_star_logprob(Y = Y,
                                              Z = Z,
                                              q_j_star = q_star_1_J_old[j,],
                                              gamma_j_star = gamma_1_J_star[j],
                                              alpha_h = alpha_h,
                                              j = j)+
                               
                               sum(log(q_star_j_new)) - sum(log(q_star_1_J_old[j,]))
                             
                             
                             accept.prob <- exp(log.accept.prob)
                             
                           }else{
                             
                             # If component j is empty
                             
                             # Simulate from the prior
                             q_star_j_new <- as.vector(extraDistr::rdirichlet(n = 1,
                                                                              alpha = alpha_h))
                             
                             # X_new
                             X_j_new <- log(q_star_j_new[1:(R-1)]/q_star_j_new[R])
                             
                             # Always accept
                             accept.prob <- 1
                             
                           }
                           
                           # Random Bernoulli
                           outcome <- rbinom(n = 1,
                                             size = 1,
                                             prob = min(1,accept.prob))
                           
                           # If outcome is to reject
                           if(outcome == 0){
                             
                             q_star_j_new <- q_star_1_J_old[j,]
                             X_j_new <- X_old[j,]
                             
                           }
                           
                           tilde_s_new.j <- tilde_s_old[[j]] + matrix(X_j_new, ncol = 1)%*%matrix(X_j_new, nrow = 1)
                           mean_x_new.j <- mean_X_old[[j]]*(1-1/n) + 1/n*matrix(X_j_new, nrow = 1)
                           covariance_new.j <- 1/(n-1)*tilde_s_new.j - n/(n-1)*t(mean_x_new.j)%*%mean_x_new.j
                           
                           list('tilde_s_new.j' = tilde_s_new.j,
                                'mean_x_new.j' = mean_x_new.j,
                                'covariance_new.j' = covariance_new.j,
                                'q_star_j_new' = q_star_j_new,
                                'outcome' = outcome)
                           
                         }
  
  # Stop cluster
  stopCluster(cl)
  
  # Prepare for output
  tilde_s_new <- lapply(1:J, 
                        function(j) loop.result$tilde_s_new.j)
  
  mean_x_new <- lapply(1:J, 
                       function(j) loop.result$mean_x_new.j)
  
  covariance_new <- lapply(1:J,
                           function(j) loop.result$covariance_new.j)
  
  q_star_j_new <- lapply(1:J,
                         function(j) matrix(loop.result$q_star_j_new, nrow = 1))
  
  q_star_1_J_new <- do.call(rbind, q_star_j_new)
  
  q_star_count <- sum(sapply(1:J,
                             function(j) loop.result$outcome))
  
  
  # Return
  return(list('tilde_s_new' = tilde_s_new,
              'mean_x_new' = mean_x_new,
              'covariance_new' = covariance_new,
              'q_star_count' = q_star_count,
              'q_star_1_J_new' = q_star_1_J_new))
  
  
}

##---------------------------- alpha_h -------------------------------------


alpha_h_logprob <- function(alpha_h,
                            q_star_1_J,
                            a_alpha,
                            b_alpha){
  
  R <- length(alpha_h)
  J <- nrow(q_star_1_J)
  
  # value <- -sum(alpha_h)
  value <- sum((a_alpha-1)*log(alpha_h)-b_alpha*alpha_h)
  
  for(j in 1:J){
    
    value <- value + lgamma(sum(alpha_h))+sum((alpha_h-1)*log(q_star_1_J[j,])-lgamma(alpha_h))
  }
  
  return(value)
}

alpha_h_mcmc <- function(alpha_h,
                         mean_X_alpha_h,
                         tilde_s_alpha_h,
                         q_star_1_J,
                         iter_num,
                         a_alpha,
                         b_alpha,
                         covariance,
                         adaptive_prop){
  
  ##-- Dimensions
  R <- length(alpha_h)
  J <- nrow(q_star_1_J)
  n <- iter_num
  
  ##-- Previous Covariance structure
  mean_X_alpha_h_old <- mean_X_alpha_h
  tilde_s_alpha_h_old <- tilde_s_alpha_h
  covariance_old <- covariance
  
  ##-- Previous alpha_h
  alpha_h_old <- alpha_h
  X_alpha_h_old <- log(alpha_h_old)
  
  ##-- Simulation
  if(n <= 100){
    
    X_alpha_h_new <- rmvnorm(n = 1,
                             mean = X_alpha_h_old,
                             sigma = diag(0.01,nrow = R, ncol = R))
  }else{
    
    X_alpha_h_new <- rmvnorm(n = 1,
                             mean = X_alpha_h_old,
                             sigma = (2.4^2)/R*(covariance_old + adaptive_prop*diag(1,nrow = R, ncol = R)))
  }
  
  ##-- alpha_h_new
  alpha_h_new <- exp(X_alpha_h_new)
  
  ##-- Acceptance probability
  accep.log.prob <- alpha_h_logprob(alpha_h = as.vector(alpha_h_new),
                                    q_star_1_J = q_star_1_J,
                                    a_alpha = a_alpha,
                                    b_alpha = b_alpha)-
    
    alpha_h_logprob(alpha_h = as.vector(alpha_h_old),
                    q_star_1_J = q_star_1_J,
                    a_alpha = a_alpha,
                    b_alpha = b_alpha)-
    
    sum(log(alpha_h_old)) + sum(log(alpha_h_new))
  
  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob = min(1,exp(accep.log.prob)))
  
  if(outcome == 0){
    
    X_alpha_h_new <- X_alpha_h_old
    alpha_h_new <- alpha_h_old
  }
  
  ##-- Update Covariance Structure
  tilde_s_alpha_h_new <- tilde_s_alpha_h_old + matrix(X_alpha_h_new, ncol=1) %*%
    matrix(X_alpha_h_new, nrow=1)
  
  mean_X_alpha_h_new <- mean_X_alpha_h_old*(1-1/n) + 1/n*matrix(X_alpha_h_new, nrow=1)
  
  covariance_new <- 1/(n-1)*tilde_s_alpha_h_new - n/(n-1)*t(mean_X_alpha_h_new)%*%mean_X_alpha_h_new
  
  ##-- Output
  return(list('alpha_h_new' = alpha_h_new,
              'accept' = outcome,
              'tilde_s_alpha_h_new' = tilde_s_alpha_h_new,
              'mean_X_alpha_h_new' = mean_X_alpha_h_new,
              'covariance_alpha_h_new' = covariance_new))
}

