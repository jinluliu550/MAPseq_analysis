# MCMC run with component probabilities dependent on x_i

#--------------------------------------- Simulation of allocations ---------------------------------------------

allocation_variables_dirmult_mcmc2 <- function(beta,
                                               delta,
                                               q_star,
                                               gamma_star,
                                               Y,
                                               x,
                                               num.cores){
  
  # Create an empty list
  M <- length(Y)
  C <- unlist(lapply(Y, ncol))
  J <- nrow(delta)
  R <- ncol(q_star)
  
  ## Initialize Z
  Z <- NULL
  allocation.prob <- NULL
  
  for(m in 1:M){
    
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)
    
    result.rbind <- foreach(c = 1:C[m],
                            .export = 'ddirmultinomial',
                            .combine = 'rbind',
                            .packages = 'extraDistr') %dopar% {
                              
                              ##-- log_p_tilde
                              log.p.tilde <- sapply(1:J, function(j){
                                
                                ddirmultinomial(Y[[m]][,c],
                                                n = sum(Y[[m]][,c]),
                                                alpha = q_star[j,]*gamma_star[j],
                                                log = TRUE)+
                                  log(delta[j,m]) + beta[j]*x[[m]][i]
                              })
                              
                              ##-- logK
                              logK <- -max(log.p.tilde)
                              
                              ##-- probability
                              allo.prob <- exp(log.p.tilde + logK)/sum(exp(log.p.tilde + logK))
                              
                              ##-- allocation
                              Z_cd <- extraDistr::rcat(n = 1,
                                                       prob = allo.prob)
                              
                              Z_cd_prob <- allo.prob[Z_cd]
                              
                              ##-- prepare for output
                              c(Z_cd, Z_cd_prob)
                              
                            }
    
    stopCluster(cl)
    
    ##-- Store Z
    Z[[m]] <- result.rbind[,1]
    
    ##-- Store allocation probability
    allocation.prob[[m]] <- result.rbind[,2]
    
  }
  
  ##-- Return Z
  return(list('Z' = Z,
              'allocation.prob' = allocation.prob))
}


#---------------------------------------- Simulation of xi ------------------------------------------

xi_mcmc <- function(delta,
                    beta,
                    x){
  
  # Dimensions
  M <- ncol(delta)
  J <- nrow(delta)
  C <- sapply(1:M,
              function(m) length(x[m]))
  
  
  output <- lapply(1:M,
                   function(m){
                     
                     rexp(n = C[m],
                          rate = colSums(matrix(rep(delta[,m], C[m])*exp(rep(beta, C[m])*rep(x[[m]], each = J)),
                                                nrow = J)))
                   })
  
  return(output)
  
}



#---------------------------------------- Simulation of delta ---------------------------------------

delta_mcmc <- function(alpha,
                       w,
                       z,
                       xi,
                       beta,
                       x){
  
  # Dimensions
  M <- length(x)
  J <- length(w)
  C <- sapply(1:M,
              function(m) length(x[[m]]))
  
  output <- lapply(1:M,
                   function(m){
                     
                     rgamma(n = J,
                            shape = alpha*w + as.numeric(table(factor(z[[m]], levels = 1:J))),
                            rate = 1+colSums(matrix(rep(xi[[m]], J)*exp(rep(beta, each = C[m])*rep(x[[m]],J)),
                                                    nrow = C[m])))
                   })
  
  # Return item: matrix with J rows and M columns
  return(do.call(cbind, output))
}





#---------------------------------------- Simulation of alpha ------------------------------------

alpha_log_prob2 <- function(w,
                            delta,
                            a_alpha,
                            b_alpha,
                            alpha){
  
  # Dimensions
  M <- ncol(delta)
  
  lprob <- dgamma(x = alpha,
                  shape = a_alpha,
                  rate = b_alpha) +
    
    sum(sapply(1:M),
        function(m){
          
          sum(-lgamma(alpha*w) + alpha*w*log(delta[,m]))
        })
  
  
  # output
  return(lprod)
}


alpha_mcmc2 <- function(w,
                        delta,
                        a_alpha,
                        b_alpha,
                        alpha,
                        X_mean,
                        M_2,
                        variance,
                        iter_num,
                        adaptive_prop){
  
  # Dimensions
  M <- ncol(delta)
  
  # Defining the inputs
  alpha_old <- alpha
  X_old <- log(alpha_old)
  X_mean_old <- X_mean
  M_2_old <- M_2
  variance_old <- variance
  
  # Defining the dimensions
  n <- iter_num
  
  if(n <= 100){
    
    X_new <- rnorm(n = 1,
                   mean = X_old,
                   sd = 1)
  }else{
    
    X_new <- rnorm(n = 1,
                   mean = X_old,
                   sd = sqrt(2.4^2*(variance_old + adaptive_prop)))
  }
  
  # alpha_new
  alpha_new <- exp(X_new)
  
  # Compute log acceptance probability
  log_acceptance <- alpha_log_prob2(w = w,
                                    delta = delta,
                                    a_alpha = a_alpha,
                                    b_alpha = b_alpha,
                                    alpha = alpha_new) -
    
    alpha_log_prob(w = w,
                   delta = delta,
                   a_alpha = a_alpha,
                   b_alpha = b_alpha,
                   alpha = alpha_old) +
    
    log(alpha_new) - log(alpha_old)
  
  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob = min(1, exp(log_acceptance)))
  
  if(outcome == 0){
    
    X_new <- X_old
    alpha_new <- alpha_old
  }
  
  # Update Covariance Structure
  X_mean_new <- (1-1/n)*X_mean_old + 1/n*X_new
  M_2_new <- M_2_old + (X_new-X_mean_old)*(X_new-X_mean_new)
  variance_new <- 1/(n-1)*M_2_new
  
  # Output
  return(list('alpha_new' = alpha_new,
              'X_mean_new' = X_mean_new,
              'M_2_new' = M_2_new,
              'variance_new' = variance_new,
              'accept' = outcome))
}




#--------------------------------------------- Simulation of beta ---------------------------------------------
b_beta_prior1 <- function(beta,
                          pi_cons,
                          s1_2,
                          s2_2){
  
  density_beta_prior <- pi_cons*dlaplace(x = beta,
                                         mu = 0,
                                         sigma = sqrt(s1_2))+
    
    (1-pi_cons)*dlaplace(x = beta,
                         mu = 0,
                         sigma = sqrt(s2_2))
  
  return(density_beta_prior)
}

# log probability
beta_log_prob <- function(beta,
                          x,
                          u,
                          xi,
                          delta,
                          pi_cons,
                          z,
                          s1_2,
                          s2_2,
                          num.cores){
  
  # Dimension
  J <- length(beta)
  M <- length(x)
  
  # log prob
  cl <- makeCluster(num.cores)
  registerDoParallel(cl)
  
  log_prob <- foreach(j = 1:J,
                      .export = 'b_beta_prior1',
                      .packages = 'extraDistr') %dopar% {
                        
                        
                        # subsets
                        subset.i <- lapply(1:M, 
                                           function(m) which(z[[m]]==j))
                        
                        subset.u <- lapply(1:M,
                                           function(m) u[[m]][subset.i[[m]], j])
                        
                        subset.xi <- lapply(1:M,
                                            function(m) xi[[m]][subset.i[[m]]])
                        
                        subset.delta <- sapply(1:M,
                                               function(m) delta[j,m])
                        
                        
                        
                        if(beta[j] < min(log(-log(unlist(subset.u))/(unlist(subset.xi)*rep(subset.delta, 
                                                                                           each = sapply(1:M,
                                                                                                         function(m) length(subset.i[[m]]))
                                                                                           )
                                                                     )
                                             )
                                         )
                           ){
                          
                          return_item <- beta[j]*length(which(unlist(x) == 1))*log(b_beta_prior1(beta = beta[j],
                                                                                                 pi_cons = pi_cons,
                                                                                                 s1_2 = s1_2,
                                                                                                 s2_2 = s2_2)
                          )
                          
                          }else{
                            
                            return_item <- 0
                            
                            }
                        
                        return(return_item)
                      }
  
  # Stop parallel computing
  stopCluster(cl)
  
  # Return
  return(unlist(log_prob))
}

beta_mcmc <- function(beta,
                      delta,
                      xi,
                      x,
                      M2,
                      z,
                      pi_cons,
                      beta_mean,
                      beta_variance,
                      iter_num,
                      adaptive_prop,
                      num.cores){
  
  # Record of previous statistics
  beta_old <- beta
  M2_old <- M2
  beta_mean_old <- beta_mean
  beta_variance_old <- beta_variance
  
  # Dimensions
  M <- ncol(delta)
  J <- nrow(delta)
  C <- sapply(1:M, 
              function(m) length(xi[[m]]))
  
  #-------------------------------------- First simulate u from an uniform distribution ---------------------------------------- 
  
  # c: List of length M, each item has C[m] rows and J columns
  c <- lapply(1:M, 
              function(m){
                
                c_IJ <- matrix(exp(-rep(xi[[m]], J)*rep(delta[,m], each = C[m])*exp(rep(beta, each = C[m])*rep(xi[[m]], J))),
                               nrow = C[m],
                               ncol = J)
                
                c_IJ
              })
  
  # u: List of length M, each item has C[m] rows and J columns
  u <- lapply(1:M,
              function(m){
                
                u_IJ <- matrix(runif(n = C[m]*J,
                                     min = 0,
                                     max = c_IJ),
                               nrow = C[m],
                               nrow = J)
                
                u_IJ
              })
  
  #-------------------------------- Second step: simulate from the posterior of beta given u
  n <- iter_num
  
  # Simulate new beta
  if(n <= 100){
    
    beta_new <- rnorm(n = J,
                      mean = beta_old,
                      sd = 1)
    
  }else{
    
    beta_new <- rnorm(n = J,
                      mean = beta_old,
                      sd = sqrt(2.4^2*(beta_variance_old + adaptive_prop))
    )
    
  }
  
  # Calculate acceptance probability
  log_acceptance <- beta_log_prob(beta = beta_new,
                                  x = x,
                                  u = u,
                                  xi = xi,
                                  delta = delta,
                                  pi_cons = pi_cons,
                                  z = z,
                                  num.cores = num.cores)-
    
    beta_log_prob(beta = beta_old,
                  x = x,
                  u = u,
                  xi = xi,
                  delta = delta,
                  pi_cons = pi_cons,
                  z = z,
                  num.cores = num.cores)
  
  
  # Random Bernoulli
  outcome <- rbinom(n = J,
                    size = 1,
                    prob = min(1, exp(log_acceptance)))
  
  if(outcome == 0){
    
    beta_new <- beta_old
  }
  
  # update variance-covariance
  beta_mean_new <- (1-1/n)*beta_mean_old + 1/n*beta_new
  M2_new <- M2_old + (beta_new-beta_mean_old)*(beta_new-beta_mean_new)
  variance_new <- 1/(n-1)*M2_new
  
  return(list('beta_new' = beta_new,
              'beta_mean_new' = beta_mean_new,
              'M2_new' = M2_new,
              'variance_new' = variance_new,
              'accept' = length(which(outcome==1)))
         )
  
}



#------------------------------------------------- Simulation of w -----------------------------------------------

w_log_prob <- function(w,
                       alpha,
                       alpha_zero,
                       delta){
  
  # Dimensions
  J <- length(w)
  M <- ncol(delta)
  
  # Log probability
  log_prob <- sum((alpha_zero/J-1)*log(w)) + sum(sapply(1:M,
                                                        function(m){
                                                          
                                                          sum((alpha*w-1)*log(delta[,m]) - delta[,m])
                                                        }))

  return(log_prob)
}

w_mcmc <- function(w,
                   delta,
                   alpha_zero,
                   alpha,
                   covariance,
                   mean_x,
                   tilde_s,
                   iter_num,
                   adaptive_prop = 0.01){
  
  
  # Dimensions
  J <- length(w)
  M <- ncol(delta)
  
  
  # Define the inputs
  w_old <- w
  covariance_old <- covariance
  mean_X_d_old <- mean_x
  tilde_s_old <- tilde_s
  
  #------------------------------------- Simulate new X
  X_old <- log(w_old[1:J-1]/w_old[J])
  
  n <- iter_num
  
  # Adaptive step
  if(n <= 100){
    
    X_new <- rmvnorm(n = 1,
                     mean = X_old,
                     sigma = diag(x=1, nrow=J-1, ncol=J-1))
    
  }else{
    
    X_new <- rmvnorm(n = 1,
                     mean = X_old,
                     sigma = 2.4^2/(J-1)*(covariance_old + adaptive_prop*diag(1, nrow = J-1, ncol = J-1)))
    
  }
  
  # Compute new w
  w_new <- c(exp(X_new)/(1+sum(exp(X_new))),1/(1+sum(exp(X_new))))
  
  
  # Compute acceptance probability
  log_acceptance <- w_log_prob(w = w_new,
                               delta = delta,
                               alpha_zero = alpha_zero,
                               alpha = alpha) -
    
    component_log_prob(w = w_old,
                       delta = delta,
                       alpha_zero = alpha_zero,
                       alpha = alpha) +
    
    sum(log(w_new)) - sum(log(w_old))
  
  
  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob=min(1, exp(log_acceptance)))
  
  if(outcome == 0){
    
    X_new <- X_old
    w_new <- w_old
  }
  
  
  
  # Update Covariance Structure
  tilde_s_new <- tilde_s_old + matrix(X_new, ncol = 1)%*%matrix(X_new, nrow = 1)
  mean_x_new <- mean_X_d_old*(1-1/n) + 1/n*matrix(X_new, nrow = 1)
  covariance_new <- 1/(n-1)*tilde_s_new - n/(n-1)*t(mean_x_new)%*%mean_x_new
  
  # Return a list of output
  return(list('w_new' = w_new,
              'tilde_s_new' = tilde_s_new,
              'mean_x_new' = mean_x_new,
              'covariance_new' = covariance_new,
              'accept' = outcome))
}


#------------------------------------------------- Simulation of alpha_zero ---------------------------------------------

alpha_zero_log_prob2 <- function(w,
                                 alpha_zero,
                                 a0,
                                 b0){
  
  # dimension
  J <- length(w)
  
  # log probability
  log_prob <- lgamma(alpha_zero) - J*lgamma(alpha_zero/J) + sum((alpha_zero/J-1)*log(w)) + (a0-1)*log(alpha_zero) - alpha_zero*b0
  
  # output
  return(log_prob)
}


alpha_zero_mcmc2 <- function(w,
                             alpha_zero,
                             a0,
                             b0,
                             X_mean,
                             M_2,
                             variance,
                             iter_num,
                             adaptive_prop){
  
  # dimension
  J <- length(w)
  n <- iter_num
  
  # previous values
  alpha_zero_old <- alpha_zero
  X_old <- log(alpha_zero_old)
  variance_old <- variance
  M_2_old <- M_2
  X_mean_old <- X_mean
  
  # Simulation
  if(n <= 100){
    
    X_new <- rnorm(n = 1,
                   mean = X_old,
                   sd = 1)
  }else{
    
    X_new <- rnorm(n = 1,
                   mean = X_old,
                   sd = sqrt(2.4^2*(variance_old + adaptive_prop)))
  }
  
  # alpha_zero_new
  alpha_zero_new <- exp(X_new)
  
  # log acceptance probability
  log_acceptance <- alpha_zero_log_prob2(w = w,
                                         alpha_zero = alpha_zero_new,
                                         a0 = a0,
                                         b0 = b0) -
    
    alpha_zero_log_prob(w = w,
                        alpha_zero = alpha_zero_old,
                        a0 = a0,
                        b0 = b0) +
    
    log(alpha_zero_new) - log(alpha_zero_old)
  
  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob = min(1, exp(log_acceptance)))
  
  if(outcome == 0){
    
    X_new <- X_old
    alpha_zero_new <- alpha_zero_old
  }
  
  # Update Covariance Structure
  X_mean_new <- (1-1/n)*X_mean_old + 1/n*X_new
  M_2_new <- M_2_old + (X_new-X_mean_old)*(X_new-X_mean_new)
  variance_new <- 1/(n-1)*M_2_new
  
  # Output
  return(list('alpha_zero_new' = alpha_zero_new,
              'X_mean_new' = X_mean_new,
              'M_2_new' = M_2_new,
              'variance_new' = variance_new,
              'accept' = outcome))
}




#---------------------------------------- Simulation of alpha_r ---------------------------------------------
alpha_r_logprob <- function(alpha_r,
                            q_star,
                            a_alpha,
                            b_alpha){
  
  R <- length(alpha_r)
  J <- nrow(q_star)
  
  # value <- -sum(alpha_r)
  value <- sum((a_alpha-1)*log(alpha_r)-b_alpha*alpha_r)
  
  for(j in 1:J){
    
    value <- value + lgamma(sum(alpha_r))+sum((alpha_r-1)*log(q_star[j,])-lgamma(alpha_r))
  }
  
  return(value)
}

alpha_r_mcmc <- function(alpha_r,
                         mean_X_alpha_r,
                         tilde_s_alpha_r,
                         q_star,
                         iter_num,
                         a_alpha,
                         b_alpha,
                         covariance,
                         adaptive_prop){
                          
                         
  
  ##-- Dimensions
  R <- length(alpha_r)
  J <- nrow(v)
  n <- iter_num
  
  ##-- Previous Covariance structure
  mean_X_alpha_r_old <- mean_X_alpha_r
  tilde_s_alpha_r_old <- tilde_s_alpha_r
  covariance_old <- covariance
  
  ##-- Previous alpha_r
  alpha_r_old <- alpha_r
  X_alpha_r_old <- log(alpha_r_old)
  
  ##-- Simulation
  if(n <= 100){
    
    X_alpha_r_new <- rmvnorm(n = 1,
                             mean = X_alpha_r_old,
                             sigma = diag(1,nrow = R, ncol = R))
  }else{
    
    X_alpha_r_new <- rmvnorm(n = 1,
                             mean = X_alpha_r_old,
                             sigma = (2.4^2)/2*(covariance_old + adaptive_prop*diag(1,nrow = R, ncol = R)))
  }
  
  ##-- alpha_r_new
  alpha_r_new <- exp(X_alpha_r_new)
  
  ##-- Acceptance probability
  accep.log.prob <- alpha_r_logprob(alpha_r = as.vector(alpha_r_new),
                                    q_star = q_star,
                                    a_alpha = a_alpha,
                                    b_alpha = b_alpha)-
    
    alpha_r_logprob(alpha_r = alpha_r_old,
                    q_star = q_star,
                    a_alpha = a_alpha,
                    b_alpha = b_alpha)-
    
    sum(log(alpha_r_old)) + sum(log(alpha_r_new))
  
  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob = min(1,exp(accep.log.prob)))
  
  if(outcome == 0){
    
    X_alpha_r_new <- X_alpha_r_old
    alpha_r_new <- alpha_r_old
  }
  
  ##-- Update Covariance Structure
  tilde_s_alpha_r_new <- tilde_s_alpha_r_old + matrix(X_alpha_r_new, ncol=1) %*%
    matrix(X_alpha_r_new, nrow=1)
  
  mean_X_alpha_r_new <- mean_X_alpha_r_old*(1-1/n) + 1/n*matrix(X_alpha_r_new, nrow=1)
  
  covariance_new <- 1/(n-1)*tilde_s_alpha_r_new - n/(n-1)*t(mean_X_alpha_r_new)%*%mean_X_alpha_r_new
  
  ##-- Output
  return(list('alpha_r_new' = alpha_r_new,
              'accept' = outcome,
              'tilde_s_alpha_r_new' = tilde_s_alpha_r_new,
              'mean_X_alpha_r_new' = mean_X_alpha_r_new,
              'covariance_alpha_r_new' = covariance_new))
}



#------------------------------------------ Simulation of gamma ------------------------------------

gamma_logprob2 <- function(Y,
                           Z,
                           q_j_star,
                           gamma_j_star,
                           a_gamma,
                           b_gamma,
                           j){
  
  
  # Bind Y by columns
  Y_bind <- do.call(cbind,Y)
  
  # Initialize the log probability - gamma prior
  log_prob <- (a_gamma-1)*log(gamma_j_star) - b_gamma*gamma_j_star
  
  # Find the cell index of the ones with allocation equal to j
  cell_index_j <- which(unlist(Z)==j)
  
  # Sum of neuron counts across all regions for these cells
  N_CM_j <- colSums(as.matrix(Y_bind[,cell_index_j]))
  
  # Number of regions
  R <- nrow(Y[[1]])
  
  # log - probability
  log_prob <- log_prob + sum(lgamma(sum(q_j_star*gamma_j_star)) + lgamma(N_CM_j + 1) - lgamma(N_CM_j + sum(q_j_star*gamma_j_star)) +
                               
                               colSums(lgamma(as.matrix(Y_bind[,cell_index_j]) + matrix(rep(q_j_star*gamma_j_star,
                                                                                            length(cell_index_j)),
                                                                                        nrow = R)
                               )) - sum(lgamma(q_j_star*gamma_j_star)) - colSums(lgamma(as.matrix(Y_bind[,cell_index_j] + 1)))
  )
  
  return(log_prob)
  
}

gamma_mcmc2 <- function(Y,
                        Z,
                        q_star,
                        gamma_star,
                        a_gamma,
                        b_gamma,
                        X_mean,
                        M_2,
                        variance,
                        iter_num,
                        adaptive_prop,
                        num.cores){
  
  J <- nrow(q_star)
  n <- iter_num
  R <- ncol(q_star)
  
  # previous values
  gamma_star_old <- gamma_star
  X_old <- log(gamma_star_old)
  variance_old <- variance
  M_2_old <- M_2
  X_mean_old <- X_mean
  
  
  # Register cores
  cl <- makeCluster(num.cores)
  registerDoParallel(cl)
  
  # For each component j
  loop.result <- foreach(j = 1:J,
                         .packages = c('stats',
                                       'base'),
                         .export = 'gamma_logprob2') %dopar% {
    
    
    # If component j is non-empty
    if(length(which(unlist(Z)==j)) != 0){
      
      if(n < 100){
        
        X_j_new <- rnorm(n = 1,
                         mean = X_old,
                         sd = 1)
        
        
      }else{
        
        X_j_new <- rnorm(n = 1,
                         mean = X_old,
                         sd = sqrt(2.4^2*(variance_old + adaptive_prop)))
        
      }
      
      # gamma_j_new
      gamma_j_new <- exp(X_j_new)
      
      # Acceptance probability
      log.accept.prob <- gamma_logprob(Y = Y,
                                       Z = Z,
                                       q_j_star = q_star[j,],
                                       gamma_j_star = gamma_j_new,
                                       a_gamma = a_gamma,
                                       b_gamma = b_gamma,
                                       j = j)-
        
        gamma_logprob(Y = Y,
                      Z = Z,
                      q_j_star = q_star[j,],
                      gamma_j_star = gamma_star_old[j],
                      a_gamma = a_gamma,
                      b_gamma = b_gamma,
                      j = j)-
        
        log(gamma_star_old[j]) + log(gamma_j_new)
      
      accept.prob <- exp(log.accept.prob)
      
    }else{
      
      # If component j is empty
      
      # Simulate from the prior
      gamma_j_new <- rgamma(n = 1,
                            shape = a_gamma,
                            rate = b_gamma)
      
      # X_new
      X_j_new <- log(gamma_j_new)
      
      # Always accept
      accept.prob <- 1
      
    }
    
    outcome <- rbinom(n = 1,
                      size = 1,
                      prob = min(1, accept.prob))
    
    # If outcome is to reject
    if(outcome == 0){
      
      gamma_j_new <- gamma_star_old[j]
      X_j_new <- X_old[j]
    }
    
    return(list('X_j_new'= X_j_new,
                'gamma_j_new' = gamma_j_new,
                'outcome' = outcome))
    
    }
  
  
  stopCluster(cl)
  
  
  
  # Update covariance structure
  X_mean_new <- sapply(1:J,
                       function(j) (1-1/n)*X_mean_old[j] + 1/n*loop.result[[j]]$X_j_new)
  
  M_2_new <- sapply(1:J,
                    function(j) M_2_old[j] + (loop.result[[j]]$X_j_new - X_mean_old[j])*(loop.result[[j]]$X_j_new - X_mean_new[j]))
  
  variance_new <- 1/(n-1)*M_2_new
  
  gamma_star_new <- sapply(1:J, 
                           function(j) loop.result[[j]]$gamma_j_new)
  
  gamma_count <- sum(sapply(1:J,
                            function(j) loop.result[[j]]$outcome))
  
  
  # Return
  return(list('gamma_star_new' = gamma_star_new,
              'variance_gamma_new' = variance_new,
              'M_2_gamma_new' = M_2_new,
              'X_mean_gamma_new' = X_mean_new,
              'gamma_count' = gamma_count))
  
}


#------------------------------------------ Simulation of q --------------------------------------------


q_star_logprob2 <- function(Y,
                            Z,
                            q_j_star,
                            gamma_j_star,
                            alpha_r,
                            j){
  
  
  # Bind Y by columns
  Y_bind <- do.call(cbind,Y)
  
  # Initialize the log probability - Dirichlet prior
  log_prob <- sum((alpha_r-1)*log(q_j_star))
  
  # Find the cell index of the ones with allocation equal to j
  cell_index_j <- which(unlist(Z)==j)
  
  # Sum of neuron counts across all regions for these cells
  N_CM_j <- colSums(as.matrix(Y_bind[,cell_index_j]))
  
  # Number of regions
  R <- nrow(Y[[1]])
  
  
  log_prob <- log_prob + sum(lgamma(sum(q_j_star*gamma_j_star)) + lgamma(N_CM_j + 1) - lgamma(N_CM_j + sum(q_j_star*gamma_j_star))) +
    
    sum(colSums(lgamma(as.matrix(Y_bind[,cell_index_j]) + matrix(rep(q_j_star*gamma_j_star,
                                                                     length(cell_index_j)),
                                                                 nrow = R)
    ))) - sum(lgamma(q_j_star*gamma_j_star))*length(cell_index_j) -
    
    sum(colSums(lgamma(as.matrix(Y_bind[,cell_index_j] + 1))))
  
  
  return(log_prob)
  
}

q_star_mcmc2 <- function(Y,
                         Z,
                         q_star,
                         gamma_star,
                         alpha_r,
                         covariance,
                         mean_x,
                         tilde_s,
                         iter_num,
                         adaptive_prop,
                         num.cores){
  
  J <- nrow(q_star)
  R <- ncol(q_star)
  n <- iter_num
  
  # Define the inputs
  q_star_old <- q_star
  covariance_old <- covariance
  mean_X_old <- mean_x
  tilde_s_old <- tilde_s
  X_old <- matrix(0, nrow = J, ncol = R-1) # Dimension = J x R-1
  
  for(j in 1:J){
    
    X_old[j,] <- log(q_star_old[j,1:(R-1)]/q_star_old[j,R])
  }
  
  
  # Register cores
  cl <- makeCluster(num.cores)
  registerDoParallel(cl)
  
  # For each component j
  loop.result <- foreach(j = 1:J,
                         .packages = c('base',
                                       'extraDistr',
                                       'mvtnorm'),
                         .export = 'q_star_logprob2'){
    
    
    # If component j is non-empty
    if(length(which(unlist(Z)==j)) != 0){
      
      if(n < 100){
        
        # Length of R - 1
        X_j_new <- as.vector(rmvnorm(n = 1,
                                     mean = X_old[j,],
                                     sigma = diag(1,nrow = R-1, ncol = R-1)))
        
      }else{
        
        # Length of R - 1
        X_j_new <- as.vector(rmvnorm(n = 1,
                                     mean = X_old[j,],
                                     sigma = (2.4^2)/2*(covariance_old[[j]] + adaptive_prop*diag(1,nrow = R-1, ncol = R-1))))
        
      }
      
      # q_star_new
      q_star_j_new <- c(exp(X_j_new)/(1+sum(exp(X_j_new))),1/(1+sum(exp(X_j_new))))
      
      # Acceptance probability
      log.accept.prob <- q_star_logprob2(Y = Y,
                                         Z = Z,
                                         q_j_star = q_star_j_new,
                                         gamma_j_star = gamma_star[j],
                                         alpha_r = alpha_r,
                                         j = j)-
        
        q_star_logprob(Y = Y,
                       Z = Z,
                       q_j_star = q_star_old[j,],
                       gamma_j_star = gamma_star[j],
                       alpha_r = alpha_r,
                       j = j)+
        
        sum(log(q_star_j_new)) - sum(log(q_star_old[j,]))
      
      
      accept.prob <- exp(log.accept.prob)
      
    }else{
      
      # If component j is empty
      
      # Simulate from the prior
      q_star_j_new <- as.vector(extraDistr::rdirichlet(n = 1,
                                                       alpha = alpha_r))
      
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
      
      q_star_j_new <- q_star_old[j,]
      X_j_new <- X_old[j,]
    }
    
    return(list('X_j_new' = X_j_new,
                'q_star_j_new' = q_star_j_new,
                'outcome' = outcome))
  }
  
  # Stop parallel computing
  stopCluster(cl)
  
  
  
  # Update covariance structure
  tilde_s_new <- lapply(1:J,
                        function(j) tilde_s_old[[j]] + matrix(loop.result[[j]]$X_j_new, ncol = 1)%*%matrix(loop.result[[j]]$X_j_new, nrow = 1))
  
  mean_x_new <- lapply(1:J,
                       function(j) (1-1/n)*mean_X_old[[j]] + 1/n*matrix(loop.result[[j]]$X_j_new, nrow = 1))
  
  covariance_new <- lapply(1:J,
                           function(j) 1/(n-1)*tilde_s_new[[j]] - n/(n-1)*t(mean_x_new[[j]])%*%mean_x_new[[j]])
  
  # q_star: matrix of dimension J x R
  q_star_new <- lapply(1:J,
                       matrix(loop.result[[j]]$q_star_j_new,
                              nrow = 1))
  
  q_star_new <- do.call(rbind, q_star_new)
  
  # counts of the total number of acceptance
  q_star_count <- sum(sapply(1:J,
                             function(j) loop.result[[j]]$outcome))
  
  
  
  
  # Return items
  return(list('tilde_s_new' = tilde_s_new,
              'mean_x_new' = mean_x_new,
              'covariance_new' = covariance_new,
              'q_star_count' = q_star_count,
              'q_star_new' = q_star_new))
  
  
}

