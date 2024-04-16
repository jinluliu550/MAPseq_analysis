
##---------------- MCMC Simulation for Dirichlet Multinomial --------------

ddirmultinomial <- function(x, n, alpha, log = FALSE){

  # ind = alpha>0
  # 
  # # zero density if x_r>0 and alpha_r=0
  # if(sum(x[!ind])>0){ log.prob = -Inf}
  # else{
  #   log.prob <- lgamma(sum(alpha[ind])) + lgamma(n+1) - lgamma(n+sum(alpha[ind])) + sum(lgamma(x[ind] + alpha[ind]) - lgamma(alpha[ind]) - lgamma(x[ind]+1))
  # }  
  log.prob <- lgamma(sum(alpha)) + lgamma(n+1) - lgamma(n+sum(alpha)) + sum(lgamma(x + alpha) - lgamma(alpha) - lgamma(x+1))
  
  if(isFALSE(log)){

    return(exp(log.prob))
    
  }else{

    return(log.prob)
  }
}


##--------------- Simulation of allocations ----------------
allocation_variables_dirmult_mcmc <- function(omega_J_M,
                                              q_star_1_J,
                                              gamma_1_J_star,
                                              Y,
                                              num.cores){

  # Create an empty list
  M <- length(Y)
  C <- unlist(lapply(Y, ncol))
  J <- nrow(omega_J_M)
  R <- ncol(q_star_1_J)

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
                                               alpha = q_star_1_J[j,]*gamma_1_J_star[j],
                                               log = TRUE)+
                                 log(omega_J_M[j,m])
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

##------------ Simulation of Data-Specific component probabilities --------------------
dataset_specific_mcmc <- function(Z,
                                  omega,
                                  alpha){

  J <- length(omega)
  M <- length(Z)

  loop.result <- lapply(1:M, function(m) {

    parameter <- as.vector(table(factor(Z[[m]], levels = c(1:J)))) + alpha*omega
    
    omega_J_M <- extraDistr::rdirichlet(n = 1,
                                             alpha = parameter)
    

    # shifted the zero values
    omega_J_M <- ifelse(omega_J_M == 0, .Machine$double.eps, omega_J_M)
    omega_J_M <- omega_J_M/sum(omega_J_M)

    return(omega_J_M)
  })

  omega_J_M <- matrix(unlist(loop.result), nrow=J, ncol=M)

  # Returning omega_J_M
  return(omega_J_M)
}



##----------- Simulation of Component probabilities ----------------
component_log_prob <- function(omega,
                               omega_J_M,
                               alpha_zero,
                               alpha){

  J <- nrow(omega_J_M)
  M <- ncol(omega_J_M)

  lprod <- sum((alpha_zero/J-1)*log(omega)) +

    sum(unlist(lapply(1:M, function(m){

      sum(alpha*omega*log(omega_J_M[,m])-lgamma(alpha*omega))
    })))

  # outputs
  return(lprod)
}

component_probabilities_mcmc <- function(omega,
                                         omega_J_M,
                                         alpha_zero,
                                         alpha,
                                         covariance,
                                         mean_x,
                                         tilde_s,
                                         iter_num,
                                         adaptive_prop = 0.01){

  J <- nrow(omega_J_M)
  M <- ncol(omega_J_M)

  # Define the inputs
  omega_old <- omega
  covariance_old <- covariance
  mean_X_d_old <- mean_x
  tilde_s_old <- tilde_s

  X_old <- log(omega_old[1:J-1]/omega_old[J]) # Length = J-1

  # Specify the iteration number
  n <- iter_num

  # Adaptive step
  if(n <= 100){

    X_new <- rmvnorm(n = 1,
                     mean = X_old,
                     sigma = diag(x=0.01, nrow=J-1, ncol=J-1))

  }else{

    X_new <- rmvnorm(n = 1,
                     mean = X_old,
                     sigma = 2.4^2/(J-1)*(covariance_old + adaptive_prop*diag(1, nrow = J-1, ncol = J-1)))

  }

  # Compute omega_new (Length = J) from X_new
  omega_new <- c(exp(X_new)/(1+sum(exp(X_new))),1/(1+sum(exp(X_new))))
  
  # shifted the zero values
  omega_new <- ifelse(omega_new == 0, .Machine$double.eps, omega_new)
  omega_new <- omega_new/sum(omega_new)

  # Compute acceptance probability
  log_acceptance <- component_log_prob(omega_new,
                                       omega_J_M,
                                       alpha_zero,
                                       alpha) -
    component_log_prob(omega_old,
                       omega_J_M,
                       alpha_zero,
                       alpha) +

    sum(log(omega_new)) - sum(log(omega_old))
  
  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob=min(1, exp(log_acceptance)))


  if(outcome == 0){

    X_new <- X_old
    omega_new <- omega_old
  }

  # Update Covariance Structure
  tilde_s_new <- tilde_s_old + matrix(X_new, ncol = 1)%*%matrix(X_new, nrow = 1)
  mean_x_new <- mean_X_d_old*(1-1/n) + 1/n*matrix(X_new, nrow = 1)
  covariance_new <- 1/(n-1)*tilde_s_new - n/(n-1)*t(mean_x_new)%*%mean_x_new

  # Output
  return(list('omega_new' = omega_new,
              'tilde_s_new' = tilde_s_new,
              'mean_x_new' = mean_x_new,
              'covariance_new' = covariance_new,
              'accept' = outcome))
}


## ----------------------- Simulation of Concentration parameter alpha -----------------------
alpha_log_prob <- function(omega_J_M,
                           omega,
                           alpha){

  M <- ncol(omega_J_M)

  lprod <- -alpha + M*lgamma(alpha) +

    sum(unlist(lapply(1:M, function(m) {
      sum(alpha*omega*log(omega_J_M[,m])-lgamma(alpha*omega))
    })))

  # output
  return(lprod)
}


alpha_mcmc <- function(omega_J_M,
                       omega,
                       alpha,
                       X_mean,
                       M_2,
                       variance,
                       iter_num,
                       adaptive_prop){

  M <- ncol(omega_J_M)

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
                   sd = 0.01)
  }else{

    X_new <- rnorm(n = 1,
                   mean = X_old,
                   sd = sqrt(2.4^2*(variance_old + adaptive_prop)))
  }

  # alpha_new
  alpha_new <- exp(X_new)

  # Compute log acceptance probability
  log_acceptance <- alpha_log_prob(omega_J_M = omega_J_M,
                                   omega = omega,
                                   alpha = alpha_new) -

    alpha_log_prob(omega_J_M = omega_J_M,
                   omega = omega,
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

##---------------------- Simulation of concentration parameter alpha_zero -----------------
alpha_zero_log_prob <- function(omega,
                                alpha_zero){

  # dimension
  J <- length(omega)

  # log probability
  lprob <- -alpha_zero + lgamma(alpha_zero) - J*lgamma(alpha_zero/J) + sum(alpha_zero/J*log(omega))

  # output
  return(lprob)
}

alpha_zero_mcmc <- function(omega,
                            alpha_zero,
                            X_mean,
                            M_2,
                            variance,
                            iter_num,
                            adaptive_prop){

  # dimension
  J <- length(omega)
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
                   sd = 0.01)
  }else{

    X_new <- rnorm(n = 1,
                   mean = X_old,
                   sd = sqrt(2.4^2*(variance_old + adaptive_prop)))
  }

  # alpha_zero_new
  alpha_zero_new <- exp(X_new)

  # log acceptance probability
  log_acceptance <- alpha_zero_log_prob(omega = omega,
                                        alpha_zero = alpha_zero_new) -

    alpha_zero_log_prob(omega = omega,
                        alpha_zero = alpha_zero_old) +

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

##------------------------------ Gamma ----------------------------------

gamma_logprob <- function(Y,
                          Z,
                          q_j_star,
                          gamma_j_star,
                          a_gamma,
                          b_gamma,
                          j){
  
  # Number of regions
  R <- nrow(Y[[1]])

  # Bind Y by columns
  Y_bind <- do.call(cbind,Y)

  # Initialize the log probability - gamma prior
  log_prob <- (a_gamma-1)*log(gamma_j_star) - b_gamma*gamma_j_star

  # Find the cell index of the ones with allocation equal to j
  cell_index_j <- which(unlist(Z)==j)

  # Sum of neuron counts across all regions for these cells
  N_CM_j <- colSums(matrix(Y_bind[,cell_index_j],
                           nrow = R,
                           ncol = length(cell_index_j)))


  # log - probability
  log_prob <- log_prob + sum(lgamma(gamma_j_star) + lgamma(N_CM_j + 1) - lgamma(N_CM_j + gamma_j_star) +

                               colSums(lgamma(matrix(Y_bind[,cell_index_j],
                                                     nrow = R,
                                                     ncol= length(cell_index_j)) + matrix(rep(q_j_star*gamma_j_star,
                                                                                            length(cell_index_j)),
                                                                                        nrow = R)
                               )) - sum(lgamma(q_j_star*gamma_j_star)) - colSums(lgamma(matrix(Y_bind[,cell_index_j] + 1,
                                                                                               nrow = R,
                                                                                               ncol = length(cell_index_j)))
                                                                                 )
  )

  return(log_prob)

}

gamma_mcmc <- function(Y,
                       Z,
                       q_star_1_J,
                       gamma_1_J_star,
                       a_gamma,
                       b_gamma,
                       X_mean,
                       M_2,
                       variance,
                       iter_num,
                       adaptive_prop){

  J <- nrow(q_star_1_J)
  n <- iter_num
  R <- ncol(q_star_1_J)

  # previous values
  gamma_1_J_star_old <- gamma_1_J_star
  X_old <- log(gamma_1_J_star_old)
  variance_old <- variance
  M_2_old <- M_2
  X_mean_old <- X_mean

  # Update
  X_mean_new <- NULL
  M_2_new <- NULL
  variance_new <- NULL
  gamma_count <- 0
  gamma_1_J_star_new <- NULL
  
  # For each component j
  for(j in 1:J){
    
    
    if(length(which(unlist(Z)==j)) != 0){
      
      if(n < 100){
        
        X_j_new <- rnorm(n = 1,
                         mean = X_old,
                         sd = 0.01)
        
        
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
                                       q_j_star = q_star_1_J[j,],
                                       gamma_j_star = gamma_j_new,
                                       a_gamma = a_gamma,
                                       b_gamma = b_gamma,
                                       j = j)-
        
        gamma_logprob(Y = Y,
                      Z = Z,
                      q_j_star = q_star_1_J[j,],
                      gamma_j_star = gamma_1_J_star_old[j],
                      a_gamma = a_gamma,
                      b_gamma = b_gamma,
                      j = j)-
        
        log(gamma_1_J_star_old[j]) + log(gamma_j_new)
      
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
 
    # Random Bernoulli
    outcome <- rbinom(n = 1,
                      size = 1,
                      prob = min(1,accept.prob))

    # If outcome is to reject
    if(outcome == 0){
      
      gamma_j_new <- gamma_1_J_star_old[j]
      X_j_new <- X_old[j]
    }
    
    X_mean_new[j] <- (1-1/n)*X_mean_old[j] + 1/n*X_j_new
    M_2_new[j] <- M_2_old[j] + (X_j_new-X_mean_old[j])*(X_j_new-X_mean_new[j])
    variance_new[j] <- 1/(n-1)*M_2_new[j]
    gamma_count <- gamma_count + outcome
    gamma_1_J_star_new[j] <- gamma_j_new
    
  }
  


  # Return
  return(list('gamma_1_J_star_new' = gamma_1_J_star_new,
              'variance_gamma_new' = variance_new,
              'M_2_gamma_new' = M_2_new,
              'X_mean_gamma_new' = X_mean_new,
              'gamma_count' = gamma_count))

}

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

q_star_mcmc <- function(Y,
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

  q_star_1_J_new <- matrix(0, J, R)
  tilde_s_new <- NULL
  mean_x_new <- NULL
  covariance_new <- NULL
  q_star_count <- 0
  
  # For each cluster
  for(j in 1:J){
    
    
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
      # shifted the zero values
      q_star_j_new <- ifelse( q_star_j_new == 0, .Machine$double.eps,  q_star_j_new)
      q_star_j_new <-  q_star_j_new/sum( q_star_j_new)
      
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
      # shifted the zero values
      q_star_j_new <- ifelse( q_star_j_new == 0, .Machine$double.eps,  q_star_j_new)
      q_star_j_new <-  q_star_j_new/sum( q_star_j_new)
      
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
    
    # Update
    tilde_s_new[[j]] <- tilde_s_old[[j]] + matrix(X_j_new, ncol = 1)%*%matrix(X_j_new, nrow = 1)
    mean_x_new[[j]] <- mean_X_old[[j]]*(1-1/n) + 1/n*matrix(X_j_new, nrow = 1)
    covariance_new[[j]] <- 1/(n-1)*tilde_s_new[[j]] - n/(n-1)*t(mean_x_new[[j]])%*%mean_x_new[[j]]
    q_star_count <- q_star_count + outcome
    q_star_1_J_new[j,] <- q_star_j_new
  }
  
  
  
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
  # shifted the zero values
  alpha_h_new <- ifelse( alpha_h_new == 0, .Machine$double.eps,  alpha_h_new)
  

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

