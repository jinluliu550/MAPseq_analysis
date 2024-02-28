
mcmc_run_all2 <- function(Y,
                          x,
                          J,
                          number_iter = 6000,
                          thinning = 5,
                          burn_in = 3000,
                          adaptive_prop = 0.1,
                          print_Z = FALSE,
                          iter_update = 100,
                          a_gamma,
                          b_gamma,
                          a0 = 1,
                          b0 = 1,
                          a_alpha = 1,
                          b_alpha = 1,
                          s,
                          num.cores = 1){
  
  
  
  ##--------------- dimensions --------------------
  
  M <- length(Y)
  
  ##-- Need to make sure that all datasets have the same number of regions
  stopifnot(length(unique(sapply(1:M, function(m) nrow(Y[[m]]))))==1)
  
  ##-- Number of regions
  R <- nrow(Y[[1]])
  
  ##-- Number of neuron cells - vector of length M
  C <- unlist(lapply(1:M, function(m) ncol(Y[[m]])))
  
  
  #------------------------ Prepare for outputs -----------------
  
  Z_output <- NULL
  xi_output <- NULL
  w_output <- NULL
  delta_output <- NULL
  alpha_output <- 0
  alpha_zero_output <- 0
  beta_output <- NULL
  alpha_r_output <- NULL
  gamma_star_output <- NULL
  q_star_output <- NULL
  allocation_prob_output <- NULL
  w_m_0_output <- NULL
  w_m_1_output <- NULL
  
  #------------------------ List of acceptance probabilities and counts --------------------------
  
  acceptance_prob_list <- data.frame(w = rep(0, number_iter-1),
                                     alpha_zero = rep(0, number_iter-1),
                                     alpha = rep(0, number_iter-1),
                                     q_star = rep(0, number_iter-1),
                                     gamma_star = rep(0, number_iter-1),
                                     alpha_r = rep(0, number_iter-1),
                                     tau = rep(0, number_iter-1))
  
  alpha_count <- 0
  alpha_zero_count <- 0
  w_count <- 0
  alpha_r_count <- 0
  q_star_count <- 0
  gamma_star_count <- 0
  tau_count <- 0
  
  
  #----------------------- Step 2: Set starting values ----------------------------
  
  # beta
  beta <- rep(1,J)
  
  # delta
  delta <- matrix(rgamma(n = J*M,
                         shape = rep(a_alpha/b_alpha*rep(1/J,J),
                                     M),
                         rate = 1),
                  nrow = J,
                  ncol = M)
  
  # q_star
  q_star <- matrix(1/R,
                   nrow = J,
                   ncol = R)
  
  # gamma_star
  gamma_star <- rep(a_gamma/b_gamma,
                    J)
  
  # alpha
  alpha <- a_alpha/b_alpha
  
  # alpha0
  alpha_zero <- a0/b0
  
  # w
  w <- rep(1/J, J)
  
  # xi
  xi <- lapply(1:M,
               function(m) rep(1,C[m]))
  
  # lambda
  lambda <- rep(0.5,J)
  
  # tau
  tau <- 5
  
  #----------------------- Step 3: Prepare for the covariance update -----------------------------
  
  # w
  mean_x_w <- matrix(log(w[1:(J-1)]/w[J]),
                     nrow = 1)
  tilde_s_w <- t(mean_x_w)%*%mean_x_w
  covariance_w <- matrix(0,
                         nrow = J-1,
                         ncol = J-1)
  
  # alpha
  X_mean_alpha <- log(alpha)
  M_2_alpha <- 0
  variance_alpha <- 0
  
  # alpha0
  X_mean_alpha0 <- log(alpha_zero)
  M_2_alpha_zero <- 0
  variance_alpha_zero <- 0
  
  # q star
  mean_x_q_star <- lapply(1:J,
                          function(j) matrix(log(q_star[j,1:(R-1)]/q_star[j,R]),
                                             nrow = 1))
  tilde_s_q_star <- lapply(1:J,
                           function(j) t(mean_x_q_star[[j]])%*%mean_x_q_star[[j]])
  
  covariance_q_star <- lapply(1:J,
                              function(j) matrix(0,
                                                 nrow = R-1,
                                                 ncol = R-1))
                                  
  
  # gamma star
  X_mean_gamma <- log(gamma_star)
  M_2_gamma_new <- rep(0,J)
  variance_gamma <- rep(0,J)
  
  
  # alpha r
  mean_X_alpha_r <- matrix(log(rep(a_r/b_r, R)),
                           nrow = 1)
  tilde_s_alpha_r <- t(mean_X_alpha_r)%*%mean_X_alpha_r
  covariance_alpha_r <- matrix(0, nrow = R, ncol = R)
  
  
  
  # tau
  mean_tau <- tau
  X_mean_tau <- log(tau)
  M_2_tau <- 0
  variance_tau <- 0
  
  
  #----------------------- Step 4: Updates -----------------------------
  
  ##-- Set initial index
  output_index <- 0
  
  
  # Iteration starts with number 2
  for(iter in 2:number_iter){
    
    if(iter >= burn_in & iter%%thinning == 0){
      
      output_index <- output_index + 1
      update <- TRUE
    }else{
      
      update <- FALSE
    }
    
    if(iter %% iter_update == 0) print(paste("iter:", iter))
    
    
    ##----------------------------------- Allocation variables ------------------------------------
    
    z_output <- allocation_variables_dirmult_mcmc2(beta = beta,
                                                   delta = delta,
                                                   q_star = q_star,
                                                   gamma_star = gamma_star,
                                                   Y = Y,
                                                   x = x,
                                                   num.cores = num.cores)
    
    Z <- z_output$Z
    allocation.prob <- z_output$allocation.prob
    
    ##-- Print table(Z) if TRUE
    if(isTRUE(print_Z)){
      for(m in 1:M) print(table(Z[[m]]))
    }
    
    
    #--------------------------------------- Delta -------------------------------------------------
    
    delta <- delta_mcmc(alpha = alpha,
                        w = w,
                        z = Z,
                        xi = xi,
                        beta = beta,
                        x = x)
                            
    
    
    #----------------------------------------- Beta ------------------------------------------------
    
    beta <- beta_mcmc(beta = beta,
                      delta = delta,
                      xi = xi,
                      x = x,
                      Z = Z,
                      tau = tau,
                      lambda = lambda,
                      num.cores = num.cores)
                             
    
    #---------------------------------------- tau ---------------------------------------------------
    
    tau_mcmc <- tau_mcmc(lambda = lambda,
                         tau = tau,
                         beta = beta,
                         s = s,
                         variance = variance_tau,
                         M_2 = M_2_tau,
                         X_mean = X_mean_tau,
                         adaptive_prop = adaptive_prop)
    
    tau <- tau_mcmc$tau
    X_mean_tau <- tau_mcmc$X_mean
    M_2_tau <- tau_mcmc$M_2
    variance_tau <- tau_mcmc$variance
    
    tau_count <- tau_count + tau_output$accept
    acceptance_prob_list$tau[iter-1] <- tau_count/(iter-1)
    
    #----------------------------------------- xi ----------------------------------------------------
    
    xi <- xi_mcmc(delta = delta,
                  beta = beta,
                  x = x)
    
    
    #---------------------------------------- w ------------------------------------------------------
    
    w_output <- w_mcmc(w = w,
                       delta = delta,
                       alpha_zero = alpha_zero,
                       alpha = alpha,
                       covariance = covariance_w,
                       mean_x = mean_x_w,
                       tilde_s = tilde_s_w,
                       iter_num = iter_num,
                       adaptive_prop = adaptive_prop)
    
    w <- w_output$w
    tilde_s_w <- w_output$tilde_s
    mean_x_w <- w_output$mean_x
    covariance_w <- w_output$covariance
    
    w_count <- w_count + w_output$accept
    acceptance_prob_list$w[iter-1] <- w_count/(iter-1)
    
    #---------------------------------------- alpha --------------------------------------------------
    
    alpha_output <- alpha_mcmc2(w = w,
                                delta = delta,
                                a_alpha = a_alpha,
                                b_alpha = b_alpha,
                                alpha = alpha,
                                X_mean = X_mean_alpha,
                                M_2 = M_2_alpha,
                                variance = variance_alpha,
                                iter_num = iter_num,
                                adaptive_prop = adaptive_prop)
    
    alpha <- alpha_output$alpha
    X_mean_alpha <- alpha_output$X_mean
    M_2_alpha <- alpha_output$M_2
    variance_alpha <- alpha_output$variance
    
    alpha_count <- alpha_count + alpha_output$accept
    acceptance_prob_list$alpha[iter-1] <- alpha_count/(iter-1)
    
    
    
    #---------------------------------------- alpha_0 -------------------------------------------------
    
    alpha_zero_output <- alpha_zero_mcmc2(w = w,
                                          alpha_zero = alpha_zero,
                                          a0 = a0,
                                          b0 = b0,
                                          X_mean = X_mean_alpha_zero,
                                          M_2 = M_2_alpha_zero,
                                          variance = variance_alpha_zero,
                                          iter_num = iter_num,
                                          adaptive_prop = adaptive_prop)
                                       
    
    alpha_zero <- alpha_zero_output$alpha_zero
    X_mean_alpha_zero <- alpha_zero_output$X_mean
    M_2_alpha_zero <- alpha_zero_output$M_2
    variance_alpha_zero <- alpha_zero_output$variance
    
    alpha_zero_count <- alpha_zero_count + alpha_zero_output$accept
    acceptance_prob_list$alpha_zero[iter-1] <- alpha_zero_count/(iter-1)
    
    
    #----------------------------------- alpha_r -------------------------------------------------------
    
    alpha_r_output <- alpha_r_mcmc(alpha_r = alpha_r,
                                   q_star = q_star,
                                   iter_num = iter_num,
                                   a_r = a_r,
                                   b_r = b_r,
                                   mean_X = mean_X_alpha_r,
                                   tilde_s = tilde_s_alpha_r,
                                   covariance = covariance_alpha_r,
                                   adaptive_prop = adaptive_prop)
    
    alpha_r <- alpha_r_output$alpha_r
    mean_X_alpha_r <- alpha_r_output$mean_X
    tilde_s_alpha_r <- alpha_r_output$tilde_s
    covariance_alpha_r <- alpha_r_output$covariance
    
    alpha_r_count <- alpha_r_count + alpha_r_output$accept
    acceptance_prob_list$alpha_r[iter-1] <- alpha_r_count/(iter-1)
    
    
    #--------------------------------------- q -------------------------------------------------------
    
    q_star_output <- q_star_mcmc2(Y = Y,
                                  Z = Z,
                                  q_star = q_star,
                                  gamma_star = gamma_star,
                                  alpha_r = alpha_r,
                                  covariance = covariance_q_star,
                                  mean_x = mean_x_q_star,
                                  tilde_s = tilde_s_q_star,
                                  iter_num = iter_num,
                                  adaptive_prop = adaptive_prop,
                                  num.cores = num.cores)
    
    q_star <- q_star_output$q_star
    mean_x_q_star <- q_star_output$mean_x
    tilde_s_q_star <- q_star_output$tilde_s
    covariance_q_star <- q_star_output$covariance
    
    # Cumulative acceptance probability
    q_star_count <- q_star_count + q_star_output$q_star_count
    acceptance_prob_list$q_star[iter-1] <- q_star_count/((iter-1)*J)
    
    
    #------------------------------------ gamma ---------------------------------------------------------
    
    gamma_star_output <- gamma_mcmc2(Y = Y,
                                     Z = Z,
                                     q_star = q_star,
                                     gamma_star = gamma_star,
                                     a_gamma = a_gamma,
                                     b_gamma = b_gamma,
                                     X_mean = X_mean_gamma,
                                     M_2 = M_2_gamma,
                                     variance = variance_gamma,
                                     iter_num = iter_num,
                                     adaptive_prop = adaptive_prop,
                                     num.cores = num.cores)
    
    gamma_star <- gamma_star_output$gamma_star
    variance_gamma <- gamma_star_output$variance
    M_2_gamma <- gamma_star_output$M_2
    X_mean_gamma <- gamma_star_output$X_mean

    gamma_star_count <- gamma_star_count + gamma_star_output$gamma_count
    acceptance_prob_list$gamma_star[iter-1] <- gamma_star_count/((iter-1)*J)
    
    
    #------------------------------- w_m -------------------------------
    
    # For x = 0, for region LEC
    w_m_0 <- matrix(0, nrow = J, ncol = M)
    for(m in 1:M){
      
      w_m_0[,m] <- delta[,m]/sum(delta[,m])
    }
    
    # For x = 1, for region MEC
    w_m_1 <- matrix(0, nrow = J, ncol = M)
    for(m in 1:M){
      
      w_m_1[,m] <- delta[,m]*exp(beta)/sum(delta[,m]*exp(beta))
    }
    
    
    #-------------------------- Step 5: Update simulated values ------------------------
    if(update == TRUE){
      
      Z_output[[output_index]] <- Z
      w_output[[output_index]] <- w
      alpha_output[output_index] <- alpha
      alpha_zero_output[output_index] <- alpha_zero
      q_star_output[[output_index]] <- q_star
      gamma_star_output[[output_index]] <- gamma_star
      alpha_r_ouput[[output_index]] <- alpha_r
      allocation_prob_output[[output_index]] <- allocation.prob
      xi_output[[output_index]] <- xi
      delta_output[[output_index]] <- delta
      beta_output[[output_index]] <- beta
      w_m_0_output[[output_index]] <- w_m_0
      w_m_1_output[[output_index]] <- w_m_1
      
    }
    
  }
  
  #------------------------------- Step 6: Return a list of items ----------------------
  
  
  my_list <- list('acceptance_prob' = acceptance_prob_list,
                  'M' = M,
                  'C' = C,
                  'R' = R,
                  'allocation_probability' = allocation_prob_output,
                  
                  'Z_output' = Z_output,
                  'w_output' = w_output,
                  'alpha_output' = alpha_output,
                  'alpha_zero_output' = alpha_zero_output,
                  'q_star_output' = q_star_output,
                  'gamma_star_output' = gamma_star_output,
                  'alpha_r_ouput' = alpha_r_ouput,
                  'xi_output' = xi_output,
                  'delta_output' = delta_output,
                  'beta_output' = beta_output,
                  'w_m_0_output' = w_m_0_output,
                  'w_m_1_output' = w_m_1_output
                  
                  
                  
                  
                  
  )
  
  
  
}