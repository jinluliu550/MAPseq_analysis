
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
                          num.cores = 1,
                          pi_cons,
                          s1_2,
                          s2_2){
  
  
  
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
                                     beta = rep(0, number_iter-1),
                                     q_star = rep(0, number_iter-1),
                                     gamma_star = rep(0, number_iter-1),
                                     alpha_r = rep(0, number_iter-1))
  
  beta_count <- 0
  alpha_count <- 0
  alpha_zero_count <- 0
  w_count <- 0
  alpha_r_count <- 0
  q_star_count <- 0
  gamma_star_count <- 0
  
  
  #----------------------- Step 2: Set starting values ----------------------------
  
  # beta
  beta_new <- rep(1,J)
  
  # delta
  delta_new <- matrix(rgamma(n = J*M,
                             shape = rep(a_alpha/b_alpha*rep(1/J,J),
                                         M),
                             rate = 1),
                      nrow = J,
                      ncol = M)
  
  # q_star
  q_star_new <- matrix(1/R,
                       nrow = J,
                       ncol = R)
  
  # gamma_star
  gamma_star <- rep(a_gamma/b_gamma,
                    J)
  
  # alpha
  alpha_new <- a_alpha/b_alpha
  
  # alpha0
  alpha_zero_new <- a0/b0
  
  # w
  w_new <- rep(1/J,
               J)
  
  # xi
  xi_new <- lapply(1:M,
                   function(m) rep(1,C[m]))
  
  
  
  #----------------------- Step 3: Prepare for the covariance update -----------------------------
  
  # w
  mean_x_w_new <- matrix(log(w_new[1:(J-1)]/w_new[J]),
                         nrow = 1)
  tilde_s_w_new <- t(mean_x_w_new)%*%mean_x_w_new
  covariance_w_new <- matrix(0,
                             nrow = J-1,
                             ncol = J-1)
  
  # alpha
  X_mean_alpha_new <- log(alpha_new)
  M_2_alpha_new <- 0
  variance_alpha_new <- 0
  
  # alpha0
  X_mean_alpha0_new <- log(alpha_zero_new)
  M_2_alpha_zero_new <- 0
  variance_alpha_zero_new <- 0
  
  # q star
  mean_x_q_star_new <- lapply(1:J,
                              function(j) matrix(log(q_star_new[j,1:(R-1)]/q_star_new[j,R]),
                                                 nrow = 1))
  tilde_s_q_star_new <- lapply(1:J,
                               function(j) t(mean_x_q_star_new[[j]])%*%mean_x_q_star_new[[j]])
  
  covariance_q_star_new <- lapply(1:J,
                                  function(j) matrix(0,
                                                     nrow = R-1,
                                                     ncol = R-1))
  
  # gamma star
  X_mean_gamma_new <- log(gamma_star_new)
  M_2_gamma_new <- rep(0,J)
  variance_gamma_new <- rep(0,J)
  
  
  # alpha r
  mean_X_alpha_r_new <- alpha_r_output$mean_X_new
  tilde_s_alpha_r_new <- t(mean_X_alpha_r_new)%*%mean_X_alpha_r_new
  covariance_alpha_r_new <- matrix(0, nrow = R, ncol = R)
  
  
  # beta
  beta_mean_new <- beta_new
  M2_beta_new <- rep(0,J)
  beta_variance_new <- rep(0,J)
  
  
  
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
    
    z_output <- allocation_variables_dirmult_mcmc2(beta = beta_new,
                                                   delta = delta_new,
                                                   q_star = q_star_new,
                                                   gamma_star = gamma_star_new,
                                                   Y = Y,
                                                   x = x,
                                                   num.cores = num.cores)
    
    Z_new <- z_output$Z
    allocation.prob <- z_output$allocation.prob
    
    ##-- Print table(Z) if TRUE
    if(isTRUE(print_Z)){
      for(m in 1:M) print(table(Z_new[[m]]))
    }
    
    
    #--------------------------------------- Delta -------------------------------------------------
    
    delta_new <- delta_mcmc(alpha = alpha_new,
                            w = w_new,
                            z = Z_new,
                            xi = xi_new,
                            beta = beta_new,
                            x = x)
    
    
    #----------------------------------------- Beta ------------------------------------------------
    
    beta_output <- beta_mcmc(beta = beta_new,
                             delta = delta_new,
                             xi = xi_new,
                             x = x,
                             M2 = M2_beta_new,
                             z = Z_new,
                             pi_cons = pi_cons,
                             s1_2 = s1_2,
                             s2_2 = s2_2,
                             beta_mean = beta_mean_new,
                             beta_variance = beta_variance_new,
                             iter_num = iter,
                             adaptive_prop = adaptive_prop,
                             num.cores = num.cores)
    
    beta_new <- beta_output$beta_new
    beta_mean_new <- beta_output$beta_mean_new
    M2_beta_new <- beta_output$M2_new
    beta_variance_new <- beta_output$variance_new
    
    # Cumulative acceptance probability
    beta_count <- beta_count + beta_output$accept
    acceptance_prob_list$beta[iter-1] <- beta_count/((iter-1)*J)
    
    
    #----------------------------------------- xi ----------------------------------------------------
    
    xi_new <- xi_mcmc(delta = delta_new,
                      beta = beta_new,
                      x = x)
    
    
    #---------------------------------------- w ------------------------------------------------------
    
    w_output <- w_mcmc(w = w_new,
                       delta = delta_new,
                       alpha_zero = alpha_zero_new,
                       alpha = alpha_new,
                       covariance = covariance_w_new,
                       mean_x = mean_x_w_new,
                       tilde_s = tilde_s_w_new,
                       iter_num = iter_num,
                       adaptive_prop = adaptive_prop)
    
    w_new <- w_output$w_new
    tilde_s_w_new <- w_output$tilde_s_new
    mean_x_w_new <- w_output$mean_x_new
    covariance_w_new <- w_output$covariance_new
    
    # Cumulative acceptance probability
    w_count <- w_count + w_output$accept
    acceptance_prob_list$w[iter-1] <- w_count/(iter-1)
    
    #---------------------------------------- alpha --------------------------------------------------
    
    alpha_output <- alpha_mcmc2(w = w_new,
                                delta = delta_new,
                                a_alpha = a_alpha,
                                b_alpha = b_alpha,
                                alpha = alpha_new,
                                X_mean = X_mean_alpha_new,
                                M_2 = M_2_alpha_new,
                                variance = variance_alpha_new,
                                iter_num = iter_num,
                                adaptive_prop = adaptive_prop)
    
    alpha_new <- alpha_output$alpha_new
    X_mean_alpha_new <- alpha_output$X_mean_new
    M_2_alpha_new <- alpha_output$M_2_new
    variance_alpha_new <- alpha_output$variance_new
    
    # Cumulative acceptance probability
    alpha_count <- alpha_count + alpha_output$accept
    acceptance_prob_list$alpha[iter-1] <- alpha_count/(iter-1)
    
    
    
    #---------------------------------------- alpha_0 -------------------------------------------------
    
    alpha_zero_output <- alpha_zero_mcmc2(w = w_new,
                                          alpha_zero = alpha_zero_new,
                                          a0 = a0,
                                          b0 = b0,
                                          X_mean = X_mean_alpha_zero_new,
                                          M_2 = M_2_alpha_zero_new,
                                          variance = variance_alpha_zero_new,
                                          iter_num = iter_num,
                                          adaptive_prop = adaptive_prop)
                                       
    
    alpha_zero_new <- alpha_zero_output$alpha_zero_new
    X_mean_alpha0_new <- alpha_zero_output$X_mean_new
    M_2_alpha_zero_new <- alpha_zero_output$M_2_new
    variance_alpha_zero_new <- alpha_zero_output$variance_new
    
    # Cumulative acceptance probability
    alpha_zero_count <- alpha_zero_count + alpha_zero_output$accept
    acceptance_prob_list$alpha_zero[iter-1] <- alpha_zero_count/(iter-1)
    
    
    #----------------------------------- alpha_r -------------------------------------------------------
    
    alpha_r_output <- alpha_r_mcmc(alpha_r = alpha_r_new,
                                   mean_X = mean_X_alpha_r_new,
                                   tilde_s = tilde_s_alpha_r_new,
                                   q_star = q_star_new,
                                   iter_num = iter_num,
                                   a_alpha = a_alpha,
                                   b_alpha = b_alpha,
                                   covariance = covariance_alpha_r_new,
                                   adaptive_prop = adaptive_prop)
    
    alpha_r_new <- alpha_r_output$alpha_r_new
    mean_X_alpha_r_new <- alpha_r_output$mean_X_new
    tilde_s_alpha_r_new <- alpha_r_output$tilde_s_new
    covariance_alpha_r_new <- alpha_r_output$covariance_new
    
    # Cumulative acceptance probability
    alpha_r_count <- alpha_r_count + alpha_r_output$accept
    acceptance_prob_list$alpha_r[iter-1] <- alpha_r_count/(iter-1)
    
    
    #--------------------------------------- q -------------------------------------------------------
    
    q_star_output <- q_star_mcmc2(Y = Y,
                                  Z = Z_new,
                                  q_star = q_star_new,
                                  gamma_star = gamma_star_new,
                                  alpha_r = alpha_r_new,
                                  covariance = covariance_q_star_new,
                                  mean_x = mean_x_q_star_new,
                                  tilde_s = tilde_s_q_star_new,
                                  iter_num = iter_num,
                                  adaptive_prop = adaptive_prop,
                                  num.cores = num.cores)
    
    q_star_new <- q_star_output$q_star_new
    mean_x_q_star_new <- q_star_output$mean_x_new
    tilde_s_q_star_new <- q_star_output$tilde_s_new
    covariance_q_star_new <- q_star_output$covariance_new
    
    # Cumulative acceptance probability
    q_star_count <- q_star_count + q_star_output$q_star_count
    acceptance_prob_list$q_star[iter-1] <- q_star_count/((iter-1)*J)
    
    
    #------------------------------------ gamma ---------------------------------------------------------
    
    gamma_star_output <- gamma_mcmc2(Y = Y,
                                     Z = Z_new,
                                     q_star = q_star_new,
                                     gamma_star = gamma_star_new,
                                     a_gamma = a_gamma,
                                     b_gamma = b_gamma,
                                     X_mean = X_mean_gamma_new,
                                     M_2 = M_2_gamma_new,
                                     variance = variance_gamma_new,
                                     iter_num = iter_num,
                                     adaptive_prop = adaptive_prop,
                                     num.cores = num.cores)
    
    gamma_star_new <- gamma_star_output$gamma_star_new
    variance_gamma_new <- gamma_star_output$variance_new
    M_2_gamma_new <- gamma_star_output$M_2_new
    X_mean_gamma_new <- gamma_star_output$X_mean_new
    
    # Cumulative acceptance probability
    gamma_star_count <- gamma_star_count + gamma_star_output$gamma_count
    acceptance_prob_list$gamma_star[iter-1] <- gamma_star_count/((iter-1)*J)
    
    
    #------------------------------- w_m -------------------------------
    
    # For x = 0, for region LEC
    w_m_0 <- matrix(0, nrow = J, ncol = M)
    for(m in 1:M){
      
      w_m_0[,m] <- delta_new[,m]/sum(delta_new[,m])
    }
    
    # For x = 1, for region MEC
    w_m_1 <- matrix(0, nrow = J, ncol = M)
    for(m in 1:M){
      
      w_m_1[,m] <- delta_new[,m]*exp(beta_new)/sum(delta_new[,m]*exp(beta_new))
    }
    
    
    #-------------------------- Step 5: Update simulated values ------------------------
    if(update == TRUE){
      
      Z_output[[output_index]] <- Z_new
      w_output[[output_index]] <- w_new 
      alpha_output[output_index] <- alpha_new
      alpha_zero_output[output_index] <- alpha_zero_new
      q_star_output[[output_index]] <- q_star_new
      gamma_star_output[[output_index]] <- gamma_star
      alpha_r_ouput[[output_index]] <- alpha_r
      allocation_prob_output[[output_index]] <- allocation.prob
      xi_output[[output_index]] <- xi_new
      delta_output[[output_index]] <- delta_new
      beta_output[[output_index]] <- beta_new
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