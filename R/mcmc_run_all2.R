
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
  
  #------------------------ List of acceptance probabilities --------------------------
  
  acceptance_prob_list <- data.frame(w = rep(0, number_iter-1),
                                     alpha_0 = rep(0, number_iter-1),
                                     alpha = rep(0, number_iter-1),
                                     beta = rep(0, number_iter-1),
                                     q_star = rep(0, number_iter-1),
                                     gamma_star = rep(0, number_iter-1),
                                     alpha_r = rep(0, number_iter-1))
  
  beta_count <- 0
  alpha_count <- 0
  alpha_0_count <- 0
  
  
  
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
    
    alpha_0_output <- alpha_zero_mcmc2(w = w_new,
                                       alpha_zero = alpha_zero_new,
                                       a0 = a0,
                                       b0 = b0,
                                       X_mean = X_mean_alpha0_new,
                                       M_2 = M_2_alpha0_new,
                                       variance = variance_alpha0_new,
                                       iter_num = iter_num,
                                       adaptive_prop = adaptive_prop)
    
    alpha_zero_new <- alpha_0_output$alpha_zero_new
    X_mean_alpha0_new <- alpha_0_output$X_mean_new
    
    
  }
  
  
}