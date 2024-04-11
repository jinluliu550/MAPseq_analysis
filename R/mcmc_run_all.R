#' MCMC Simulation
#'

mcmc_run_all <- function(Y,
                         J,
                         number_iter = 6000,
                         thinning = 5,
                         burn_in = 3000,
                         adaptive_prop = 0.001,
                         print_Z = FALSE,
                         iter_update = 100,
                         a_gamma,
                         b_gamma,
                         a_alpha = 1,
                         b_alpha = 1,
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
  omega_J_M_output <- NULL
  omega_output <- NULL
  alpha_output <- 0
  alpha_zero_output <- 0
  allocation_prob_output <- NULL
  q_star_1_J_output <- NULL
  gamma_star_1_J_output <- NULL
  alpha_h_output <- NULL


  #--------------------------------------------------------------------------------

  # Acceptance probability
  acceptance_prob_list <- data.frame(omega = rep(0, number_iter-1),
                                     alpha = rep(0, number_iter-1),
                                     alpha_zero = rep(0, number_iter-1),
                                     q_star = rep(0, number_iter-1),
                                     gamma_star = rep(0, number_iter-1),
                                     alpha_h = rep(0, number_iter-1))


  #----------------------- Step 2: Set starting values ----------------------------

  Z_new <- lapply(1:M,
                  function(m) sample(1:J, size = C[m], replace = TRUE))
  omega_J_M_new <- matrix(1/J,
                          nrow = J,
                          ncol = M)
  omega_new <- rep(1/J,
                   J)
  
  alpha_new <- 1
  alpha_zero_new <- 1

  # Based on Z_new
  q_star_1_J_new <- lapply(1:J,
                          function(j){
                            
                            if(length(which(unlist(Z_new)==j)) != 0){
                              
                              data.j.prior <- do.call(cbind, Y)[,which(unlist(Z_new)==j)]
                              
                              data.j.mean <- rowMeans(data.j.prior)
                              
                              mx <- matrix(data.j.mean/sum(data.j.mean),
                                           nrow = 1)
                              
                            }else{
                              
                              mx <- matrix(rep(1:R, R), nrow = 1)
                            }
                            
                            mx
                            
                          })
  
  q_star_1_J_new <- do.call(rbind, q_star_1_J_new)
  


  # Let gamma be a large number
  gamma_1_J_star_new <- rep(a_gamma/b_gamma,
                            J)

  alpha_h_new <- rep(a_alpha/b_alpha,
                     R)


  # Total count of acceptance
  omega_count <- 0
  alpha_count <- 0
  alpha_zero_count <- 0
  q_star_count <- 0
  gamma_star_count <- 0
  alpha_h_count <- 0

  #----------------------- Step 3: Prepare for the covariance update -----------------------------

  # omega
  mean_X_component_new <- matrix(log(omega_new[1:J-1]/omega_new[J]),
                                 nrow = 1)
  tilde_s_component_new <- t(mean_X_component_new)%*%mean_X_component_new
  covariance_component_new <- matrix(0,
                                     nrow = J-1,
                                     ncol = J-1)

  # alpha
  mean_X_alpha_new <- log(alpha_new)
  M_2_alpha_new <- 0
  variance_alpha_new <- 0

  # alpha_zero
  mean_X_alpha_zero_new <- log(alpha_zero_new)
  M_2_alpha_zero_new <- 0
  variance_alpha_zero_new <- 0

  # q star
  mean_x_q_new <- lapply(1:J, function(j) matrix(log(q_star_1_J_new[j,1:(R-1)]/q_star_1_J_new[j,R]),
                                                 nrow = 1))
  tilde_s_q_new <- lapply(1:J, function(j) t(mean_x_q_new[[j]])%*%mean_x_q_new[[j]])
  covariance_q_new <- lapply(1:J, function(j) matrix(0, nrow = R-1, ncol = R-1))

  # gamma star
  X_mean_gamma_new <- log(gamma_1_J_star_new)
  M_2_gamma_new <- rep(0,J)
  variance_gamma_new <- rep(0,J)

  # alpha_h
  mean_X_alpha_h_new <- matrix(log(alpha_h_new),
                               nrow = 1)
  tilde_s_alpha_h_new <- t(mean_X_alpha_h_new)%*%mean_X_alpha_h_new
  covariance_alpha_h_new <- matrix(0, nrow = R, ncol = R)


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

    allocation_output <- allocation_variables_dirmult_mcmc(omega_J_M = omega_J_M_new,
                                                           q_star_1_J = q_star_1_J_new,
                                                           gamma_1_J_star = gamma_1_J_star_new,
                                                           Y = Y,
                                                           num.cores = num.cores)

    Z_new <- allocation_output$Z
    allocation.prob <- allocation_output$allocation.prob

    ##-- Print table(Z) if TRUE
    if(isTRUE(print_Z)){
      for(m in 1:M) print(table(Z_new[[m]]))
    }
    
    print(paste0('J' = length(unique(unlist(Z_new)))))

    ##------------------------ Dataset specific component probabilities ------------------------------

    dataset_specfic_output <- dataset_specific_mcmc(Z = Z_new,
                                                    omega = omega_new,
                                                    alpha = alpha_new)

    omega_J_M_new <- dataset_specfic_output

    ##------------------------ Component probabilities -----------------------------------------------

    component_output <- component_probabilities_mcmc(omega = omega_new,
                                                     omega_J_M = omega_J_M_new,
                                                     alpha_zero = alpha_zero_new,
                                                     alpha = alpha_new,
                                                     covariance = covariance_component_new,
                                                     mean_x = mean_X_component_new,
                                                     tilde_s = tilde_s_component_new,
                                                     iter_num = iter,
                                                     adaptive_prop = adaptive_prop)

    omega_new <- component_output$omega_new
    tilde_s_component_new <- component_output$tilde_s_new
    mean_X_component_new <- component_output$mean_x_new
    covariance_component_new <- component_output$covariance_new
    omega_count <- omega_count + component_output$accept
    acceptance_prob_list$omega[iter-1] <- omega_count/(iter-1)

    ##-------------------------- Update alpha ----------------------------------

    alpha_output_sim <- alpha_mcmc(omega_J_M = omega_J_M_new,
                                   omega = omega_new,
                                   alpha = alpha_new,
                                   X_mean = mean_X_alpha_new,
                                   M_2 = M_2_alpha_new,
                                   variance = variance_alpha_new,
                                   iter_num = iter,
                                   adaptive_prop = adaptive_prop)

    alpha_new <- alpha_output_sim$alpha_new
    mean_X_alpha_new <- alpha_output_sim$X_mean_new
    M_2_alpha_new <- alpha_output_sim$M_2_new
    variance_alpha_new <- alpha_output_sim$variance_new
    alpha_count <- alpha_count + alpha_output_sim$accept
    acceptance_prob_list$alpha[iter-1] <- alpha_count/(iter-1)

    ##----------------------- Update alpha_zero --------------------------------

    alpha_zero_output_sim <- alpha_zero_mcmc(omega = omega_new,
                                             alpha_zero = alpha_zero_new,
                                             X_mean = mean_X_alpha_zero_new,
                                             M_2 = M_2_alpha_zero_new,
                                             variance = variance_alpha_zero_new,
                                             iter_num = iter,
                                             adaptive_prop = adaptive_prop)

    alpha_zero_new <- alpha_zero_output_sim$alpha_zero_new
    mean_X_alpha_zero_new <- alpha_zero_output_sim$X_mean_new
    M_2_alpha_zero_new <- alpha_zero_output_sim$M_2_new
    variance_alpha_zero_new <- alpha_zero_output_sim$variance_new
    alpha_zero_count <- alpha_zero_count + alpha_zero_output_sim$accept
    acceptance_prob_list$alpha_zero[iter-1] <- alpha_zero_count/(iter-1)

    ##--------------------------- alpha_h ---------------------------------

    alpha_h_output_sim <- alpha_h_mcmc(alpha_h = alpha_h_new,
                                       mean_X_alpha_h = mean_X_alpha_h_new,
                                       tilde_s_alpha_h = tilde_s_alpha_h_new,
                                       q_star_1_J = q_star_1_J_new,
                                       iter_num = iter,
                                       a_alpha = a_alpha,
                                       b_alpha = b_alpha,
                                       covariance = covariance_alpha_h_new,
                                       adaptive_prop = adaptive_prop)

    alpha_h_new <- alpha_h_output_sim$alpha_h_new
    tilde_s_alpha_h_new <- alpha_h_output_sim$tilde_s_alpha_h_new
    mean_X_alpha_h_new <- alpha_h_output_sim$mean_X_alpha_h_new
    covariance_alpha_h_new <- alpha_h_output_sim$covariance_alpha_h_new
    alpha_h_count <- alpha_h_count + alpha_h_output_sim$accept
    acceptance_prob_list$alpha_h[iter-1] <- alpha_h_count/(iter-1)

    ##--------------------------- q_star ----------------------------------

    q_output_sim <- q_star_mcmc(Y = Y,
                                Z = Z_new,
                                q_star_1_J = q_star_1_J_new,
                                gamma_1_J_star = gamma_1_J_star_new,
                                alpha_h = alpha_h_new,
                                covariance = covariance_q_new,
                                mean_x = mean_x_q_new,
                                tilde_s = tilde_s_q_new,
                                iter_num = iter,
                                adaptive_prop = adaptive_prop)

    covariance_q_new <- q_output_sim$covariance_new
    mean_x_q_new <- q_output_sim$mean_x_new
    tilde_s_q_new <- q_output_sim$tilde_s_new
    q_star_1_J_new <- q_output_sim$q_star_1_J_new
    q_star_count <- q_output_sim$q_star_count + q_star_count
    acceptance_prob_list$q_star[iter-1] <- q_star_count/((iter-1)*J)


    ##-------------------------- gamma ---------------------------------

    gamma_output_sim <- gamma_mcmc(Y = Y,
                                   Z = Z_new,
                                   q_star_1_J = q_star_1_J_new,
                                   gamma_1_J_star = gamma_1_J_star_new,
                                   a_gamma = a_gamma,
                                   b_gamma = b_gamma,
                                   X_mean = X_mean_gamma_new,
                                   M_2 = M_2_gamma_new,
                                   variance = variance_gamma_new,
                                   iter_num = iter,
                                   adaptive_prop = adaptive_prop)

    gamma_1_J_star_new <- gamma_output_sim$gamma_1_J_star_new
    variance_gamma_new <- gamma_output_sim$variance_gamma_new
    M_2_gamma_new <- gamma_output_sim$M_2_gamma_new
    X_mean_gamma_new <- gamma_output_sim$X_mean_gamma_new
    gamma_star_count <- gamma_star_count + gamma_output_sim$gamma_count
    acceptance_prob_list$gamma_star[iter-1] <- gamma_star_count/((iter-1)*J)
    
    print(paste0('max gamma = ', max(gamma_1_J_star_new), 'min gamma = ', min(gamma_1_J_star_new)))


    #-------------------------- Step 5: Update simulated values ------------------------
    if(update == TRUE){

      Z_output[[output_index]] <- Z_new
      omega_J_M_output[[output_index]] <- omega_J_M_new
      omega_output[[output_index]] <- omega_new
      alpha_output[output_index] <- alpha_new
      alpha_zero_output[output_index] <- alpha_zero_new
      q_star_1_J_output[[output_index]] <- q_star_1_J_new
      allocation_prob_output[[output_index]] <- allocation.prob
      alpha_h_output[[output_index]] <- alpha_h_new
      gamma_star_1_J_output[[output_index]] <- gamma_1_J_star_new
    }

  }


  my_list <- list('acceptance_prob' = acceptance_prob_list,
                  'M' = M,
                  'C' = C,
                  'R' = R,
                  'allocation_probability' = allocation_prob_output,
                  'Z_output' = Z_output,
                  'omega_J_M_output' = omega_J_M_output,
                  'omega_output' = omega_output,
                  'alpha_output' = alpha_output,
                  'alpha_zero_output' = alpha_zero_output,
                  'alpha_h_output' = alpha_h_output,
                  'q_star_1_J_output' = q_star_1_J_output,
                  'gamma_star_1_J_output' = gamma_star_1_J_output




  )

  return(my_list)

}



