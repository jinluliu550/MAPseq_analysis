
#' MCMC sampling of gamma and projection strength
#'

mcmc_run_post <- function(mcmc_run_all_output,
                          Z,
                          thinning = 5,
                          burn_in,
                          number_iter,
                          Y,
                          regions.name = NULL,
                          a_gamma,
                          b_gamma,
                          adaptive_prop = 0.1,
                          axis.title.size = 10,
                          axis.text.size = 10,
                          plot.title.size = 10,
                          cred.int = 0.9){
  
  
  
  # Dimensions
  J <- length(unique(unlist(Z)))
  R <- mcmc_run_all_output$R
  M <- length(Y)
  C <- sapply(1:M,
              function(m) ncol(Y[[m]]))
  
  if(is.null(regions.name)){
    
    regions.name <- paste('region', 1:R)
  }
  
  # Posterior mean of alpha h
  alpha_h_all <- do.call(rbind,
                         mcmc_run_all_output$alpha_h_output)
  mean_alpha_h <- colSums(alpha_h_all)/nrow(alpha_h_all)
  
  # Initial value of q star
  q_star_1_J_new <- matrix(1/R,
                           nrow = J,
                           ncol = R)
  
  # Initial covariance structure
  mean_x_q_new <- lapply(1:J, function(j) matrix(log(q_star_1_J_new[j,1:(R-1)]/q_star_1_J_new[j,R]),
                                                 nrow = 1))
  tilde_s_q_new <- lapply(1:J, function(j) t(mean_x_q_new[[j]])%*%mean_x_q_new[[j]])
  covariance_q_new <- lapply(1:J, function(j) matrix(0, nrow = R, ncol = R))
  
  
  # Let gamma start with a large number
  gamma_1_J_star_new <- rep(a_gamma/b_gamma, J)
  
  # gamma star covariance structure
  X_mean_gamma_new <- log(gamma_1_J_star_new)
  M_2_gamma_new <- rep(0,J)
  variance_gamma_new <- rep(0,J)
  
  
  # projection probabilities
  q_star_1_J_output <- NULL
  
  # gamma star
  gamma_star_1_J_output <- NULL
  
  # Set initial index
  output_index <- 0
  
  # Iteration starts with number 2
  for(iter in 2:number_iter){
    
    if(iter >= burn_in & iter%%thinning == 0){
      
      output_index <- output_index + 1
      update <- TRUE
    }else{
      
      update <- FALSE
    }
    
    ##-- Print for every 100th iteration
    if(iter %% 100 == 0) print(paste("iter:", iter))
    
    # Posterior sample of gamma
    gamma_output_sim <- gamma_mcmc(Y = Y,
                                   Z = Z,
                                   q_star_1_J = q_star_1_J_new,
                                   gamma_1_J_star = gamma_1_J_star_new,
                                   a_gamma = a_gamma,
                                   b_gamma = b_gamma,
                                   X_mean = X_mean_gamma_new,
                                   M_2 = M_2_gamma_new,
                                   variance = variance_gamma_new,
                                   iter_num = iter,
                                   adaptive_prop = adaptive_prop)
    
    # Update covariance structure of gamma
    gamma_1_J_star_new <- gamma_output_sim$gamma_1_J_star_new
    variance_gamma_new <- gamma_output_sim$variance_gamma_new
    M_2_gamma_new <- gamma_output_sim$M_2_gamma_new
    X_mean_gamma_new <- gamma_output_sim$X_mean_gamma_new
    
    # Posterior sample of projection probability
    q_output_sim <- q_star_mcmc(Y = Y,
                                Z = Z,
                                q_star_1_J = q_star_1_J_new,
                                gamma_1_J_star = gamma_1_J_star_new,
                                alpha_h = mean_alpha_h,
                                covariance = covariance_q_new,
                                mean_x = mean_x_q_new,
                                tilde_s = tilde_s_q_new,
                                iter_num = iter,
                                adaptive_prop = adaptive_prop,
                                num.cores = num.cores)
    
    # Update covariance sturcture of q
    covariance_q_new <- q_output_sim$covariance_new
    mean_x_q_new <- q_output_sim$mean_x_new
    tilde_s_q_new <- q_output_sim$tilde_s_new
    q_star_1_J_new <- q_output_sim$q_star_1_J_new
    
    # If update is TRUE
    if(isTRUE(update)){
      
      q_star_1_J_output[[output_index]] <- q_star_1_J_new
      gamma_star_1_J_output[[output_index]] <- gamma_1_J_star_new
    }
    
  }
  
  ##-- Compute posterior mean
  proj_prob_mean <- Reduce("+", q_star_1_J_output)/length(q_star_1_J_output)
  gamma_mean <- Reduce("+", gamma_star_1_J_output)/length(gamma_star_1_J_output)
  
  ##-- Credible intervals
  proj_prob_lower <- matrix(0, nrow = J, ncol = R)
  proj_prob_upper <- matrix(0, nrow = J, ncol = R)
  
  for(j in 1:J){
    for(r in 1:R){
      
      proj_prob_lower[j,r] <- quantile(sapply(1:length(q_star_1_J_output), function(i) q_star_1_J_output[[i]][j,r]), (1-cred.int)/2)
      proj_prob_upper[j,r] <- quantile(sapply(1:length(q_star_1_J_output), function(i) q_star_1_J_output[[i]][j,r]), 1-(1-cred.int)/2)
      
    }
  }
  
  
  
  #----------------------- Label clusters
  
  # How many areas are each cluster projecting
  neuron_projection_df <- list('q_star_1_J_output' = q_star_1_J_output,
                               'gamma_star_1_J_output' = gamma_star_1_J_output,
                               'Z' = Z,
                               'J' = J)
  
  
  q_tilde_001 <- q_tilde(neuron_projection_df = neuron_projection_df,
                         epsilon = 0.01)
  
  
  unicast.index <- (1:J)[apply(q_tilde_001,
                               1,
                               function(x) length(which(x>0.5)) == 1)]
  
  
  broadcast.index <- (1:J)[apply(q_tilde_001,
                                 1,
                                 function(x) length(which(x>0.5)) == R)]
  
  bicast.and.more.index <- (1:J)[-unique(c(unicast.index,
                                           broadcast.index))]
  
  # Summary table of clusters
  cluster.label.summary <- data.frame(cluster = 1:J,
                                      
                                      unicast = ifelse((1:J) %in% unicast.index, 'Yes', 'No'),
                                      bicast.and.more = ifelse((1:J) %in% bicast.and.more.index, 'Yes', 'No'),
                                      broadcast = ifelse((1:J) %in% broadcast.index, 'Yes', 'No')
  )
  
  # Relabel
  strong.proj <- apply(proj_prob_mean,
                       1,
                       function(x) which(x == max(x))[1])
  
  
  df0 <- lapply(1:ncol(proj_prob_mean),
                function(i){
                  
                  # Clusters with the strongest projection to region i
                  strong.proj.i <- which(strong.proj==i)
                  
                  if(length(strong.proj.i) != 0){
                    
                    qj_estimate_i <- matrix(proj_prob_mean[strong.proj.i,],
                                            nrow = length(strong.proj.i))
                    
                    output <- strong.proj.i[order(-qj_estimate_i[,i])]
                  }else{
                    
                    output <- NULL
                  }
                  
                  
                }
  )
  
  original.label <- unlist(df0)
  
  
  
  
  
  # Relabel the summary table
  cluster.label.summary <- cluster.label.summary[original.label,]
  cluster.label.summary$cluster <- 1:J
  
  
  # Change to the new labeling system
  q_star_1_J_output <- lapply(1:length(q_star_1_J_output),
                              function(t) q_star_1_J_output[[t]][original.label,])
  
  gamma_star_1_J_output <- lapply(1:length(gamma_star_1_J_output),
                                  function(t) gamma_star_1_J_output[[t]][original.label])
  
  proj_prob_mean <- proj_prob_mean[original.label,]
  proj_prob_lower <- proj_prob_lower[original.label,]
  proj_prob_upper <- proj_prob_upper[original.label,]
  
  gamma_mean <- gamma_mean[original.label]
  
  
  
  Z <- lapply(1:M,
              function(m){
                
                sapply(1:C[m], function(c) which(original.label == Z[[m]][c]))
              })
  
  
  
  # Summary table of projection probabilities
  loop.result <- lapply(1:J, function(j){
    
    data.frame(cluster = paste('cluster', j),
               region = regions.name,
               projection.mean = proj_prob_mean[j,],
               projection.lower = proj_prob_lower[j,],
               projection.upper = proj_prob_upper[j,])
  })
  
  estimated.projection.df <- do.call(rbind, loop.result)
  
  
  #---------------------------------------- Plot of projection strengths ------------------------------------
  
  
  # Label the number of projecting regions
  cluster_label_vec <- sapply(1:J,
                              function(j){
                                
                                colnames(cluster.label.summary)[2:4][which(cluster.label.summary[j,2:4]=='Yes')[1]]
                              })
  
  estimated.projection.df$class <- rep(cluster_label_vec, each = R)
  estimated.projection.df$class <- factor(estimated.projection.df$class,
                                          levels = colnames(cluster.label.summary)[2:4])
  
  estimated.projection.df$cluster <- factor(estimated.projection.df$cluster,
                                            levels = paste('cluster', 1:length(unique(unlist(Z)))))
  
  
  # All clusters
  estimated.pp.plot <- ggplot2::ggplot(estimated.projection.df,
                                       mapping = aes(x = factor(region, levels = regions.name),
                                                     y = projection.mean,
                                                     group = cluster,
                                                     color = class)) +
    geom_line()+
    geom_point()+
    geom_errorbar(aes(ymin = projection.lower,
                      ymax = projection.upper),
                  width = 0.1)+
    theme_bw()+
    ylim(c(0,1))+
    ylab('projection strength')+
    xlab('region')+
    facet_wrap(~cluster)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(axis.text=element_text(size=axis.text.size),
          axis.title=element_text(size=axis.title.size),
          plot.title = element_text(size=plot.title.size))
  
  
  
  q_tilde_001 <- q_tilde_001[original.label,]
  
  colnames(q_tilde_001) <- regions.name
  rownames(q_tilde_001) <- paste('cluster', 1:J)
  
  
  
  q_tilde_001_table <- data.frame(region = rep(regions.name,
                                               each = J),
                                  cluster = rep(paste('cluster', 1:J),
                                                R),
                                  probability = as.vector(q_tilde_001))
  
  q_tilde_001_table$region <- factor(q_tilde_001_table$region,
                                     levels = regions.name)
  
  q_tilde_001_table$cluster <- factor(q_tilde_001_table$cluster,
                                      levels = paste('cluster', 1:J))
  
  q_tilde_plot <- q_tilde_001_table %>%
    ggplot(mapping = aes(x = cluster,
                         y = region,
                         fill = probability))+
    geom_tile()+
    theme_bw()+
    scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                         values = scales::rescale(c(0, 0.25, 1)),
                         limits = c(0,1))+
    xlab('cluster')+
    ylab('regions')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  ##-- Output values
  return(list('J' = J,
              'q_star_1_J_output' = q_star_1_J_output,
              'proj_prob_mean' = proj_prob_mean,
              'proj_prob_lower' = proj_prob_lower,
              'proj_prob_upper' = proj_prob_upper,
              'Z' = Z,
              'gamma_mean' = gamma_mean,
              'gamma_star_1_J_output' = gamma_star_1_J_output,
              'estimated.projection.df' = estimated.projection.df,
              'cluster.label.summary' = cluster.label.summary,
              'estimated.pp.plot' = estimated.pp.plot,
              'q_tilde_001' = q_tilde_001,
              'q_tilde_plot' = q_tilde_plot)
  )
}

