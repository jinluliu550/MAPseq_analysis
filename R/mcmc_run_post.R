
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
                          cred.int = 0.9,
                          
                          # Brain regions project to the front, middle and back of the brain
                          front.regions = NULL,
                          middle.regions = NULL,
                          back.regions = NULL){
  
  # If any of the label are present, we carry out relabeling
  if(all(c(is.null(front.regions),
           is.null(middle.regions),
           is.null(back.regions))
  )){
    
    relabel <- FALSE
    
  }else{
    
    relabel <- TRUE
  }
  
  
  
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
                                adaptive_prop = adaptive_prop)
    
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
  neuron_projection_df <- NULL
  neuron_projection_df$q_star_1_J_output <- q_star_1_J_output
  neuron_projection_df$gamma_star_1_J_output <- gamma_star_1_J_output
  neuron_projection_df$Z <- Z
  neuron_projection_df$J <- J
  
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
  
  if(isTRUE(relabel)){
    
    
    # Area of projection
    
    # Only project to the back
    if(!is.null(back.regions)){
      
      
      project.to.back <- (1:J)[apply(q_tilde_001,
                                     1,
                                     function(x) (any(x[back.regions] > 0.5)) & (all(x[-back.regions] < 0.5)))]
    }else{
      
      project.to.back <- NULL
    }
    
    
    # Only project to the front
    if(!is.null(front.regions)){
      
      project.to.front <- (1:J)[apply(q_tilde_001,
                                      1,
                                      function(x) (any(x[front.regions] > 0.5)) & (all(x[-front.regions] < 0.5)))]
      
    }else{
      
      project.to.front <- NULL
    }
    
    
    # Only project to the middle
    if(!is.null(middle.regions)){
      
      project.to.middle <- (1:J)[apply(q_tilde_001,
                                       1,
                                       function(x) (any(x[middle.regions] > 0.5)) & (all(x[-middle.regions] < 0.5)))]
      
    }else{
      
      project.to.middle <- NULL
    }
    
    
    # Project to back and middle
    if(!is.null(middle.regions) & !is.null(back.regions)){
      
      project.to.back.middle <- (1:J)[apply(q_tilde_001,
                                            1,
                                            function(x) (any(x[middle.regions] > 0.5)) & (any(x[back.regions] > 0.5))
                                            & (all(x[front.regions] < 0.5))
      )
      ]
      
    }else{
      
      project.to.back.middle <- NULL
    }
    
    # Project to front and middle
    if(!is.null(front.regions) & !is.null(middle.regions)){
      
      project.to.front.middle <- (1:J)[apply(q_tilde_001,
                                             1,
                                             function(x) (any(x[front.regions] > 0.5)) & (any(x[middle.regions] > 0.5))
                                             & (all(x[back.regions] < 0.5))
      )
      ]
      
    }else{
      
      project.to.front.middle <- NULL
    }
    
    # Project to front and back
    if(!is.null(front.regions) & !is.null(back.regions)){
      
      project.to.front.back <- (1:J)[apply(q_tilde_001,
                                           1,
                                           function(x) (any(x[front.regions] > 0.5)) & (any(x[back.regions] > 0.5))
                                           & (all(x[middle.regions] < 0.5))
      )
      ]
      
    }else{
      
      project.to.front.back <- NULL
    }
    
    
    # Project to all three
    if(!is.null(front.regions) & !is.null(middle.regions) & !is.null(back.regions)){
      
      project.to.all.area <- (1:J)[apply(q_tilde_001,
                                         1,
                                         function(x) (any(x[front.regions] > 0.5)) & (any(x[middle.regions] > 0.5))
                                         & (any(x[back.regions] > 0.5))
      )
      ]
      
    }else{
      
      project.to.all.area <- NULL
    }
    
    
    # Cluster labels
    cluster.label.summary <- data.frame(cluster = 1:J,
                                        
                                        unicast = ifelse((1:J) %in% unicast.index, 'Yes', 'No'),
                                        bicast.and.more = ifelse((1:J) %in% bicast.and.more.index, 'Yes', 'No'),
                                        broadcast = ifelse((1:J) %in% broadcast.index, 'Yes', 'No'),
                                        
                                        project.to.back = ifelse((1:J) %in% project.to.back, 'Yes', 'No'),
                                        project.to.front = ifelse((1:J) %in% project.to.front, 'Yes', 'No'),
                                        project.to.middle = ifelse((1:J) %in% project.to.middle, 'Yes', 'No'),
                                        project.to.back.middle = ifelse((1:J) %in% project.to.back.middle, 'Yes', 'No'),
                                        project.to.front.middle = ifelse((1:J) %in% project.to.front.middle, 'Yes', 'No'),
                                        project.to.all.area = ifelse((1:J) %in% project.to.all.area, 'Yes', 'No'),
                                        project.to.front.back = ifelse((1:J) %in% project.to.front.back, 'Yes', 'No')
                                        
    )
    
    
    # Reorder original label
    original.label <- c(project.to.front,
                        project.to.front.middle,
                        project.to.middle,
                        project.to.back.middle,
                        project.to.back,
                        project.to.front.back,
                        project.to.all.area)
    
    
    
    
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
    
    
  }else{
    
    # In the case with no reordering of clusters
    cluster.label.summary <- data.frame(cluster = 1:J,
                                        
                                        unicast = ifelse((1:J) %in% unicast.index, 'Yes', 'No'),
                                        bicast.and.more = ifelse((1:J) %in% bicast.and.more.index, 'Yes', 'No'),
                                        broadcast = ifelse((1:J) %in% broadcast.index, 'Yes', 'No')
                                        
    )
    
  }
  
  # Summary of projection probabilities
  loop.result <- lapply(1:J, function(j){
    
    data.frame(cluster = paste('cluster', j),
               region = regions.name,
               projection.mean = proj_prob_mean[j,],
               projection.lower = proj_prob_lower[j,],
               projection.upper = proj_prob_upper[j,])
  })
  
  estimated.projection.df <- do.call(rbind, loop.result)
  
  
  #---------------------------------------- Plot of projection strengths ------------------------------------
  
  
  
  cluster_label_vec <- sapply(1:J,
                              function(j){
                                
                                colnames(cluster.label.summary)[2:4][which(cluster.label.summary[j,2:4]=='Yes')[1]]
                              })
  
  # Label the number of projecting regions
  estimated.projection.df$class <- rep(cluster_label_vec, each = R)
  estimated.projection.df$class <- factor(estimated.projection.df$class,
                                          levels = colnames(cluster.label.summary)[2:4])
  
  # Manually set the color
  col <- c('red','green','blue')
  unique.class <- unique(estimated.projection.df$class)
  col.i <- sapply(1:length(unique.class),
                  function(i) col[which(colnames(cluster.label.summary)[2:4] == unique.class[i])])
  
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
          plot.title = element_text(size=plot.title.size))+
    scale_color_manual(values = col.i)
  
  
  #----------------------------------------- If the location of the brain regions are given -----------------------
  
  if(isTRUE(relabel)){
    
    # Create an output space
    plot.by.region <- list()
    
    #-------------------------------------- Clusters projecting to the back ----------------------------------------
    
    if(!is.null(project.to.back) & (length(project.to.back) != 0)){
      
      
      unique.class <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.back]))
      
      unique.class <- sort(unique(unique.class$class))
      
      
      col.i <- sapply(1:length(unique.class),
                      function(i) col[which(colnames(cluster.label.summary)[2:4] == unique.class[i])])
      
      plot.by.region$back <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.back])) %>%
        ggplot(mapping = aes(x = factor(region, levels = regions.name),
                             y = projection.mean,
                             group = cluster,
                             color = class)) +
        geom_line()+
        geom_point()+
        geom_errorbar(mapping = aes(ymin = projection.lower,
                                    ymax = projection.upper),
                      width = 0.1)+
        theme_bw()+
        ylim(c(0,1))+
        ylab('projection probabilities')+
        ggtitle('Clusters projecting to the back')+
        xlab('region')+
        facet_wrap(~cluster)+
        
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              plot.title = element_text(size=10),
              legend.text = element_text(size=10))+
        scale_color_manual(values = col.i)
      
    }else{
      
      plot.by.region$back <- NA
    }
    
    
    
    
    #-------------------------------------- Clusters projecting to the front ----------------------------------------
    
    if(!is.null(project.to.front) & (length(project.to.front) != 0)){
      
      
      
      unique.class <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.front])) %>%
        select(class) %>%
        pull() %>%
        unique()%>%
        sort()
      
      
      col.i <- sapply(1:length(unique.class),
                      function(i) col[which(colnames(cluster.label.summary)[2:4] == unique.class[i])])
      
      
      plot.by.region$front <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.front])) %>%
        ggplot(mapping = aes(x = factor(region, levels = regions.name),
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
        ggtitle('Clusters projecting to the front')+
        xlab('region')+
        facet_wrap(~cluster)+
        
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              plot.title = element_text(size=10),
              legend.text = element_text(size=10))+
        scale_color_manual(values = col.i)
      
    }else{
      
      plot.by.region$front <- NA
    }
    
    
    
    
    #-------------------------------------- Clusters projecting to the middle ----------------------------------------
    
    if(!is.null(project.to.middle) & (length(project.to.middle) != 0)){
      
      
      
      unique.class <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.middle])) %>%
        select(class) %>%
        pull() %>%
        unique()%>%
        sort()
      
      
      col.i <- sapply(1:length(unique.class),
                      function(i) col[which(colnames(cluster.label.summary)[2:4] == unique.class[i])])
      
      
      plot.by.region$middle <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.middle])) %>%
        ggplot(mapping = aes(x = factor(region, levels = regions.name),
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
        ggtitle('Clusters projecting to the middle')+
        xlab('region')+
        facet_wrap(~cluster)+
        
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              plot.title = element_text(size=10),
              legend.text = element_text(size=10))+
        scale_color_manual(values = col.i)
      
      
    }else{
      
      plot.by.region$middle <- NA
    }
    
    
    
    
    #-------------------------------------- Clusters projecting to the front and middle ----------------------------------------
    
    if(!is.null(project.to.front.middle) & (length(project.to.front.middle) != 0)){
      
      
      unique.class <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.front.middle])) %>%
        select(class) %>%
        pull() %>%
        unique()%>%
        sort()
      
      
      col.i <- sapply(1:length(unique.class),
                      function(i) col[which(colnames(cluster.label.summary)[2:4] == unique.class[i])])
      
      
      plot.by.region$front_middle <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.front.middle])) %>%
        ggplot(mapping = aes(x = factor(region, levels = regions.name),
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
        ggtitle('Clusters projecting to the middle and front')+
        xlab('region')+
        facet_wrap(~cluster)+
        
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              plot.title = element_text(size=10),
              legend.text = element_text(size=10))+
        scale_color_manual(values = col.i)
      
    }else{
      
      plot.by.region$front_middle <- NA
    }
    
    
    
    
    #-------------------------------------- Clusters projecting to the middle and back ----------------------------------------
    
    if(!is.null(project.to.back.middle) & (length(project.to.back.middle) != 0)){
      
      
      unique.class <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.back.middle])) %>%
        select(class) %>%
        pull() %>%
        unique()%>%
        sort()
      
      
      col.i <- sapply(1:length(unique.class),
                      function(i) col[which(colnames(cluster.label.summary)[2:4] == unique.class[i])])
      
      
      plot.by.region$middle_back <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.back.middle])) %>%
        ggplot(mapping = aes(x = factor(region, levels = regions.name),
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
        ggtitle('Clusters projecting to the middle and back')+
        xlab('region')+
        facet_wrap(~cluster)+
        
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              plot.title = element_text(size=10),
              legend.text = element_text(size=10))+
        scale_color_manual(values = col.i)
      
    }else{
      
      plot.by.region$middle_back <- NA
    }
    
    
    
    
    
    #-------------------------------------- Clusters projecting to the front and back ----------------------------------------
    if(!is.null(project.to.front.back) & (length(project.to.front.back) != 0)){
      
      
      unique.class <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.front.back])) %>%
        select(class) %>%
        pull() %>%
        unique()%>%
        sort()
      
      
      col.i <- sapply(1:length(unique.class),
                      function(i) col[which(colnames(cluster.label.summary)[2:4] == unique.class[i])])
      
      
      plot.by.region$front_back <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.front.back])) %>%
        ggplot(mapping = aes(x = factor(region, levels = regions.name),
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
        ggtitle('Clusters projecting to front and back')+
        xlab('region')+
        facet_wrap(~cluster)+
        
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              plot.title = element_text(size=10),
              legend.text = element_text(size=10))+
        scale_color_manual(values = col.i)
      
    }else{
      
      plot.by.region$front_back <- NA
    }
    
    
    
    
    
    #---------------------------------------------- Project to all three regions ----------------------------------
    if(!is.null(project.to.all.area) & (length(project.to.all.area) != 0)){
      
      unique.class <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.all.area])) %>%
        select(class) %>%
        pull() %>%
        unique()%>%
        sort()
      
      
      col.i <- sapply(1:length(unique.class),
                      function(i) col[which(colnames(cluster.label.summary)[2:4] == unique.class[i])])
      
      
      plot.by.region$all_regions <- estimated.projection.df %>%
        filter(cluster %in% paste('cluster', (1:J)[original.label %in% project.to.all.area])) %>%
        ggplot(mapping = aes(x = factor(region, levels = regions.name),
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
        ggtitle('Clusters projecting to front, middle and back')+
        xlab('region')+
        facet_wrap(~cluster)+
        
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              plot.title = element_text(size=10),
              legend.text = element_text(size=10))+
        scale_color_manual(values = col.i)
      
    }else{
      
      plot.by.region$all_regions <- NA
    }
    
    
    
  }else{
    
    # In the case when the clusters are not reordered
    plot.by.region <- NA
  }
  
  
  #--------------------------------------- Plot of q_tilde ----------------------------------
  
  if(relabel == TRUE){
    
    q_tilde_001 <- q_tilde_001[original.label,]
  }
  
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
              'plot.by.region' = plot.by.region,
              'q_tilde_001' = q_tilde_001,
              'q_tilde_plot' = q_tilde_plot,
              'relabel' = relabel)
  )
}

