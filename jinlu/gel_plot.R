gel_plot <- function(Y){
  
  M <- length(Y)
  
  plot.output <- NULL
  
  for(m in 1:M){
    
    Y_m <- Y[[m]]
    
    Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
    
    # Index of the strongest projecting region
    strongest_projecting_region <- apply(Y_m, 2, function(x) which(x == max(x))[1])
    
    df_list <- lapply(1:nrow(Y_m),
                      function(r){
                        
                        if(length(which(strongest_projecting_region == r)) != 0){
                          
                          
                         matrix0 <- matrix(Y_m_prop[,which(strongest_projecting_region == r)],
                                           nrow = nrow(Y_m_prop))
                         
                         matrix0[,order(-Y_m_prop[r,which(strongest_projecting_region == r)])]

                          
                        }else{
                          
                        }
                        
                      })
    
    df <- do.call(cbind, df_list)
    
    df1 <- data.frame(projection_strength = as.vector(df),
                      region = rownames(Y_m),
                      cell_index = rep(1:ncol(df), each = nrow(Y_m)))
    if(m != M){
      
      plot.output[[m]] <- df1 %>%
        ggplot(mapping = aes(x = factor(region, levels = rownames(Y_m)),
                             y = cell_index,
                             fill = projection_strength))+
        geom_tile()+
        theme_bw()+
        scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                             values = scales::rescale(c(0, 0.25, 1)),
                             limits = c(0,1))+
        guides(fill=guide_legend(title="projection strength"))+
        xlab('region')+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              plot.title = element_text(size=12),
              legend.position = 'none',
              panel.grid = element_blank(),
              panel.border = element_blank())+
        ylab('neurons')+
        ggtitle(paste('mouse', m))+
        ylim(c(0, ncol(Y_m)))
      
    }else{
      
      
      plot.output[[m]] <- df1 %>%
        ggplot(mapping = aes(x = factor(region, levels = rownames(Y_m)),
                             y = cell_index,
                             fill = projection_strength))+
        geom_tile()+
        theme_bw()+
        scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                             values = scales::rescale(c(0, 0.25, 1)),
                             limits = c(0,1))+
        guides(fill=guide_legend(title="projection strength"))+
        xlab('region')+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              plot.title = element_text(size=12),
              panel.grid = element_blank(),
              panel.border = element_blank())+
        ylab('neurons')+
        ggtitle(paste('mouse', m))+
        ylim(c(0, ncol(Y_m)))
    }
    
    
  }
  
  plot.output
}

# Histograms of brain regions
histogram_empirical <- function(Y,
                                region1,
                                region2){
  
  # region 1
  list_region1 <- lapply(1:length(Y),
                  function(m){
                    
                    
                    Y_m <- Y[[m]]
                    
                    Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
                    
                    data.frame(mouse = m,
                               region = region1,
                               projection_strength = as.vector(Y_m_prop[region1,]),
                               neuron = factor(1:ncol(Y_m)))
                    
                  })
  
  df_region1 <- do.call(rbind, list_region1)
  
  # region 2
  list_region2 <- lapply(1:length(Y),
                         function(m){
                           
                           
                           Y_m <- Y[[m]]
                           
                           Y_m_prop <- apply(Y_m, 2, function(x) x/sum(x))
                           
                           data.frame(mouse = m,
                                      region = region2,
                                      projection_strength = as.vector(Y_m_prop[region2,]),
                                      neuron = factor(1:ncol(Y_m)))
                           
                         })
  
  df_region2 <- do.call(rbind, list_region2)
  
  
  df <- rbind(df_region1,
              df_region2)
  
  df$mouse <- factor(paste('mouse', df$mouse),
                     levels = paste('mouse', 1:length(Y)))
  
  # Histogram
  histogram <- ggplot(df, aes(x=projection_strength))+
    geom_histogram(binwidth = 0.05)+
    facet_grid(vars(region), vars(mouse))+
    xlab('y/N')+
    theme_bw()
  
  # Scatter plot
  df2 <- df %>%
    pivot_wider(names_from = region,
                values_from = projection_strength) 
  
  plot(x = as.vector(df2[,3]),
       y = as.vector(df2[,4]))
  
}
