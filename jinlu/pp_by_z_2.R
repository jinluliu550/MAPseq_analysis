
pp.standard.ordering2 <- function(Y,
                                  Z,
                                  regions.name = NULL,
                                  mouse.index,
                                  title = ''){
  
  
  J <- max(unlist(Z))
  R <- nrow(Y[[1]])
  M <- length(Y)
  
  if(is.null(regions.name)){
    
    regions.name <- paste('region', 1:R)
  }
  
  Y.cbind <- do.call(cbind, Y)
  Y.prop.cbind <- matrix(as.vector(Y.cbind)/rep(colSums(Y.cbind), each = R),
                         nrow = R)
  
  Y.prop <- lapply(1:M, function(m) matrix(Y[[m]]/rep(colSums(Y[[m]]), each = R),
                                           nrow = R))
  
  
  # Re-ordered projection strength
  df0 <- lapply(1:J,
                function(j){
                  
                  
                  mx0 <- matrix(Y.prop.cbind[,which(unlist(Z)==j)],
                                nrow = R)
                  
                  strong.proj.j <- which(rowSums(mx0) == max(rowSums(mx0)))[1]
                  
                  df0_list <- lapply(1:M, function(m){
                    
                    if(length(which(Z[[m]]==j)) != 0){
                      
                      mx0_m <- matrix(Y.prop[[m]][,which(Z[[m]]==j)],
                                      nrow = R)
                      
                      mx0_m[,order(mx0_m[strong.proj.j,])]
                    }
                    
                  })
                  
                  do.call(cbind, df0_list)
                  
                })
  
  # Mouse index
  mouse.index.df <- lapply(1:J,
                           function(j){
                  
                  df0_list <- lapply(1:M, function(m){
                    
                    
                    if(length(which(Z[[m]]==j)) != 0){
                      
                      
                      rep(m, length(which(Z[[m]]==j)))
                    }
                    
                  })
                  
                  unlist(df0_list)
                  
                })
  
  mouse.index.df <- unlist(mouse.index.df)
  
  
  
  
  df0 <- do.call(cbind, df0)
  
  N.j <- sapply(1:J,
                function(j) length(which(unlist(Z)==j)))
  
  # Convert to data frame for plotting
  Y.prop.df <- data.frame(region = rep(regions.name,
                                       ncol(Y.prop.cbind)),
                          
                          neuron = rep(1:ncol(Y.prop.cbind),
                                       each = length(regions.name)),
                          
                          projection.strength = as.vector(df0),
                          
                          mouse = factor(rep(mouse.index.df, each = R),
                                         levels = 1:M))
  
  Y.prop.df %>%
    ggplot(mapping = aes(x = factor(region, levels = regions.name),
                         y = neuron,
                         fill = mouse))+
    
    geom_tile(mapping = aes(alpha = projection.strength))+
    theme_bw()+
    xlab('region')+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size=12),
          plot.background = element_blank() ,
          panel.grid.major = element_blank() ,
          panel.grid.minor = element_blank() ,
          panel.border = element_blank() ,
          panel.background = element_blank()) +
    #draws x and y axis line  
    theme(axis.line = element_line(color = 'black'))+
    ylab('neurons')+
    
    geom_hline(yintercept = cumsum(N.j)[-length(cumsum(N.j))],
               color = 'blue',
               linetype = 'dashed',
               linewidth = 0.05)+
    
    ggtitle(title)
  
  
  
}