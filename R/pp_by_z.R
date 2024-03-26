#' Heat-map to show projection strength
#'

pp.standard.ordering <- function(Y,
                                 Z,
                                 regions.name = NULL,
                                 title = ''){


  J <- max(unlist(Z))
  R <- nrow(Y[[1]])

  if(is.null(regions.name)){

    regions.name <- paste('region', 1:R)
  }

  Y.cbind <- do.call(cbind, Y)
  Y.prop.cbind <- matrix(as.vector(Y.cbind)/rep(colSums(Y.cbind), each = R),
                         nrow = R)

  # Re-ordered projection strength
  df0 <- lapply(1:J,
                function(j){

                  matrix(Y.prop.cbind[,which(unlist(Z)==j)],
                         nrow = R)
                })

  df0 <- do.call(cbind, df0)

  # Number of neurons in each cluster
  N.j <- sapply(1:J, function(j) length(which(unlist(Z)==j)))

  # Convert to data frame for plotting
  Y.prop.df <- data.frame(region = rep(regions.name,
                                       ncol(Y.prop.cbind)),

                          neuron = rep(1:ncol(Y.prop.cbind),
                                       each = length(regions.name)),

                          projection.strength = as.vector(df0))


  Y.prop.df %>%
    ggplot(mapping = aes(x = factor(region, levels = regions.name),
                         y = neuron,
                         fill = projection.strength))+
    geom_tile()+
    theme_bw()+
    scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                         values = scales::rescale(c(0, 0.25, 1)),
                         limits = c(0,1))+
    guides(fill=guide_legend(title="projection strength"))+
    xlab('region')+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size=12))+
    ylab('neurons')+

    geom_hline(yintercept = cumsum(N.j)[-length(cumsum(N.j))],
               color = 'blue',
               linetype = 'dashed',
               linewidth = 0.05)+

    ggtitle(title)


}

