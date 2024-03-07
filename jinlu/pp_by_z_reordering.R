
#' Heat-map to show projection stength of each neuron with reordered clusters
#'

pp.reordered <- function(qj_estimate,
                         Y,
                         Z,
                         regions.name = NULL){

  R <- nrow(Y[[1]])

  if(is.null(regions.name)){

    regions.name <- paste('region', 1:R)
  }

  strong.proj <- apply(qj_estimate,
                       1,
                       function(x) which(x == max(x)))

  df0 <- lapply(1:ncol(qj_estimate),
                function(i){

                  # Clusters with the strongest projection to region i
                  strong.proj.i <- which(strong.proj==i)

                  if(length(strong.proj.i) != 0){

                    qj_estimate_i <- matrix(qj_estimate[strong.proj.i,],
                                            nrow = length(strong.proj.i))

                    output <- strong.proj.i[order(-qj_estimate_i[,i])]
                  }else{

                    output <- NULL
                  }


                  }
                )

  # Number of neurons in each cluster
  N.j <- sapply(rev(unlist(df0)), function(i) length(which(unlist(Z)==i)))

  # Reordered Y
  Y.cbind <- do.call(cbind,
                     Y)

  Y.reordered <- lapply(rev(unlist(df0)),

                        function(i){

                          Y.cbind[,which(unlist(Z)==i)]
                        })

  Y.reordered <- do.call(cbind,
                         Y.reordered)

  Y.prop.reordered <- matrix(as.vector(Y.reordered)/rep(colSums(Y.reordered),
                                                        each = length(regions.name)),

                             nrow = length(regions.name))

  Y.prop.by.z <- data.frame(region = rep(regions.name,
                                         ncol(Y.reordered)),

                            cell = rep(1:ncol(Y.reordered),
                                       each = length(regions.name)),

                            cellwise.pp = as.vector(Y.prop.reordered))


  Y.prop.by.z %>%
    ggplot(mapping = aes(x = factor(region, levels = regions.name),
                         y = cell,
                         fill = cellwise.pp))+
    geom_tile()+
    theme_bw()+
    scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                         values = scales::rescale(c(0, 0.25, 1)),
                         limits = c(0,1))+
    guides(fill=guide_legend(title="Projection probability"))+
    xlab('region')+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size=12))+
    ylab('neurons')+

    geom_hline(yintercept = cumsum(N.j)[-length(cumsum(N.j))],
               color = 'blue',
               linetype = 'dashed',
               size = 0.05)+
    guides(fill=guide_legend(title="Projection strength"))


}

