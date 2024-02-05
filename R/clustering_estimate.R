#' Function to obtain clustering estimate which minimizes the variation of information
#'

opt.clustering <- function(mcmc_run_all_output,
                           post_similarity){

  C <- mcmc_run_all_output$C
  M <- mcmc_run_all_output$M

  ##-- Unlist for each iteration and put into a matrix
  Z_output_all <- lapply(1:length(mcmc_run_all_output$Z_output),

                         function(i) matrix(unlist(mcmc_run_all_output$Z_output[[i]]),
                                            nrow = 1)
                         )


  mcmc.z.matrix <- do.call(rbind, Z_output_all)

  opt.clust <- minVI(psm = post_similarity$psm.combined,
                     cls.draw = mcmc.z.matrix,
                     method = 'all')

  opt.clust <- opt.clust$cl[1,]

  # Cumulative sum of C
  C_cumsum <- c(0,cumsum(C))

  opt.clust <- lapply(1:M, function(m) opt.clust[(C_cumsum[m]+1):C_cumsum[m+1]])

  return(opt.clust)
}
