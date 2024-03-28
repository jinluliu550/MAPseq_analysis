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


dlso_cluster_estimate <- function(mcmc_run_all_output){
  
  # trace of z
  z_trace <- mcmc_run_all_output$Z_output
  
  # trace of z in a matrix
  z_trace_mx <- lapply(1:length(z_trace),
                       function(t) matrix(unlist(z_trace[[t]]),
                                          nrow = 1))
  
  
  z_trace_mx <- do.call(rbind, z_trace_mx)
  
  z_point_estimate <- salso::dlso(truth = z_trace_mx)
  
  # Cut into data
  M <- mcmc_run_all_output$M
  C <- mcmc_run_all_output$C
  
  C_cumsum <- c(0, cumsum(C))
  
  z_point_estimate <- lapply(1:M,
                             function(m) z_point_estimate[(C_cumsum[m]+1):C_cumsum[m+1]])
  
  # Reorder
  z_point_estimate_r <- lapply(1:M,
                               function(m){
                                 
                                 sapply(1:C[m],
                                        function(c){
                                          
                                          which(sort(unique(unlist(z_point_estimate))) == z_point_estimate[[m]][c])
                                        })
                               })
  
  return(z_point_estimate_r)
}



# Function which compares minVI and dlso

opt.clustering.comb <- function(mcmc_run_all_output,
                                post_similarity,
                                run.on.pc = T){
  
  
  
  #------------------------------------------------ Point estimate of z and trace of z from MCMC ------------------------------------------------
  
  # Result from minvi
  opt.clust.minvi <- opt.clustering(mcmc_run_all_output = mcmc_run_all_output,
                                         post_similarity = post_similarity)
  
  # Result from dlso
  opt.clust.dlso <- dlso_cluster_estimate(mcmc_run_all_output = mcmc_run_all_output)
  
  # Trace of z
  z.trace <- mcmc_run_all_output$Z_output
  
  #---------------------------------- Difference between each sample and the minvi estimated allocation -----------------------------------------
  
  if(run.on.pc == FALSE){
    
    cl <- makeCluster(num.cores,
                      type = "FORK",
                      .packages = 'clevr')
  }else{
    
    cl <- makeCluster(num.cores,
                      .packages = 'clevr')
  }
  
  vi.minvi <- pblapply(1:length(z.trace),
                          cl = cl,
                          FUN = function(t){
                            
                            # Similarity
                            vi.minvi.t <- variation_info(unlist(opt.clust.minvi),
                                                         unlist(z.trace[[t]]))
                            
                            vi.minvi.t
                          }
  )
  
  stopCluster(cl)
  
  #--------------------------------- Difference between each sample and dlso estimated allocation ------------------------------------------------
  
  if(run.on.pc == FALSE){
    
    cl <- makeCluster(num.cores,
                      type = "FORK",
                      .packages = 'clevr')
  }else{
    
    cl <- makeCluster(num.cores,
                      .packages = 'clevr')
  }
  
  vi.dlso <-  pblapply(1:length(z.trace),
                       cl = cl,
                       FUN = function(t){
                         
                         # Similarity
                         vi.dlso.t <- variation_info(unlist(opt.clust.dlso),
                                                      unlist(z.trace[[t]]))
                         
                         vi.dlso.t
                       }
  )
  
  stopCluster(cl)
  
  #-------------------------------- Return item ----------------------------------
  
  # In the case when minvi gives smaller VI than dlso
  if(mean(unlist(vi.minvi)) < mean(unlist(vi.dlso))){
    
    # we export the minvi point estimate of z
    return(opt.clust.minvi)
    
  }else{
    
    # otherwise export the dlso point estimate of z
    return(opt.clust.dlso)
  }
  
}