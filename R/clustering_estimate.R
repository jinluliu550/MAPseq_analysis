#' Function to obtain clustering estimate which minimizes the variation of information
#'

#-- Update the trace of z such that all the cluster labels are consecutive
z_trace_updated <- function(mcmc_run_all_output){
  
  # Dimensions
  C <- mcmc_run_all_output$C
  M <- mcmc_run_all_output$M
  
  # Updated allocations
  z_trace_updated <- lapply(1:length(mcmc_run_all_output$Z_output),
                            function(t){
                              
                              lapply(1:M,
                                     function(m){
                                       
                                       sapply(1:C[m],
                                              function(c) which(sort(unique(unlist(mcmc_run_all_output$Z_output[[t]][[m]][c]))) == 
                                                                  mcmc_run_all_output$Z_output[[t]][[m]][c]))
                                     })
                            })
  
  return(z_trace_updated)
}


opt.clustering <- function(z_trace,
                           post_similarity,
                           max.k=NULL){

  M <- length(z_trace[[1]])
  C <- sapply(1:M, function(m) length(z_trace[[1]][[m]]))
  

  ##-- Unlist for each iteration and put into a matrix
  Z_output_all <- lapply(1:length(z_trace),

                         function(t) matrix(unlist(z_trace[[t]]),
                                            nrow = 1)
                         )


  mcmc.z.matrix <- do.call(rbind, Z_output_all)

  opt.clust <- minVI(psm = post_similarity$psm.combined,
                     cls.draw = mcmc.z.matrix,
                     method = 'all',
                     max.k=max.k)

  opt.clust <- opt.clust$cl[1,]

  # Cumulative sum of C
  C_cumsum <- c(0,cumsum(C))

  opt.clust <- lapply(1:M, function(m) opt.clust[(C_cumsum[m]+1):C_cumsum[m+1]])

  return(opt.clust)
}


salso_cluster_estimate <- function(z_trace,
                                   max.k=NULL){
  
  
  M <- length(z_trace[[1]])
  C <- sapply(1:M, function(m) length(z_trace[[1]][[m]]))
  C_cumsum <- c(0, cumsum(C))
  
  
  ##-- Unlist for each iteration and put into a matrix
  Z_output_all <- lapply(1:length(z_trace),
                         
                         function(i) matrix(unlist(z_trace[[i]]),
                                            nrow = 1)
  )
  
  
  mcmc.z.matrix <- do.call(rbind, Z_output_all)
  
  ##-- Point estimate of z
  if(is.null(max.k)){max.k=0}
  z_point_estimate <- salso::salso(x = mcmc.z.matrix, 
                                   maxNClusters = max.k)
  
  ##-- Convert to a list
  z_point_estimate <- lapply(1:M,
                             function(m) z_point_estimate[(C_cumsum[m]+1):C_cumsum[m+1]])
  
  
  return(z_point_estimate)
}



# Function which compares minVI and dlso

opt.clustering.comb <- function(z_trace,
                                post_similarity,
                                max.k=NULL){
  
  
  
  #------------------------------------------------ Point estimate of z and trace of z from MCMC ------------------------------------------------
  
  # Result from minvi
  opt.clust.minvi <- opt.clustering(z_trace = z_trace,
                                    post_similarity = post_similarity,
                                    max.k=max.k)
  
  # Result from dlso
  opt.clust.dlso <- salso_cluster_estimate(z_trace = z_trace,
                                           max.k=max.k)
  
  
  #---------------------------------- Difference between each sample and the minvi estimated allocation -----------------------------------------
  
  vi.min <- sapply(1:length(z_trace),
                   function(t){
                     
                     vi.minvi.t <- clevr::variation_info(true = unlist(opt.clust.minvi),
                                                         pred = unlist(z_trace[[t]]), base=2)
                     
                     vi.minvi.t
                   })
  
  
  #--------------------------------- Difference between each sample and dlso estimated allocation ------------------------------------------------
  
  vi.dlso <- sapply(1:length(z_trace),
                    function(t){
                      
                      # Similarity
                      vi.dlso.t <- variation_info(unlist(opt.clust.dlso),
                                                  unlist(z_trace[[t]]), base=2)
                      
                      vi.dlso.t
                    })
  
  
  #-------------------------------- Return item ----------------------------------
  
  # In the case when minvi gives smaller VI than dlso
  if(mean(vi.min) < mean(vi.dlso)){
    
    # we export the minvi point estimate of z
    return(opt.clust.minvi)
    
  }else{
    
    # otherwise export the dlso point estimate of z
    return(opt.clust.dlso)
  }
  
}