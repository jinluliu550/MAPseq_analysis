#' Posterior similarity matrix
#'

similarity_matrix <- function(mcmc_run_all_output){

  # Input from MCMC
  M <- mcmc_run_all_output$M
  C <- mcmc_run_all_output$C

  # Allocations
  allocation_list <- mcmc_run_all_output$Z_output



  ##-- Cumulative sums of cells
  C_cum <- c(0, cumsum(C))

  ##-------------------- Similarity between cells from the same dataset and across datasets --------------

  psm.within <- NULL
  for(m in 1:M) psm.within[[m]] <- NULL

  ##-- For all possible combinations of datasets
  for(m1 in rev(1:M)){
    for(m2 in 1:m1){

      ##-- For each iteration
      loop.result <- lapply(1:length(allocation_list), function(iter){

        ##-- Allocations of cells in dataset m1 and iteration = iter
        allocation_list_m1_iter <- allocation_list[[iter]][[m1]]
        allocation_list_m2_iter <- allocation_list[[iter]][[m2]]

        psm.empty <- matrix(0, nrow = C[m1], ncol = C[m2])
        for(i in 1:C[m1]){
          psm.empty[i,] <- psm.empty[i,] + ifelse(allocation_list_m2_iter == allocation_list_m1_iter[i],
                                                  1,
                                                  0)
        }

        psm.empty
      })

      ##-- Sum up the matrices
      loop.result <- Reduce('+', loop.result)
      psm.within[[m1]][[m2]] <- loop.result/length(allocation_list)

    }
  }

  ##-------------------------------- Computing a combined matrix -------------------------

  combined_matrix <- matrix(0, nrow = max(C_cum), ncol = max(C_cum))

  for(m1 in 1:M){
    for(m2 in 1:m1){

      combined_matrix[(C_cum[m1]+1):C_cum[m1+1], (C_cum[m2]+1):C_cum[m2+1]] <- psm.within[[m1]][[m2]]
    }

    if(m1 != M){
      for(m2 in (m1+1):M){

        combined_matrix[(C_cum[m1]+1):C_cum[m1+1], (C_cum[m2]+1):C_cum[m2+1]] <- t(psm.within[[m2]][[m1]])
      }
    }
  }

  ##-------------------------------- Return both -----------------------------------------

  return(list('psm.within' = psm.within,
              'psm.combined' = combined_matrix))

}


