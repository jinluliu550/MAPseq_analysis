
# Compare different methods

#------------------------------- k-means ------------------------------------

k_20_difference <- variation_info(unlist(clust20_r),
                                  unlist(clust20_r_5))/
  
  log(base = 2, x = sum(C))




k_30_difference <- variation_info(unlist(clust30_r),
                                  unlist(clust30_r_5))/
  
  log(base = 2, x = sum(C))



k_40_difference <- variation_info(unlist(clust40_r),
                                  unlist(clust40_r_5))/
  
  log(base = 2, x = sum(C))



#------------------------------ Bayesian ---------------------------------


bayesian_difference <- variation_info(unlist(hans_Z_reordered_5),
                                      unlist(hans_Z_reordered))/
  
  log(base = 2, x = sum(C))



#------------------------------ Binomial ---------------------------------

binomial_difference <- variation_info(unlist(hans_binomial_reorder5$allocation),
                                      unlist(hans_binomial_reorder1$allocation))/
  
  log(base = 2, x = sum(C))


binomial_difference_pe <- variation_info(unlist(hans_binomial_reorder5_pe$allocation),
                                         unlist(hans_binomial_reorder1_pe$allocation))/
  
  log(base = 2, x = sum(C))
