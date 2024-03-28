EC8_new <- read.csv('./data/EC8_new/all_brains_setE_JL.csv')

# Convert to list
digit.list <- list(one = 1,
                   two = 2,
                   three = 3,
                   four = 4,
                   five = 5,
                   six = 6)

EC8_new_list <- lapply(1:6, 
                       function(m) EC8_new %>%
                         filter(brain == names(digit.list)[m]))

EC8_new <- lapply(1:6,
                  function(m) t(round(EC8_new_list[[m]][,1:8], 0)))

EC8_EC_label <- lapply(1:6,
                       function(m) EC8_new_list[[m]]$EC)

mcmc_all_EC8 <- mcmc_run_all(Y = EC8_new,
                             J = 150,
                             number_iter = 8000,
                             thinning = 5,
                             burn_in = 3000,
                             adaptive_prop = 0.1,
                             print_Z = TRUE,
                             
                             
                             a_gamma = 500,
                             b_gamma = 10,
                             a_alpha = 1/5,
                             b_alpha = 1/2,
                             num.cores = 10)

