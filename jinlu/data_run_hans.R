Hans_data1 <- read.csv('./data/Han-data/han_brain4_gh.csv')
Hans_data2 <- read.csv('./data/Han-data/han_brain5_gh.csv')
Hans_data3 <- read.csv('./data/Han-data/han_brain6_gh.csv')

data_Hans <- list(t(Hans_data1),
                  t(Hans_data2),
                  t(Hans_data3))


mcmc_all_hans <- mcmc_run_all(Y = data_Hans,
                             J = 20,
                             number_iter = 8000,
                             thinning = 5,
                             burn_in = 3000,
                             adaptive_prop = 0.1,
                             print_Z = TRUE,
                             
                             
                             a_gamma = 300,
                             b_gamma = 10,
                             a_alpha = 1/5,
                             b_alpha = 1/2,
                             num.cores = 10)


