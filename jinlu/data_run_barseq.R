# Running on the barSeq data


# Load data
data1 <- read.csv('./data/BRAIN_C9.csv')
data2 <- read.csv('./data/BRAIN_C14.csv')
data3 <- read.csv('./data/BRAIN_C28.csv')

# Round to integers
data1 <- round(data1, 0)
data2 <- round(data2, 0)
data3 <- round(data3, 0)

data_barseq <- list(t(data1),
                    t(data2),
                    t(data3))

# Save data
# save(data_barseq, file = './data/bar_seq.RData')

# MCMC run
mcmc_all_barseq <- mcmc_run_all(Y = data_barseq,
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


