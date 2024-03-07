# Running on the barSeq data


# Load data
data1 <- read.csv('./data/BRAIN_C9.csv')
data2 <- read.csv('./data/BRAIN_C14.csv')
data3 <- read.csv('./data/BRAIN_C28.csv')

# Load data
load('data/psm_barseq.RData')
load('data/mcmc_all_barseq.RData')

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



psm_barseq <- similarity_matrix(mcmc_run_all_output = mcmc_all_barseq,
                                num.cores = 10,
                                run.on.pc = FALSE)


# Plot of posterior similarity matrix
plotpsm(psm.ind = psm_barseq$psm.within,
        psm.tot = psm_barseq$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)


plotpsm(psm.ind = psm_barseq$psm.within,
        psm.tot = psm_barseq$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')


opt.clust0.barseq <- opt.clustering(mcmc_run_all_output = mcmc_all_barseq,
                                    post_similarity = psm_barseq)

mcmc_unique_barseq <- mcmc_run_post(mcmc_run_all_output = mcmc_all_barseq,
                                    Z = opt.clust0.barseq,
                                    thinning = 5,
                                    burn_in = 1500,
                                    number_iter = 3000,
                                    Y = data_barseq,
                                    a_gamma = 500,
                                    b_gamma = 10,
                                    regions.name = rownames(data_barseq[[1]]))
                             
