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

# MCMC of unique parameters
mcmc_unique_barseq <- mcmc_run_post(mcmc_run_all_output = mcmc_all_barseq,
                                    Z = opt.clust0.barseq,
                                    thinning = 5,
                                    burn_in = 1500,
                                    number_iter = 3000,
                                    Y = data_barseq,
                                    a_gamma = 500,
                                    b_gamma = 10,
                                    regions.name = rownames(data_barseq[[1]]))
                             

# Number of neurons in each cluster
opt.clustering.frequency(clustering = mcmc_unique_barseq$Z)


# Plot of estimated projection strength
mcmc_unique_barseq$estimated.pp.plot




# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc_barseq <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_barseq,
                                          mcmc_run_post_output = mcmc_unique_barseq,
                                          thinning = 5,
                                          burn_in = 1500,
                                          number_iter = 3000)
                                   


omega_JM_mcmc_barseq$omega_JM_plot


# Difference in omega_{j,m}
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc_barseq)
difference_omega_JM$probability_plot


# Heat-map of projection strengths of each neuron
pp.standard.ordering(Y = data_barseq,
                     Z = mcmc_unique_barseq$Z,
                     regions.name = rownames(data_barseq[[1]]))



# Posterior predictive checks
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_barseq,
                      Y = data_barseq,
                      N = 3,
                      regions.name = rownames(data_barseq[[1]]))

ppc_multiple$zero.plot
ppc_multiple$non.zero.plot


# Distribution of N_{i,m} for neurons within the same cluster
data_barseq_N <- lapply(1:length(data_barseq),
                        function(m) colSums(data_barseq[[m]]))

df <- data.frame(N = unlist(data_barseq_N),
                 motif = unlist(mcmc_unique_barseq$Z))

ggplot(df)+
  geom_boxplot(mapping = aes(x = factor(motif, levels = 1:max(unlist(mcmc_unique_barseq$Z))),
                             y = N))+
  theme_bw()+
  xlab('Cluster')

#-------------------------------------------------------------------------------------------------

# K-means clustering
data_barseq_cbind <- do.call(cbind, data_barseq)

df <- scale(t(data_barseq_cbind))

total.withinss <- data.frame(k = 1:50,
                             total.ss = sapply(1:50,
                                               function(i){
                                                 
                                                 
                                                 km.res <- kmeans(df, i, nstart = 25)
                                                 km.res$tot.withinss
                                               }))

ggplot(total.withinss,
       mapping = aes(x = k, y = total.ss))+
  geom_line(color = 'blue')+
  geom_point(color = 'blue')+
  xlab('Total number of clusters')+
  ylab('Total within sum of square')+
  theme_bw()


