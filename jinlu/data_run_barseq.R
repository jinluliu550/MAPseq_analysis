# Running on the barSeq data


# Load data
data1 <- read.csv('./data/barseq/BRAIN_C9.csv')
data2 <- read.csv('./data/barseq/BRAIN_C14.csv')
data3 <- read.csv('./data/barseq/BRAIN_C28.csv')


# Round to integers
data1 <- round(data1, 0)
data2 <- round(data2, 0)
data3 <- round(data3, 0)

data_barseq <- list(t(data1),
                    t(data2),
                    t(data3))

M <- length(data_barseq)
C <- sapply(1:M, function(m) ncol(data_barseq[[m]]))


#Initialize z
data_barseq_cbind <- do.call(cbind, data_barseq)

df <- t(data_barseq_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))


C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(data_barseq[[m]]))))

# k-means
k_mean_clust_60 <- kmeans(df, 60, iter.max = 100, nstart = 25)$cluster

clust60 <- lapply(1:M,
                  function(m) k_mean_clust_60[(C_cumsum[m]+1):C_cumsum[m+1]])



# MCMC run
mcmc_all_barseq <- mcmc_run_all(Y = data_barseq,
                                J = 80,
                                number_iter = 20000,
                                thinning = 5,
                                burn_in = 5000,
                                adaptive_prop = 0.0001,
                                print_Z = TRUE,
                                a_gamma = 30,
                                b_gamma = 1,
                                a_alpha = 1/5,
                                b_alpha = 1/2,
                                Z.init = clust60)


Zmat = matrix(unlist(mcmc_all_barseq$Z_output), length(mcmc_all_barseq$Z_output), sum(C),byrow = TRUE)

# Number of occupied components
k = apply(Zmat,1,function(x){length(unique(x))})
plot(k, type = 'l')


# Posterior similarity matrix
psm_barseq = similarity_matrix(mcmc_run_all_output = mcmc_all_barseq)


# Reordered posterior samples of z
barseq_z_reordered <- z_trace_updated(mcmc_run_all_output = mcmc_all_barseq)


# optimal clustering
barseq_Z <- opt.clustering.comb(z_trace = barseq_z_reordered,
                                post_similarity = psm_barseq,
                                max.k = max(k))


#-- Convert to a list
C_cumsum <- c(0, cumsum(C))
barseq_Z <- lapply(1:M,
                 function(m) barseq_Z[(C_cumsum[m]+1):C_cumsum[m+1]])


# Plot of posterior similarity matrix
png(file = './plots/barseq/heatmap_psm_1.png',
    width = 664,
    height = 664)

plotpsm(psm.ind = psm_barseq$psm.within,
        psm.tot = psm_barseq$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)

dev.off()


png(file = './plots/barseq/heatmap_psm_2.png',
    width = 664,
    height = 664)

plotpsm(psm.ind = psm_barseq$psm.within,
        psm.tot = psm_barseq$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')

dev.off()

# MCMC unique
mcmc_unique_barseq <- mcmc_run_post(mcmc_run_all_output = mcmc_all_barseq,
                                  Z = barseq_Z,
                                  thinning = 5,
                                  burn_in = 2000,
                                  number_iter = 12000,
                                  Y = data_barseq,
                                  a_gamma = 30,
                                  b_gamma = 1,
                                  regions.name = rownames(data_barseq[[1]]))

barseq_Z_reordered <- mcmc_unique_barseq$Z


# Number of neurons by cluster and mosue
png(file = './plots/barseq/number_of_neuron_by_m.png',
    width = 2500,
    height = 700)

opt.clustering.frequency(clustering = mcmc_unique_barseq$Z)

dev.off()


# Plot of estimated projection strength
png(file = './plots/barseq/estimated_pp.png',
    width = 3500,
    height = 2000)

mcmc_unique_barseq$estimated.pp.plot

dev.off()


# q tilde
png(file = './plots/barseq/q_tilde.png',
    width = 3000,
    height = 600)

mcmc_unique_barseq$q_tilde_plot

dev.off()


# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_barseq,
                                   mcmc_run_post_output = mcmc_unique_barseq,
                                   thinning = 5,
                                   burn_in = 2000,
                                   number_iter = 12000)



png(file = './plots/barseq/w_jm.png',
    width = 3000,
    height = 1500)

omega_JM_mcmc$omega_JM_plot

dev.off()

# Examining the difference of omega_jm between any 2 mice:
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc)


png(file = './plots/barseq/w_jm_difference.png',
    width = 1000,
    height = 400)

difference_omega_JM$probability_plot

dev.off()


png(file = './plots/barseq/w_jm_difference_eg.png',
    width = 1000,
    height = 400)

difference_in_omega_jm_plot(difference_in_omega_jm_output = difference_omega_JM,
                            mcmc_run_omega_output = omega_JM_mcmc,
                            N = 10)

dev.off()


# Heatmap of projection stength of neurons in each cluster
png(file = './plots/barseq/heatmap_neuron.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = data_barseq,
                     Z = mcmc_unique_barseq$Z,
                     regions.name = rownames(data_barseq[[1]]))

dev.off()


# Distribution of N_{i,m} for neurons within the same cluster
data_barseq_N <- lapply(1:length(data_barseq),
                      function(m) colSums(data_barseq[[m]]))

df <- data.frame(N = unlist(data_barseq_N),
                 motif = unlist(mcmc_unique_barseq$Z))

png(file = './plots/barseq/N_sum_by_cluster.png',
    width = 1200,
    height = 600)

ggplot(df)+
  geom_boxplot(mapping = aes(x = factor(motif, levels = 1:max(unlist(mcmc_unique_barseq$Z))),
                             y = N))+
  theme_bw()+
  xlab('Cluster')

dev.off()


# Posterior predictive check with multiple replicates
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_barseq,
                      Y = data_barseq,
                      N = 3,
                      regions.name = rownames(data_barseq[[1]]))


png(file = './plots/barseq/ppc_zero.png',
    width = 1200,
    height = 700)

ppc_multiple$zero.plot

dev.off()


png(file = './plots/barseq/ppc_nonzero.png',
    width = 1200,
    height = 700)

ppc_multiple$non.zero.plot

dev.off()


# Posterior predictive check with single replicated data
ppc_single <- ppc_single_f(mcmc_run_all_output = mcmc_all_barseq,
                           Y = data_barseq,
                           regions.name = rownames(data_barseq[[1]]))

for(m in 1:M){
  
  png(file = paste0('./plots/barseq/ppc_single_mouse_', m, '.png'),
      width = 1200,
      height = 700)
  
  print(ppc_single[[m]])
  
  dev.off()
}


#Initialize z
data_Hans_cbind <- do.call(cbind, data_barseq)

df <- t(data_Hans_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))

k_mean_clust_40 <- kmeans(df, 40, nstart = 25)$cluster

clust40 <- lapply(1:M,
                  function(m) k_mean_clust_40[(C_cumsum[m]+1):C_cumsum[m+1]])

clust40_r <- k_means_reorder(Y = data_barseq,
                             Z = clust40)

k_mean_clust_60 <- kmeans(df, 60, nstart = 25)$cluster

clust60 <- lapply(1:M,
                  function(m) k_mean_clust_60[(C_cumsum[m]+1):C_cumsum[m+1]])


clust60_r <- k_means_reorder(Y = data_barseq,
                             Z = clust60)


k_mean_clust_80 <- kmeans(df, 80, nstart = 25)$cluster

clust80 <- lapply(1:M,
                  function(m) k_mean_clust_80[(C_cumsum[m]+1):C_cumsum[m+1]])


clust80_r <- k_means_reorder(Y = data_barseq,
                             Z = clust80)


# Heatmap of projection strength of neurons in each cluster
png(file = './plots/barseq/k_means_40.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = data_barseq,
                     Z = clust40_r,
                     regions.name = rownames(data_barseq[[1]]))

dev.off()

png(file = './plots/barseq/k_means_60.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = data_barseq,
                     Z = clust60_r,
                     regions.name = rownames(data_barseq[[1]]))

dev.off()

png(file = './plots/barseq/k_means_80.png',
    width = 900,
    height = 900)

pp.standard.ordering(Y = data_barseq,
                     Z = clust80_r,
                     regions.name = rownames(data_barseq[[1]]))

dev.off()