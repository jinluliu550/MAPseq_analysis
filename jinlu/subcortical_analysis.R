all_brains_set_sub <- read_excel("data/EC8_new/all_brains_set_sub.xlsx")


# Convert to list
digit.list <- list(one = 1,
                   two = 2,
                   three = 3,
                   four = 4,
                   five = 5,
                   six = 6)

EC8_subnew_list <- lapply(1:6, 
                       function(m) all_brains_set_sub %>%
                         filter(brain == names(digit.list)[m]))

EC8_sub_new <- lapply(1:6,
                  function(m) t(round(EC8_subnew_list[[m]][,c(1:8, 11, 10, 9)], 0)))






df_sub <- t(do.call(cbind, EC8_sub_new))
df_sub = t(apply(df_sub, 1, function(x){return(x/sum(x))}))

C <- sapply(1:6, function(m) ncol(EC8_sub_new[[m]]))

C_cumsum <- c(0, cumsum(sapply(1:6, function(m) ncol(EC8_sub_new[[m]]))))

mouse.index <- c(rep(1, C[1]),
                 rep(2, C[2]),
                 rep(3, C[3]),
                 rep(4, C[4]),
                 rep(5, C[5]),
                 rep(6, C[6]))

# K-means
k_mean_clust_60 <- kmeans(df_sub, 60, nstart = 25)$cluster

clust60 <- lapply(1:6,
                  function(m) k_mean_clust_60[(C_cumsum[m]+1):C_cumsum[m+1]])


clust60_r <- k_means_reorder(Y = EC8_sub_new,
                             Z = clust60)

# Heatmap of neurons
pp.standard.ordering2(Y = EC8_sub_new,
                      Z = clust60_r,
                      mouse.index = mouse.index)

# Run MCMC
mcmc_all_EC8_subcortical <- mcmc_run_all(Y = EC8_sub_new,
                                         J = 100,
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
                            




Zmat = matrix(unlist(mcmc_all_EC8_subcortical$Z_output), length(mcmc_all_EC8_subcortical$Z_output), sum(C),byrow = TRUE)

# Number of occupied components
k = apply(Zmat,1,function(x){length(unique(x))})
plot(k, type = 'l')

# Posterior similarity matrix
psm_EC8 = similarity_matrix(mcmc_run_all_output = mcmc_all_EC8_subcortical)

# Reordered posterior samples of z
EC8_z_reordered <- z_trace_updated(mcmc_run_all_output = mcmc_all_EC8_subcortical)

# optimal clustering
EC8_Z <- opt.clustering.comb(z_trace = EC8_z_reordered,
                             post_similarity = psm_EC8,
                             max.k = max(k))

#-- Convert to a list
C_cumsum <- c(0, cumsum(C))
EC8_Z <- lapply(1:6,
                function(m) EC8_Z[(C_cumsum[m]+1):C_cumsum[m+1]])

length(unique(unlist(EC8_Z)))


# Plot of posterior similarity matrix
psm_plot <- plotpsm(psm.ind = psm_EC8$psm.within,
                    psm.tot = psm_EC8$psm.combined,
                    method = 'complete')

png(file = './plots/EC8_sub/heatmap_psm_1.png',
    width = 664,
    height = 664)

psm_plot$plot.ind

dev.off()


png(file = './plots/EC8_sub/heatmap_psm_2.png',
    width = 664,
    height = 664)

psm_plot$plot.tot

dev.off()


# MCMC unique
mcmc_unique_EC8_sub <- mcmc_run_post(mcmc_run_all_output = mcmc_all_EC8_subcortical,
                                 Z = EC8_Z,
                                 thinning = 5,
                                 burn_in = 2000,
                                 number_iter = 12000,
                                 Y = EC8_sub_new,
                                 a_gamma = 30,
                                 b_gamma = 1,
                                 regions.name = rownames(EC8_sub_new[[1]]))

EC8_Z_reordered <- mcmc_unique_EC8_sub$Z

# Heat-map

png(file = './plots/EC8_sub/heatmap_neuron.png',
    width = 900,
    height = 900)

pp.standard.ordering2(Y = EC8_sub_new,
                      Z = mcmc_unique_EC8_sub$Z,
                      mouse.index = mouse.index,
                      regions.name = rownames(EC8_sub_new[[1]]))

dev.off()


mcmc_unique_EC8_sub$estimated.pp.plot
