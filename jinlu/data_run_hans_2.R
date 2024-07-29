
M <- length(data_Hans)
C <- sapply(1:M, function(m) ncol(data_Hans[[m]]))

# Add threshold
# data_Hans <- lapply(1:M,
#                       function(m) apply(data_Hans[[m]], c(1,2), function(x) ifelse(x >= 5, x, 0)))

gel_plot_data_Hans <- gel_plot(Y = data_Hans)

png(file = './plots/Hans/gelplot_hans.png',
    width = 2000,
    height = 500)

ggarrange(gel_plot_data_Hans[[1]],
          gel_plot_data_Hans[[2]],
          gel_plot_data_Hans[[3]],
          gel_plot_data_Hans[[4]],
          nrow = 1,
          widths = c(1,1,1,1.3))

dev.off()


mouse.index <- c(rep(1, C[1]),
                 rep(2, C[2]),
                 rep(3, C[3]),
                 rep(4, C[4]))


#Initialize z
data_Hans_cbind <- do.call(cbind, data_Hans)

df <- t(data_Hans_cbind)
df = t(apply(df, 1, function(x){return(x/sum(x))}))


C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(data_Hans[[m]]))))

# k-means
k_mean_clust_20 <- kmeans(df, 20, iter.max = 100, nstart = 25)$cluster

clust20 <- lapply(1:M,
                  function(m) k_mean_clust_20[(C_cumsum[m]+1):C_cumsum[m+1]])

# Run without threshold
mcmc_all_hans <- mcmc_run_all(Y = data_Hans,
                              J = 40,
                              number_iter = 20000,
                              thinning = 5,
                              burn_in = 5000,
                              adaptive_prop = 0.0001,
                              print_Z = TRUE,
                              a_gamma = 20,
                              b_gamma = 1,
                              a_alpha = 1/5,
                              b_alpha = 1/2,
                              Z.init = clust20)


                                

