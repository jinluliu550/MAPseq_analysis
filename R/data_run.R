# Run the model on the MAPseq data

# pre-processing
LEC <- read.csv('./data/lec_newSet_800_plusBrain.csv')
MEC <- read.csv('./data/mec_newSet_800_plusBrain.csv')


digit.list <- list(one = 1,
                   two = 2,
                   three = 3,
                   four = 4,
                   five = 5,
                   six = 6)

# Change labels from text to numbers
LEC$brain <- sapply(1:nrow(LEC),
                    function(i) which(names(digit.list) == LEC$brain[i]))

MEC$brain <- sapply(1:nrow(MEC),
                    function(i) which(names(digit.list) == MEC$brain[i]))

# round data
data <- round(rbind(LEC, MEC), 0)

data_by_mouse <- lapply(1:6,
                        function(m) data %>%
                          filter(brain == m) %>%
                          select(-brain) %>%
                          as.matrix())

data_by_mouse <- lapply(1:6,
                        function(m) t(data_by_mouse[[m]]))


# Label of LEC and MEC
LEC$EC_label <- 'LEC'
MEC$EC_label <- 'MEC'

data2 <- rbind(LEC,MEC)
data_by_mouse2 <- lapply(1:6,
                         function(m) data2 %>%
                           filter(brain == m) %>%
                           select(EC_label) %>%
                           pull())


# Load everything
load('data/mcmc_all_sample.RData')
load('data/psm.RData')
load('data/opt.clust0.RData')
load('data/mcmc_unique.RData')



# MCMC run
mcmc_all_sample <- mcmc_run_all(Y = data_by_mouse,
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


psm <- similarity_matrix(mcmc_run_all_output = mcmc_all_sample,
                         num.cores = 10,
                         run.on.pc = FALSE)

# Plot of posterior similarity matrix
plotpsm(psm.ind = psm$psm.within,
        psm.tot = psm$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)


plotpsm(psm.ind = psm$psm.within,
        psm.tot = psm$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')


# 118 clusters
opt.clust0 <- opt.clustering(mcmc_run_all_output = mcmc_all_sample,
                             post_similarity = psm)


# MCMC of unique parameters
mcmc_unique <- mcmc_run_post(mcmc_run_all_output = mcmc_all_sample,
                             Z = opt.clust0,
                             thinning = 5,
                             burn_in = 1500,
                             number_iter = 3000,
                             Y = data_by_mouse,
                             a_gamma = 50,
                             b_gamma = 1,
                             front.regions = c(1,2),
                             middle.regions = c(3,4,5),
                             back.regions = c(6,7,8),
                             regions.name = rownames(data_by_mouse[[1]]))

# Number of neurons in each cluster
opt.clustering.frequency(clustering = mcmc_unique$Z)

# Plot of estimated projection strength
mcmc_unique$estimated.pp.plot

mcmc_unique$plot.by.region$front

# Plot of q tilde
mcmc_unique$q_tilde_plot

# unique to one mouse
unique.to.one.mouse <- sapply(1:118,
                              function(j){
                                
                                length(which(sapply(1:6, 
                                                    function(m) length(which(mcmc_unique$Z[[m]]==j))) != 0))
                              })

length(which(unique.to.one.mouse==1))

# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_sample,
                                   mcmc_run_post_output = mcmc_unique,
                                   thinning = 5,
                                   burn_in = 1500,
                                   number_iter = 3000)


omega_JM_mcmc$omega_JM_plot


# Examining the difference of omega_jm between any 2 mice:
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc)

difference_omega_JM$probability_plot

difference_in_omega_jm_plot(difference_in_omega_jm_output = difference_omega_JM,
                            mcmc_run_omega_output = omega_JM_mcmc,
                            N = 20)

sapply(1:M, 
       function(m){
         
         length(which((difference_omega_JM$significant_obs$data1 == m)|
                        difference_omega_JM$significant_obs$data2 == m))
       })


# Heat-map of projection strengths of each neuron
pp.reordered(qj_estimate = mcmc_unique$proj_prob_mean,
             Y = data_by_mouse,
             Z = mcmc_unique$Z,
             regions.name = rownames(data_by_mouse[[1]]))

# Posterior predictive checks
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_sample,
                      Y = data_by_mouse,
                      N = 3,
                      regions.name = rownames(data_by_mouse[[1]]))

ppc_multiple$zero.plot
ppc_multiple$non.zero.plot

ppc_single_f(mcmc_run_all_output = mcmc_all_sample,
             Y = data_by_mouse)
