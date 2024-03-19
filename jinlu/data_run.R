# Run the model on the MAPseq data

# pre-processing
LEC <- read.csv('./data/EC/lec_newSet_800_plusBrain.csv')
MEC <- read.csv('./data/EC/mec_newSet_800_plusBrain.csv')

R <- 9
M <- 6

# Load data
load("./data/EC/mcmc_all_sample.RData")
load("./data/EC/psm_EC.RData")
load("./data/EC/opt.clust0.EC.RData")
load("./data/EC/mcmc_unique_EC.RData")

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


data_by_mouse <- lapply(1:M,
                        function(m) data %>%
                          filter(brain == m) %>%
                          select(-brain) %>%
                          as.matrix())

data_by_mouse <- lapply(1:M,
                        function(m) t(data_by_mouse[[m]]))




# Label of LEC and MEC of each mouse
LEC$EC_label <- 'LEC'
MEC$EC_label <- 'MEC'

data_by_mouse2 <- lapply(1:M,
                         function(m) rbind(LEC,MEC) %>%
                           filter(brain == m) %>%
                           select(EC_label) %>%
                           pull())

# A summary data frame
data3 <- data %>%
  mutate(neuron = c(paste('LEC neuron', 1:nrow(LEC)),
                    paste('MEC neuron', 1:nrow(MEC))))

data3_rowsums <- rowSums(data3[,1:R])

for(i in 1:R){
  
  data3[,i] <- data3[,i]/data3_rowsums
}


neuron_in_each_mouse <- lapply(1:M, 
                               function(m) data3 %>%
                                 filter(brain == m) %>%
                                 select(neuron) %>%
                                 pull())


# Allocations
allocation_bayesian <- data.frame(neuron = unlist(neuron_in_each_mouse),
                                  allocation_bayes = unlist(mcmc_unique$Z))




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


psm_EC <- similarity_matrix(mcmc_run_all_output = mcmc_all_sample,
                         num.cores = 10,
                         run.on.pc = FALSE)

# Plot of posterior similarity matrix
plotpsm(psm.ind = psm_EC$psm.within,
        psm.tot = psm_EC$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)


plotpsm(psm.ind = psm_EC$psm.within,
        psm.tot = psm_EC$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')

# Optimal clustering
opt.clust0.EC <- dlso_cluster_estimate(mcmc_run_all_output = mcmc_all_sample)



# MCMC of unique parameters
mcmc_unique_EC <- mcmc_run_post(mcmc_run_all_output = mcmc_all_sample,
                             Z = opt.clust0.EC,
                             thinning = 5,
                             burn_in = 1500,
                             number_iter = 3000,
                             Y = data_by_mouse,
                             a_gamma = 500,
                             b_gamma = 10,
                             regions.name = rownames(data_by_mouse[[1]]))

# Number of neurons in each cluster - 148 clusters
png(file = './plots/EC/number_of_neuron.png',
    width = 2500,
    height = 700)

opt.clustering.frequency(clustering = mcmc_unique_EC$Z)

dev.off()

# Number of LEC and MEC neurons in each cluster
png(file = './plots/EC/number_of_neuron_by_EC.png',
    width = 2500,
    height = 700)

opt.clustering.frequency1(clustering = mcmc_unique_EC$Z,
                          EC_label = data_by_mouse2)

dev.off()

# Proportion of LEC and MEC in each cluster
png(file = './plots/EC/proportion_of_EC.png',
    width = 2500,
    height = 700)

opt.clustering.frequency2(clustering = mcmc_unique_EC$Z,
                          EC_label = data_by_mouse2)

dev.off()

# Plot of estimated projection strength
png(file = './plots/EC/estimated_pp.png',
    width = 3500,
    height = 2000)

mcmc_unique_EC$estimated.pp.plot

dev.off()

# Plot of q tilde
png(file = './plots/EC/q_tilde.png',
    width = 3000,
    height = 600)

mcmc_unique_EC$q_tilde_plot

dev.off()

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

# Projection strength of each neuron in each cluster, colour-coded by the injection site
projection_by_EC(Y = data_by_mouse,
                 EC_label = data_by_mouse2,
                 Z = mcmc_unique$Z,
                 region_name = rownames(data_by_mouse[[1]]))

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
pp.standard.ordering(Y = data_by_mouse,
                     Z = mcmc_unique$Z,
                     regions.name = rownames(data_by_mouse[[1]]))


# Distribution of N_{i,m} for neurons within the same cluster
data_EC_N <- lapply(1:length(data_by_mouse),
                        function(m) colSums(data_by_mouse[[m]]))

df <- data.frame(N = unlist(data_EC_N),
                 motif = unlist(mcmc_unique$Z))

ggplot(df)+
  geom_boxplot(mapping = aes(x = factor(motif, levels = 1:max(unlist(mcmc_unique$Z))),
                             y = N))+
  theme_bw()+
  xlab('Cluster')

# Posterior predictive checks
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_sample,
                      Y = data_by_mouse,
                      N = 3,
                      regions.name = rownames(data_by_mouse[[1]]))

ppc_multiple$zero.plot
ppc_multiple$non.zero.plot

ppc_single_f(mcmc_run_all_output = mcmc_all_sample,
             Y = data_by_mouse)
