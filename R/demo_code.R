# Demonstration

# The input data must be a list, each item within the list is a matrix. Each matrix must have equal number of rows,
# but could have different number of columns. For example in the brain connectivity case, each matrix has 8 rows,
# one row representing each brain region. Each column of the matrix represents a specific neuron.

# To obtain MCMC samples of all parameters of interest: We run the following code. In the following run, we assume
# there are 3 components, and a_gamma = 50, b_gamma = 1, a_alpha = 1, b_alpha = 1, 1500 iterations of burn-in, thinning = 5,
# and 3000 total number of iterations. In addition, the function will run on a single core and print out the summary of allocation
# after each iteration.


mcmc_all_sample <- mcmc_run_all(Y = Y,
                                J = 3,
                                number_iter = 3000,
                                thinning = 5,
                                burn_in = 1000,
                                adaptive_prop = 0.1,
                                print_Z = TRUE,

                                # For the example data, this combination of a_gamma and b_gamma seems to work really well.
                                # Those values need to be tuned for different datasets.
                                a_gamma = 50,
                                b_gamma = 1,
                                a_alpha = 1,
                                b_alpha = 1,
                                num.cores = 1)




# Acceptance probabilities
mcmc_all_sample$acceptance_prob

# List of allocation probability of each neuron in each iteration
mcmc_all_sample$allocation_probability

# List of MCMC samples
mcmc_all_sample$Z_output
mcmc_all_sample$omega_J_M_output
mcmc_all_sample$omega_output
mcmc_all_sample$alpha_output
mcmc_all_sample$alpha_zero_output
mcmc_all_sample$alpha_h_output
mcmc_all_sample$q_star_1_J_output
mcmc_all_sample$gamma_star_1_J_output

#------------------------------------------------------------------------------------------------------------------

# To obtain posterior similarity matrix. we run the following code. The function returns a list;
# The first item shows the posterior similarity of neurons for the same mouse, or between a pair of two mice.
#' For example: item [[a]][[b]] is the posterior similarity between neurons of mouse a and b. Note that to avoid duplication,
#' index a must be smaller or equal to b. The second item shows the posterior similarity of neurons from all available mice.

psm <- similarity_matrix(mcmc_run_all_output = mcmc_all_sample)

# To have the posterior similarity plot with neurons from difference mice separated, we run the following code.
# The diagonals shows the neuron similarity for each mouse, and the off-diagonals shows the neuron similarity
# between difference mice. hclust() is used to shows clustering property.

plotpsm(psm.ind = psm$psm.within,
        psm.tot = psm$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        plot.type = 'ind',
        cex.main = 1.5)

# Another version of the posterior similarity plot is to ignore the source of the neuron, i.e.: neurons from all mouse are
# put together.

plotpsm(psm.ind = psm$psm.within,
        psm.tot = psm$psm.combined,
        xlab = 'neurons',
        ylab = 'neurons',
        cex.lab = 1.5,
        main = 'Posterior Similarity Matrix',
        cex.main = 1.5,
        plot.type = 'tot')

# To obtain the clustering estimate which minimizes the variation of information, we have

opt.clust0 <- opt.clustering(mcmc_run_all_output = mcmc_all_sample,
                             post_similarity = psm)

# The optimal clustering here is a list of length M, and each item of the list has length
# equal to the number of neurons in m[th] mouse.



#------------------------------------------------------------------------------------------------------------------


# To avoid label switching problem, we need to obtain point estimate of the projection strength and gamma
# using the following function.

mcmc_unique <- mcmc_run_post(mcmc_run_all_output = mcmc_all_sample,
                             Z = opt.clust0,
                             thinning = 5,
                             burn_in = 1500,
                             number_iter = 3000,
                             Y = Y,
                             a_gamma = 50,
                             b_gamma = 1,
                             front.regions = 1,
                             middle.regions = c(2,3))


# This function has the following outputs:

# MCMC samples of projection strengths
mcmc_unique$q_star_1_J_output

# Mean, lower and upper limit of the 90 percent credible interval
mcmc_unique$proj_prob_mean
mcmc_unique$proj_prob_lower
mcmc_unique$proj_prob_upper

# Summary of estimated projection strengths in a data frame
mcmc_unique$estimated.projection.df

# Mean and trace of gamma
mcmc_unique$gamma_mean
mcmc_unique$gamma_star_1_J_output

# Summary of clusters (classify the cluster by the number and location of projecting regions)
mcmc_unique$cluster.label.summary


# Plot of estimated projection strengths - all clusters
# Those are color-coded by the number of projecting regions
mcmc_unique$estimated.pp.plot

# In the case when the location of brain regions are given
mcmc_unique$plot.by.region


# Plot for the probability that the projection strength is greater
# than epsilon for each cluster and region
mcmc_unique$q_tilde_plot

# Table of the probability that the projection strength is greater than
# epsilon for each cluster and region
mcmc_unique$q_tilde_001

# Re-labelled clusters - Note we will be using this cluster from here onward;
# Clusters are reordered depending on the location of the projecting region
opt.clust <- mcmc_unique$Z

# Bar chart to summarize the number of clusters and the number of neurons
# from each mouse in each cluster
opt.clustering.frequency(clustering = opt.clust)

#------------------------------------------------------------------------------------

# Function below obtains posterior samples of mouse-specific component probabilities
omega_JM_mcmc <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_all_sample,
                                   mcmc_run_post_output = mcmc_unique,
                                   thinning = 5,
                                   burn_in = 1500,
                                   number_iter = 3000)


# Posterior samples
omega_JM_mcmc$omega_J_M_output

# Posterior mean, lower and upper bound of 90 percent credible interval
omega_JM_mcmc$mean_omega_J_M
omega_JM_mcmc$lower_omega_J_M
omega_JM_mcmc$upper_omega_J_M

# Summary of posterior estimated mouse-specific component probabilities in a data frame
omega_JM_mcmc$omega_JM_estimate_df

# Box plot of estimated values for all clusters
omega_JM_mcmc$omega_JM_plot

# Box plot of estimated values with clusters separated depending on the location of projecting regions.
omega_JM_mcmc$plot.by.region


# Examining the difference of omega_jm between any 2 mice:
difference_omega_JM <- difference_in_omega_jm(mcmc_run_omega_output = omega_JM_mcmc)

# This function has the following outputs

# Difference between all possible pairs of mice and clusters
difference_omega_JM$all_obs

# Difference between all pairs with significant difference
difference_omega_JM$significant_obs

# A heat-map to show difference of omega_jm between any 2 mice and for each cluster
difference_omega_JM$probability_plot

# To compare the distribution of omega_jm of pair of mice with
# significantly difference in omega_jm, we use the following function to choose
# N examples. The code below will show 2 plots of histograms to compare the
# distribution of omega_jm between 2 mouse pairings.
difference_in_omega_jm_plot(difference_in_omega_jm_output = difference_omega_JM,
                            mcmc_run_omega_output = omega_JM_mcmc,
                            N = 2)

#-------------------------------------------------------------------------------------------

# We can also visualize the projection strength of each neuron by cluster.
# Each row represent a neuron, and each column represent a brain region.
# Neurons from different clusters are separated by horizontal lines.

# We can choose either to not reorder the clusters,
pp.standard.ordering(Y = Y,
                     Z = opt.clust)

# or we could reorder the clusters by selecting the dominant brain
# region (the strongest projecting) region of each cluster, to show a step pattern.
pp.reordered(qj_estimate = mcmc_unique$proj_prob_mean,
             Y = Y,
             Z = opt.clust)

#--------------------------------------------------------------------------------------------

# For posterior predictive checks, we could either use multiple replicated data or single
# replicated data.

# For the case with multiple replicated data, we compare the number of zeros
# and the distribution of non-zero neuron counts in each brain region.
ppc_multiple <- ppc_f(mcmc_run_all_output = mcmc_all_sample,
                      Y = Y,
                      N = 3)

# Projection strength of the observed data
ppc_multiple$Y_prop

# Projection strength of each of the replicated data
ppc_multiple$Y_rep_prop

# Parameters used to simulate each of the replicated data
ppc_multiple$theta

# Histogram to show the total number of zeros
ppc_multiple$zero.plot

# Boxplot to show the distribution of non-zero counts
ppc_multiple$non.zero.plot

# For the case with single replicated data, we can compare the distribution of neuron
# counts for each mouse and brain region using histograms.
ppc_single_f(mcmc_run_all_output = mcmc_all_sample,
             Y = Y)
