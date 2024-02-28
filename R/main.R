
##---------- Load all required packages --------------------

suppressPackageStartupMessages({

  library(mvtnorm)
  library(extraDistr)
  library(SciViews)
  library(doParallel)
  library(foreach)
  library(ggplot2)
  library(tidyverse)
  library(ggpubr)
  library(fields)
  library(Matrix)
  library(Rtsne)
  library(ggpubr)
  library(ggrepel)
  library(parallel)
  library(progress)
  library(pheatmap)
  library(RColorBrewer)
  library(readxl)
  library(truncnorm)
})

`%notin%` <- Negate(`%in%`)


# MCMC
source('R/mcmc_steps.R')
source('R/mcmc_run_all.R')
source('R/mcmc_run_post.R')
source('R/mcmc_run_omega_JM.R')
source('R/cluster_label.R')

# Posterior similarity matrix and plot
source('R/similarity_matrix.R')
source('R/posterior_similarity_plot.R')

# Optimal clustering
source('R/optimal_clustering.R')
source('R/clustering_estimate.R')
source('R/clustering_estimate_plot.R')
source('R/difference_in_omega_jm.R')
source('R/difference_in_omega_jm_plot.R')


# Projection probability by cluster
source('R/pp_by_z_no_reordering.R')
source('R/pp_by_z_reordering.R')

# Posterior predictive checks
source('R/ppc.R')
source('R/ppc_single.R')

