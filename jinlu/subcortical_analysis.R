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
                  function(m) t(round(EC8_subnew_list[[m]][,1:11], 0)))



df_sub <- t(do.call(cbind, EC8_sub_new))
df_sub = t(apply(df_sub, 1, function(x){return(x/sum(x))}))


C_cumsum <- c(0, cumsum(sapply(1:M, function(m) ncol(EC8_sub_new[[m]]))))

# K-means
k_mean_clust_70 <- kmeans(df_sub, 70, nstart = 25)$cluster

clust70 <- lapply(1:6,
                  function(m) k_mean_clust_70[(C_cumsum[m]+1):C_cumsum[m+1]])


clust70_r <- k_means_reorder(Y = EC8_sub_new,
                             Z = clust70)


pp.standard.ordering(Y = EC8_sub_new,
                     Z = clust70_r,
                     regions.name = rownames(EC8_sub_new[[1]]))

# Neurons which project to the subcortical regions
EC8_sub_only <- lapply(1:6,
                       function(m) EC8_subnew_list[[m]] %>% filter(OLF > 0 | vStr > 0 | dStr > 0))