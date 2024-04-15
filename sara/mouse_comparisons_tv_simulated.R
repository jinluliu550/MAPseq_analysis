## Load data
#load("~/Documents/GitHub/MAPseq_analysis/data/EC8_new/mcmc_all_EC8.RData")

#mcmc_results = mcmc_all_EC8
mcmc_results = mcmc_simulated

library(ggplot2)
library(tidyr)

# Compute the total variation distance between the mixing measures for each mouse pair

mytv_dist = function(x,ind){
  xdim = dim(mcmc_results$omega_J_M_output[[1]])
  y = matrix(x[,ind],xdim[1],xdim[2])
  return(0.5*colSums(abs(x-y)))
}

m =  3 # index of mouse
tv_dist = lapply(mcmc_results$omega_J_M_output,mytv_dist, ind = m )
tv_dist <- data.frame(matrix(unlist(tv_dist), nrow=length(tv_dist), byrow=TRUE))
#names(tv_dist) = c("1", "2", "3", "4", "5", "6")
names(tv_dist) = as.factor(c(1:mcmc_results$M))
tv_dist = tv_dist[-m]

tv_dist = tv_dist %>% 
  pivot_longer(
    cols = everything(),
    names_to = "Mouse", 
    values_to = "TV"
  )


#group.colors <- c('1'= "#FF3000",'2'= "#FF9900",'3'= "#FFDB6D",'4'= "#009999",'5'= "#333BFF",'6'= "#9633FF")
group.colors <- c('1'= "#FF3000",'2'= "#009999",'3'= "#333BFF",'4'= "#9633FF")

ggplot(tv_dist) +
  geom_density(aes(x=TV,col=Mouse, fill=Mouse), alpha=0.25) +
  ggtitle(paste('Mouse',m,'comparisons')) + 
  theme_bw()+
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors)


# Compute the posterior mean
tv_mean = matrix(0, mcmc_results$M,mcmc_results$M)
for (m in c(1:mcmc_results$M)){
  tv_dist_m = lapply(mcmc_results$omega_J_M_output,mytv_dist, ind = m )
  tv_dist_m = data.frame(matrix(unlist(tv_dist_m), nrow=length(tv_dist_m), byrow=TRUE))
  tv_mean[m,] = colMeans(tv_dist_m)
}
tv_mean = data.frame(tv_mean, row.names = as.factor(c(1:mcmc_results$M)))
names(tv_mean) = as.factor(c(1:mcmc_results$M))
print(tv_mean)

tv_mean = data.frame(tv_mean, "Mouse 1" = as.factor(c(1:mcmc_results$M)))
names(tv_mean)[1:mcmc_results$M] = as.factor(c(1:mcmc_results$M))
tv_mean <-  pivot_longer(tv_mean,
  cols = !Mouse.1,
  names_to = "Mouse.2", 
  values_to = "TV"
)

ggplot(tv_mean, aes(x = Mouse.1, y = Mouse.2, fill = TV)) +
  geom_tile() +
  labs(x = "Mouse",
       y = "Mouse") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.4) +
  geom_text(aes(label = round(tv_mean$TV,3)), color = "black", size = 4) +
  coord_fixed()

##### Compute HPD interval
library(coda)
# Compute HPD intervals
tv_lower = matrix(0, mcmc_results$M,mcmc_results$M)
tv_upper = matrix(0, mcmc_results$M,mcmc_results$M)
for (m in c(1:mcmc_results$M)){
  tv_dist_m = lapply(mcmc_results$omega_J_M_output,mytv_dist, ind = m )
  tv_dist_m = as.mcmc(matrix(unlist(tv_dist_m), nrow=length(tv_dist_m), byrow=TRUE))
  tv_hpd =  HPDinterval((tv_dist_m))
  tv_lower[m,] = tv_hpd[,1]
  tv_upper[m,] = tv_hpd[,2]
}

tv_lower = data.frame(tv_lower, row.names = as.factor(c(1:mcmc_results$M)))
names(tv_lower) = as.factor(c(1:mcmc_results$M))
print(tv_lower)

tv_lower = data.frame(tv_lower, "Mouse 1" = as.factor(c(1:mcmc_results$M)))
names(tv_lower)[1:mcmc_results$M] = as.factor(c(1:mcmc_results$M))
tv_lower <-  pivot_longer(tv_lower,
                         cols = !Mouse.1,
                         names_to = "Mouse.2", 
                         values_to = "TV"
)

tv_upper = data.frame(tv_upper, row.names = as.factor(c(1:mcmc_results$M)))
names(tv_upper) = as.factor(c(1:mcmc_results$M))
print(tv_upper)

tv_upper = data.frame(tv_upper, "Mouse 1" = as.factor(c(1:mcmc_results$M)))
names(tv_upper)[1:mcmc_results$M] = as.factor(c(1:mcmc_results$M))
tv_upper <-  pivot_longer(tv_upper,
                          cols = !Mouse.1,
                          names_to = "Mouse.2", 
                          values_to = "TV"
)

labs = paste0(round(tv_mean$TV,3),' [',round(tv_lower$TV,3),',', round(tv_upper$TV,3),']')
ggplot(tv_mean, aes(x = Mouse.1, y = Mouse.2, fill = TV)) +
  geom_tile() +
  labs(x = "Mouse",
       y = "Mouse") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.4) +
  geom_text(aes(label = labs), color = "black", size = 4) +
  coord_fixed()
