## Load data
load("~/Documents/GitHub/MAPseq_analysis/data/EC8_new/mcmc_all_EC8.RData")

library(ggplot2)
library(tidyr)

# Compute the total variation distance between the mixing measures for each mouse pair

mytv_dist = function(x,ind){
  xdim = dim(mcmc_all_EC8$omega_J_M_output[[1]])
  y = matrix(x[,ind],xdim[1],xdim[2])
  return(0.5*colSums(abs(x-y)))
}

m =  6 # index of mouse
tv_dist = lapply(mcmc_all_EC8$omega_J_M_output,mytv_dist, ind = m )
tv_dist <- data.frame(matrix(unlist(tv_dist), nrow=length(tv_dist), byrow=TRUE))
names(tv_dist) = c("1", "2", "3", "4", "5", "6")
tv_dist = tv_dist[-m]

tv_dist = tv_dist %>% 
  pivot_longer(
    cols = everything(),
    names_to = "Mouse", 
    values_to = "TV"
  )


group.colors <- c('1'= "#FF3000",'2'= "#FF9900",'3'= "#FFDB6D",'4'= "#009999",'5'= "#333BFF",'6'= "#9633FF")

ggplot(tv_dist) +
  geom_density(aes(x=TV,col=Mouse, fill=Mouse),, alpha=0.25) +
  ggtitle(paste('Mouse',m,'comparisons')) + 
  theme_bw()+
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors)

# Compute the posterior mean
tv_mean = matrix(0, 6,6)
for (m in c(1:6)){
  tv_dist_m = lapply(mcmc_all_EC8$omega_J_M_output,mytv_dist, ind = m )
  tv_dist_m = data.frame(matrix(unlist(tv_dist_m), nrow=length(tv_dist_m), byrow=TRUE))
  tv_mean[m,] = colMeans(tv_dist_m)
}
tv_mean = data.frame(tv_mean, row.names = c("1", "2", "3", "4", "5", "6"))
names(tv_mean) = c("1", "2", "3", "4", "5", "6")
print(tv_mean)

tv_mean = data.frame(tv_mean, "Mouse 1" = c("1", "2", "3", "4", "5", "6"))
names(tv_mean)[1:6] = c("1", "2", "3", "4", "5", "6")
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
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.4)

