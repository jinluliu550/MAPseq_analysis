## Compute conditional probabilities of projection
library(coda)

# Compute the marginal probability of zero counts for each region for a DirMult mixture
ppzero <- function(n, alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  pp = matrix(0,R,1)
  for (j in 1:J){
    pp = pp+w[j]*apply(matrix(1:R,R,1), 1, function(r){ddirmultinomial(c(0,n),n,c(alpha[j,r],sum(alpha[j,])-alpha[j,r]))})
  }
  return(pp)
}

# Compute the marginal probability of zero counts for each pair of region for a DirMult mixture
ppzero_pair <- function(n, alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  pp = matrix(0,R,R)
  for (j in 1:J){
    for (r1 in 1:R){
      pp[r1,] = pp[r1,]+w[j]*apply(matrix(1:R,R,1), 1, function(r2){ddirmultinomial(c(0,0,n),n,c(alpha[j,r1],alpha[j,r2],sum(alpha[j,])-alpha[j,r1]-alpha[j,r2]))})
    }
  }
  return(pp)
}

# Compute the marginal probability of non-zero counts for each region and the conditional
# proability of non-zero counts for region one given another
projection_prob= function(m,N,mcmc_all_out){
  
  R= mcmc_all_out$R
  I = length(mcmc_all_out$alpha_output)
  
  #output
  pprob_pair_array = array(0,dim= c(I,R,R))
  pprob_mat = matrix(0,I,R)
  
  for(i in 1:I){
     wjm = mcmc_all_out$omega_J_M_output[[i]][,m]
     qstarjm = mcmc_all_out$q_star_1_J_output[[i]]
     gammastarjm = mcmc_all_out$gamma_star_1_J_output[[i]]
     pzero_i = ppzero(N,qstarjm*gammastarjm, wjm)
     pzero_pair_i = ppzero_pair(N,qstarjm*gammastarjm, wjm)
     pprob_mat[i,] = 1- pzero_i
     pprob_pair_temp =  (1-(matrix(pzero_i, R,R) + t(matrix(pzero_i, R,R)) - pzero_pair_i))/matrix(1-pzero_i, R,R)
     diag(pprob_pair_temp) = 1
     pprob_pair_array[i,,] = pprob_pair_temp 
  }
  pprob = apply(pprob_mat,2,mean)
  pprob_pair = apply(pprob_pair_array,c(2,3),mean)
  pprob_ci = apply(pprob_mat,2,function(x){HPDinterval(as.mcmc(x))})
  pprobpair_ci = apply(pprob_pair_array,c(2,3),function(x){HPDinterval(as.mcmc(x))})
  return(list(pprob = pprob, pprob_pair = pprob_pair, pprob_ci = pprob_ci, pprob_pair_ci = pprobpair_ci))
}

# pprob_pair + pprob_pair_temp 

#Mouse 1
pp_m1 = projection_prob(1,100,mcmc_all_EC8)

#Mouse 2
pp_m2 = projection_prob(2,100,mcmc_all_EC8)

#Mouse 3
pp_m3 = projection_prob(3,100,mcmc_all_EC8)

#Mouse 4
pp_m4 = projection_prob(4,100,mcmc_all_EC8)

#Mouse 5
pp_m5 = projection_prob(5,100,mcmc_all_EC8)

#Mouse 6
pp_m6 = projection_prob(6,100,mcmc_all_EC8)

## Plots

#Mouse 1
pp1 = data.frame(pp_m1$pprob_pair, row.names = rownames(EC8_new[[1]]))
names(pp1) = rownames(EC8_new[[1]])
print(pp1)

pp1 = data.frame(pp1, "Region B" = factor(rownames(EC8_new[[1]]), levels = rownames(EC8_new[[1]])))
names(pp1)[1:R] = rownames(EC8_new[[1]])
pp1 <-  pivot_longer(pp1,
                         cols = !Region.B,
                         names_to = "Region.A", 
                         values_to = "CP"
)
pp1$Region.A = factor(pp1$Region.A, levels = rownames(EC8_new[[1]]))

ggplot(pp1, aes(x = Region.A, y = Region.B, fill = CP)) +
  geom_tile() +
  labs(x = "Region A",
       y = "Region B",
       fill= "P(B|A)") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(CP,3)), color = "black", size = 4) +
  coord_fixed()

#Mouse 2
pp1 = data.frame(pp_m2$pprob_pair, row.names = rownames(EC8_new[[1]]))
names(pp1) = rownames(EC8_new[[1]])
print(pp1)

pp1 = data.frame(pp1, "Region B" = factor(rownames(EC8_new[[1]]), levels = rownames(EC8_new[[1]])))
names(pp1)[1:R] = rownames(EC8_new[[1]])
pp1 <-  pivot_longer(pp1,
                     cols = !Region.B,
                     names_to = "Region.A", 
                     values_to = "CP"
)
pp1$Region.A = factor(pp1$Region.A, levels = rownames(EC8_new[[1]]))

ggplot(pp1, aes(x = Region.A, y = Region.B, fill = CP)) +
  geom_tile() +
  labs(x = "Region A",
       y = "Region B",
       fill= "P(B|A)") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(CP,3)), color = "black", size = 4) +
  coord_fixed()

#Mouse 3
pp1 = data.frame(pp_m3$pprob_pair, row.names = rownames(EC8_new[[1]]))
names(pp1) = rownames(EC8_new[[1]])
print(pp1)

pp1 = data.frame(pp1, "Region B" = factor(rownames(EC8_new[[1]]), levels = rownames(EC8_new[[1]])))
names(pp1)[1:R] = rownames(EC8_new[[1]])
pp1 <-  pivot_longer(pp1,
                     cols = !Region.B,
                     names_to = "Region.A", 
                     values_to = "CP"
)
pp1$Region.A = factor(pp1$Region.A, levels = rownames(EC8_new[[1]]))

ggplot(pp1, aes(x = Region.A, y = Region.B, fill = CP)) +
  geom_tile() +
  labs(x = "Region A",
       y = "Region B",
       fill= "P(B|A)") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(CP,3)), color = "black", size = 4) +
  coord_fixed()

#Mouse 4
pp1 = data.frame(pp_m4$pprob_pair, row.names = rownames(EC8_new[[1]]))
names(pp1) = rownames(EC8_new[[1]])
print(pp1)

pp1 = data.frame(pp1, "Region B" = factor(rownames(EC8_new[[1]]), levels = rownames(EC8_new[[1]])))
names(pp1)[1:R] = rownames(EC8_new[[1]])
pp1 <-  pivot_longer(pp1,
                     cols = !Region.B,
                     names_to = "Region.A", 
                     values_to = "CP"
)
pp1$Region.A = factor(pp1$Region.A, levels = rownames(EC8_new[[1]]))

ggplot(pp1, aes(x = Region.A, y = Region.B, fill = CP)) +
  geom_tile() +
  labs(x = "Region A",
       y = "Region B",
       fill= "P(B|A)") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(CP,3)), color = "black", size = 4) +
  coord_fixed()

#Mouse 5
pp1 = data.frame(pp_m5$pprob_pair, row.names = rownames(EC8_new[[1]]))
names(pp1) = rownames(EC8_new[[1]])
print(pp1)

pp1 = data.frame(pp1, "Region B" = factor(rownames(EC8_new[[1]]), levels = rownames(EC8_new[[1]])))
names(pp1)[1:R] = rownames(EC8_new[[1]])
pp1 <-  pivot_longer(pp1,
                     cols = !Region.B,
                     names_to = "Region.A", 
                     values_to = "CP"
)
pp1$Region.A = factor(pp1$Region.A, levels = rownames(EC8_new[[1]]))

ggplot(pp1, aes(x = Region.A, y = Region.B, fill = CP)) +
  geom_tile() +
  labs(x = "Region A",
       y = "Region B",
       fill= "P(B|A)") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(CP,3)), color = "black", size = 4) +
  coord_fixed()

#Mouse 6
pp1 = data.frame(pp_m6$pprob_pair, row.names = rownames(EC8_new[[1]]))
names(pp1) = rownames(EC8_new[[1]])
print(pp1)

pp1 = data.frame(pp1, "Region B" = factor(rownames(EC8_new[[1]]), levels = rownames(EC8_new[[1]])))
names(pp1)[1:R] = rownames(EC8_new[[1]])
pp1 <-  pivot_longer(pp1,
                     cols = !Region.B,
                     names_to = "Region.A", 
                     values_to = "CP"
)
pp1$Region.A = factor(pp1$Region.A, levels = rownames(EC8_new[[1]]))

ggplot(pp1, aes(x = Region.A, y = Region.B, fill = CP)) +
  geom_tile() +
  labs(x = "Region A",
       y = "Region B",
       fill= "P(B|A)") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(CP,3)), color = "black", size = 4) +
  coord_fixed()

group.colors <- c('1'= "#FF3000",'2'= "#FF9900",'3'= "#FFDB6D",'4'= "#009999",'5'= "#333BFF",'6'= "#9633FF")

df = data.frame(p = c(pp_m1$pprob,pp_m2$pprob,pp_m3$pprob,pp_m4$pprob,pp_m5$pprob,pp_m6$pprob))
df$Region = factor(rep(rownames(EC8_new[[1]]), M), levels = rownames(EC8_new[[1]]))
df$Mouse = c(rep('1', R),rep('2', R),rep('3', R),rep('4', R),rep('5', R),rep('6', R))
df$lower = c(pp_m1$pprob_ci[1,],pp_m2$pprob_ci[1,],pp_m3$pprob_ci[1,],pp_m4$pprob_ci[1,],pp_m5$pprob_ci[1,],pp_m6$pprob_ci[1,])
df$upper = c(pp_m1$pprob_ci[2,],pp_m2$pprob_ci[2,],pp_m3$pprob_ci[2,],pp_m4$pprob_ci[2,],pp_m5$pprob_ci[2,],pp_m6$pprob_ci[2,])

ggplot(df, mapping = aes(x = Region, y= p, color = Mouse, group = Mouse)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_color_manual(values=group.colors) +
  ylab('P(y>0)') +
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




estimated.pp.plot <- ggplot2::ggplot(estimated.projection.df,
                                     mapping = aes(x = factor(region, levels = regions.name),
                                                   y = projection.med,
                                                   group = cluster,
                                                   color = class)) +
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin = projection.lower,
                    ymax = projection.upper),
                width = 0.1)+
  theme_bw()+
  ylim(c(0,1))+
  ylab('projection strength')+
  xlab('region')+
  facet_wrap(~cluster)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=axis.text.size),
        axis.title=element_text(size=axis.title.size),
        plot.title = element_text(size=plot.title.size))
