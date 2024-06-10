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


# CI for the difference in the proability of non-zero counts for region one given another
# between two mice
projection_prob_d= function(m1,m2,N,mcmc_all_out){
  
  R= mcmc_all_out$R
  I = length(mcmc_all_out$alpha_output)
  
  #output
  pprob_pair_array_1 = array(0,dim= c(I,R,R))
  pprob_mat_1 = matrix(0,I,R)
  
  pprob_pair_array_2 = array(0,dim= c(I,R,R))
  pprob_mat_2 = matrix(0,I,R)
  
  for(i in 1:I){
    wjm1 = mcmc_all_out$omega_J_M_output[[i]][,m1]
    wjm2 = mcmc_all_out$omega_J_M_output[[i]][,m2]
    qstarjm = mcmc_all_out$q_star_1_J_output[[i]]
    gammastarjm = mcmc_all_out$gamma_star_1_J_output[[i]]
    pzero_i1 = ppzero(N,qstarjm*gammastarjm, wjm1)
    pzero_i2 = ppzero(N,qstarjm*gammastarjm, wjm2)
    pzero_pair_i1 = ppzero_pair(N,qstarjm*gammastarjm, wjm1)
    pzero_pair_i2 = ppzero_pair(N,qstarjm*gammastarjm, wjm2)
    pprob_mat_1[i,] = 1- pzero_i1
    pprob_mat_2[i,] = 1- pzero_i2
    pprob_pair_temp1 =  (1-(matrix(pzero_i1, R,R) + t(matrix(pzero_i1, R,R)) - pzero_pair_i1))/matrix(1-pzero_i1, R,R)
    pprob_pair_temp2 =  (1-(matrix(pzero_i2, R,R) + t(matrix(pzero_i2, R,R)) - pzero_pair_i2))/matrix(1-pzero_i2, R,R)
    diag(pprob_pair_temp1) = 1
    diag(pprob_pair_temp2) = 1
    pprob_pair_array_1[i,,] = pprob_pair_temp1
    pprob_pair_array_2[i,,] = pprob_pair_temp2 
  }
  pprob_d = apply(pprob_mat_1-pprob_mat_2,2,mean)
  pprob_pair_d = apply(pprob_pair_array_1-pprob_pair_array_2,c(2,3),mean)
  pprob_d_ci = apply(pprob_mat_1-pprob_mat_2,2,function(x){HPDinterval(as.mcmc(x))})
  pprobpair_d_ci = apply(pprob_pair_array_1-pprob_pair_array_2,c(2,3),function(x){HPDinterval(as.mcmc(x))})
  return(list(pprob = pprob_d, pprob_pair = pprob_pair_d, pprob_ci = pprob_d_ci, pprob_pair_ci = pprobpair_d_ci))
}

# Compute expected number of counts/N for each region for a DirMult mixture
epstrength <- function(alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  eps = matrix(0,R,1)
  for (j in 1:J){
    eps = eps+w[j]*alpha[j,]/sum(alpha[j,])
  }
  return(eps)
}

# Sample from DirMult mixture
rpstrength <- function(n, alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  j = sample(c(1:J),1, prob = w)
  ps = rdirmnom(1,n,alpha[j,])/n
  return(ps)
}


# Compute the variance of counts/N for each region for a DirMult mixture
vpstrength <- function(n,alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  eps = matrix(0,R,1)
  vps = matrix(0,R,1)
  for (j in 1:J){
    gamma0 = sum(alpha[j,])
    eps = eps+w[j]*alpha[j,]/gamma0
    vps = vps+w[j]*(1/n*(alpha[j,]/gamma0)*(1-alpha[j,]/gamma0)*(n+gamma0)/(1+gamma0)+(alpha[j,]/gamma0)^2)
  }
  vps = vps - eps^2
  return(vps)
}


# Compute the covariance of y1/N and y2?N counts for each pair of region for a DirMult mixture
cpstrength <- function(n, alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  cps = matrix(0,R,R)
  eps = matrix(0,R,1)
  for (j in 1:J){
    gamma0 = sum(alpha[j,])
    eps = eps+w[j]*alpha[j,]/gamma0
    for (r1 in 1:R){
      cps[r1,] =  cps[r1,] + w[j]*(-1/n*(alpha[j,r1]/gamma0)*(alpha[j,]/gamma0)*(n+gamma0)/(1+gamma0)+(alpha[j,r1]/gamma0)*(alpha[j,]/gamma0))
    }
  }
  cps = cps - eps%*%t(eps)
  diag(cps) = vpstrength(n,alpha, w)
  return(cps)
}

# Compute the posterior expected number of counts/N for each region 
post_epstrength = function(m,N,mcmc_all_out){
  
  R= mcmc_all_out$R
  I = length(mcmc_all_out$alpha_output)
  
  #output
  eps_mat = matrix(0,I,R)
  rps_mat = matrix(0,I,R)
  
  for(i in 1:I){
    wjm = mcmc_all_out$omega_J_M_output[[i]][,m]
    qstarjm = mcmc_all_out$q_star_1_J_output[[i]]
    gammastarjm = mcmc_all_out$gamma_star_1_J_output[[i]]
    eps_mat[i,] = epstrength(qstarjm*gammastarjm, wjm)
    rps_mat[i,] = rpstrength(N,qstarjm*gammastarjm, wjm)
  }
  eps = apply(eps_mat,2,mean)
  eps_ci = apply(eps_mat,2,function(x){HPDinterval(as.mcmc(x))})
  rps_ci = apply(rps_mat,2,function(x){HPDinterval(as.mcmc(x))})
  return(list(eps = eps, eps_ci = eps_ci, rps_ci = rps_ci))
}

# Compute the posterior covariance of counts/N for each pair of regions
post_cpstrength= function(m,N,mcmc_all_out){
  
  R= mcmc_all_out$R
  I = length(mcmc_all_out$alpha_output)
  
  #output
  cps_array = array(0,dim= c(I,R,R))
  eps_mat = matrix(0,I,R)
  eps2_array = array(0,dim= c(I,R,R))
  
  for(i in 1:I){
    wjm = mcmc_all_out$omega_J_M_output[[i]][,m]
    qstarjm = mcmc_all_out$q_star_1_J_output[[i]]
    gammastarjm = mcmc_all_out$gamma_star_1_J_output[[i]]
    eps_mat[i,] = epstrength(qstarjm*gammastarjm, wjm)
    eps2_array[i,,] = eps_mat[i,]%*%t(eps_mat[i,]) 
    cps_array[i,,] = cpstrength(N,qstarjm*gammastarjm, wjm) 
  }
  eps = apply(eps_mat,2,mean)
  eps2 = apply(eps2_array,c(2,3),mean)
  cps = apply(cps_array,c(2,3),mean) + eps2 - eps%*%t(eps)
  cps_ci = apply(cps_array,c(2,3),function(x){HPDinterval(as.mcmc(x))})
  return(list(cps = cps, cps_ci = cps_ci))
}


#Mouse 1
pp_m1 = projection_prob(1,100,mcmc_all_hans)
ps_m1 = post_epstrength(1,100,mcmc_all_hans)
cps_m1 = post_cpstrength(1,100,mcmc_all_hans)

#Mouse 2
pp_m2 = projection_prob(2,100,mcmc_all_hans)
ps_m2 = post_epstrength(2,100,mcmc_all_hans)
cps_m2 = post_cpstrength(2,100,mcmc_all_hans)

#Mouse 3
pp_m3 = projection_prob(3,100,mcmc_all_hans)
ps_m3 = post_epstrength(3,100,mcmc_all_hans)
cps_m3 = post_cpstrength(3,100,mcmc_all_hans)


## Plots

R <- nrow(data_Hans[[1]])

#Mouse 1
pp1 = data.frame(pp_m1$pprob_pair, row.names = rownames(data_Hans[[1]]))
names(pp1) = rownames(data_Hans[[1]])
print(pp1)

pp1 = data.frame(pp1, "Region B" = factor(rownames(data_Hans[[1]]), levels = rownames(data_Hans[[1]])))
names(pp1)[1:R] = rownames(data_Hans[[1]])
pp1 <-  pivot_longer(pp1,
                     cols = !Region.B,
                     names_to = "Region.A", 
                     values_to = "CP"
)
pp1$Region.A = factor(pp1$Region.A, levels = rownames(data_Hans[[1]]))


ggplot(pp1, aes(x = Region.A, y = Region.B, fill = CP)) +
  geom_tile() +
  labs(x = "Region A",
       y = "Region B",
       fill= "P(A|B)") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(CP,3)), color = "black", size = 4) +
  coord_fixed()



#Mouse 2
pp1 = data.frame(pp_m2$pprob_pair, row.names = rownames(data_Hans[[1]]))
names(pp1) = rownames(data_Hans[[1]])
print(pp1)

pp1 = data.frame(pp1, "Region B" = factor(rownames(data_Hans[[1]]), levels = rownames(data_Hans[[1]])))
names(pp1)[1:R] = rownames(data_Hans[[1]])
pp1 <-  pivot_longer(pp1,
                     cols = !Region.B,
                     names_to = "Region.A", 
                     values_to = "CP"
)
pp1$Region.A = factor(pp1$Region.A, levels = rownames(data_Hans[[1]]))

ggplot(pp1, aes(x = Region.A, y = Region.B, fill = CP)) +
  geom_tile() +
  labs(x = "Region A",
       y = "Region B",
       fill= "P(A|B)") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(CP,3)), color = "black", size = 4) +
  coord_fixed()

#Mouse 3
pp1 = data.frame(pp_m3$pprob_pair, row.names = rownames(data_Hans[[1]]))
names(pp1) = rownames(data_Hans[[1]])
print(pp1)

pp1 = data.frame(pp1, "Region B" = factor(rownames(data_Hans[[1]]), levels = rownames(data_Hans[[1]])))
names(pp1)[1:R] = rownames(data_Hans[[1]])
pp1 <-  pivot_longer(pp1,
                     cols = !Region.B,
                     names_to = "Region.A", 
                     values_to = "CP"
)
pp1$Region.A = factor(pp1$Region.A, levels = rownames(data_Hans[[1]]))

ggplot(pp1, aes(x = Region.A, y = Region.B, fill = CP)) +
  geom_tile() +
  labs(x = "Region A",
       y = "Region B",
       fill= "P(A|B)") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(CP,3)), color = "black", size = 4) +
  coord_fixed()


# Projection probability across mice

group.colors <- c('1'= "#FF3000",'2'= "#FF9900",'3'= "#FFDB6D")

df = data.frame(p = c(pp_m1$pprob,pp_m2$pprob,pp_m3$pprob))
df$Region = factor(rep(rownames(data_Hans[[1]]), M), levels = rownames(data_Hans[[1]]))
df$Mouse = c(rep('1', R),rep('2', R),rep('3', R))
df$lower = c(pp_m1$pprob_ci[1,],pp_m2$pprob_ci[1,],pp_m3$pprob_ci[1,])
df$upper = c(pp_m1$pprob_ci[2,],pp_m2$pprob_ci[2,],pp_m3$pprob_ci[2,])

ggplot(df, mapping = aes(x = Region, y= p, color = Mouse, group = Mouse)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_color_manual(values=group.colors) +
  ylab('P(y>0)') +
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Projection strength across mice

group.colors <- c('1'= "#FF3000",'2'= "#FF9900",'3'= "#FFDB6D")

df = data.frame(p = c(ps_m1$eps,ps_m2$eps,ps_m3$eps))
df$Region = factor(rep(rownames(data_Hans[[1]]), M), levels = rownames(data_Hans[[1]]))
df$Mouse = c(rep('1', R),rep('2', R),rep('3', R))
df$lower = c(ps_m1$eps_ci[1,],ps_m2$eps_ci[1,],ps_m3$eps_ci[1,])
df$upper = c(ps_m1$eps_ci[2,],ps_m2$eps_ci[2,],ps_m3$eps_ci[2,])

ggplot(df, mapping = aes(x = Region, y= p, color = Mouse, group = Mouse)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_color_manual(values=group.colors) +
  ylab('E[y/N]') +
  geom_errorbar(aes(ymin = lower, ymax = upper),width = 0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


