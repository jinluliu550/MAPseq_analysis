##############################################
## Brain Connectivity - MAPseq
## Generate data with different projection patterns and mouse differences
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 

##############################################
## Simulate data 
##############################################

M = 4 # number of mice
n = c(200,120, 100, 150) # number of neurons
R = 3 # number of regions
J = 2 + 4 + 1 #number of motifs (unicasting, bicasting, multicasting)

# Projection motifs
q = matrix(0,J, R)
q[1,] = c(1,0,0)
q[2,] = c(0,0,1)
# bicasting to regions 1 and 2 but with different strengths (combined in binomial)
q[3,] = c(.8,.2,0)
q[4,] = c(.3,.7,0)
# bicasting to regions 3 and 2 but with different strengths (combined in binomial)
q[5,] = c(0,.2,.8)
q[6,] = c(0,.7,0.3)
# multicast
q[7,] = c(0.3,0.4,0.3)
df = data.frame(motif = rep(as.factor(c(1:J)),R), region = rep(as.factor(1:R),each=J), prob = as.vector(q))
ggplot(df, aes(x=region, y=motif, fill=prob)) +
  geom_tile() +
  theme_void()+
  scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                       values = scales::rescale(c(0, 0.5, 1)),
                       limits = c(0,1)) +
  guides(fill=guide_legend(title="Projection\nprobability"))+
  xlab('region')+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.title.y = element_text(size = 12,angle = 90))+
  labs(y = "motif")

# Dispersion
gamma_star = matrix(50,J, 1)

# Mouse proportions (no difference between 2 and 3, but 1 and 4 are different)
w = matrix(0, J, M)
w[,1] = c(.3,0,.3,.2,0,.1,0.1)
w[,2] = c(0.05,.05,.1,.2,.1,.2,.3)
w[,3] = c(0.05,.05,.1,.2,.1,.2,.3)
w[,4] = c(0,.3,0,.1,.3,.2,0.1)
df = data.frame(motif = rep(as.factor(c(1:J)),M), mouse = rep(as.factor(1:M),each=J), w = as.vector(w))
ggplot(df, aes(x=mouse, y=motif, fill=w)) +
  geom_tile() +
  theme_void()+
  scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                       values = scales::rescale(c(0, 0.15, 0.3)),
                       limits = c(0,0.4)) +
  guides(fill=guide_legend(title="Motif\nproportion"))+
  xlab('mouse')+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.title.y = element_text(size = 12,angle = 90))+
  labs(y = "motif")


##############################################
## Simulate data from Dir-multinomial mixture
##############################################

library(rBeta2009)

rdirichlet_sparse = function(n,a){
  K = length(a)
  ind = a>0
  p = matrix(0, n, K)
  if(sum(ind)>1){
    p[,ind] = rdirichlet(n,a[ind])
    if(sum(ind)==2){
      p[,which(ind)[2]] = 1- p[,which(ind)[2]]}
  }else{
    p[,ind]=1
  }
  return(p)
}

# Simulate allocations and counts
Z = list(M)
Y = list(M)
N = list(M)
for (m in c(1:M)){
  Z[[m]] =  apply(rmultinom(n[m],1, w[,m]),2,which.max)
  Y[[m]] = matrix(0,n[m],R)
  N[[m]] = floor(runif(n[m],50,201))
  for (j in c(1:J)){
    njm = sum(Z[[m]]==j)
    if(njm>0){
      qtilde = rdirichlet_sparse(njm,gamma_star[j]*q[j,])
      Y[[m]][Z[[m]]==j,] = t(apply(cbind(qtilde,N[[m]][Z[[m]]==j]), 1, function(x){rmultinom(1,x[length(x)],x[-length(x)])}))
    }
  }
}

########################
## Plot data

#sort neurons

sort_neurons = function(p, R){
  max_strength = apply(p,2,max)
  max_region = apply(p,2,which.max)
  s_ind = sort(max_region, index.return=T)$ix
  for (r in c(1:R)){
    if (sum(max_region==r)>1){
      sr_ind = sort(max_strength[max_region==r],decreasing = TRUE,index.return=T)$ix
      s_ind[max_region[s_ind]==r] = s_ind[max_region[s_ind]==r][sr_ind]
    }
  }
  return(p[,s_ind])
}

phat = t(Y[[1]])/t(matrix(colSums(t(Y[[1]])),n[1], R))
row.names(phat) = as.factor(c(1:R))
sphat = sort_neurons(phat,R)
df1 <- data.frame(region = rep(row.names(phat), n[[1]]), cell = rep(1:n[[1]], each = R), ps = as.vector(sphat))
gp1 = ggplot(df1, mapping = aes(x = factor(region, levels = row.names(phat)), y = cell, fill = ps))+
  geom_tile()+
  theme_void()+
  scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                       values = scales::rescale(c(0, 0.5, 1)),
                       limits = c(0,1))+
  guides(fill=guide_legend(title="Projection\nstrength"))+
  xlab('region')+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.title.y = element_text(size = 12,angle = 90))+
  labs(title = paste0("Mouse ",1), y = "neurons")

phat2 = t(Y[[2]])/t(matrix(colSums(t(Y[[2]])),n[2], R))
row.names(phat2) = as.factor(c(1:R))
sphat2 = sort_neurons(phat2,R)
df2 <- data.frame(region = rep(row.names(phat), n[[2]]), cell = rep(1:n[[2]], each = R), ps = as.vector(sphat2))
gp2 = ggplot(df2, mapping = aes(x = factor(region, levels = row.names(phat)), y = cell, fill = ps))+
  geom_tile()+
  theme_void()+
  scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                       values = scales::rescale(c(0, 0.5, 1)),
                       limits = c(0,1))+
  guides(fill=guide_legend(title="Projection\nstrength"))+
  xlab('region')+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.title.y = element_text(size = 12,angle = 90))+
  labs(title = paste0("Mouse ",2), y = "neurons")

phat3 = t(Y[[3]])/t(matrix(colSums(t(Y[[3]])),n[3], R))
row.names(phat3) = as.factor(c(1:R))
sphat3 = sort_neurons(phat3,R)
df3 <- data.frame(region = rep(row.names(phat3), n[[3]]), cell = rep(1:n[[3]], each = R), ps = as.vector(sphat3))
gp3 = ggplot(df3, mapping = aes(x = factor(region, levels = row.names(phat)), y = cell, fill = ps))+
  geom_tile()+
  theme_void()+
  scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                       values = scales::rescale(c(0, 0.5, 1)),
                       limits = c(0,1))+
  guides(fill=guide_legend(title="Projection\nstrength"))+
  xlab('region')+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.title.y = element_text(size = 12,angle = 90))+
  labs(title = paste0("Mouse ",3), y = "neurons")

phat4 = t(Y[[4]])/t(matrix(colSums(t(Y[[4]])),n[4], R))
row.names(phat4) = as.factor(c(1:R))
sphat4 = sort_neurons(phat4,R)
df4 <- data.frame(region = rep(row.names(phat4), n[[4]]), cell = rep(1:n[[4]], each = R), ps = as.vector(sphat4))
gp4 = ggplot(df4, mapping = aes(x = factor(region, levels = row.names(phat)), y = cell, fill = ps))+
  geom_tile()+
  theme_void()+
  scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                       values = scales::rescale(c(0, 0.5, 1)),
                       limits = c(0,1))+
  guides(fill=guide_legend(title="Projection\nstrength"))+
  xlab('region')+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.title.y = element_text(size = 12,angle = 90))+
  labs(title = paste0("Mouse ",4), y = "neurons")

library(ggpubr)
ggarrange(gp1, gp2, gp3, gp4, ncol=4, nrow=1, common.legend = TRUE, legend="right")

###########################################
# Compute true total variation distance
mytv_dist = function(x,ind){
  xdim = dim(x)
  y = matrix(x[,ind],xdim[1],xdim[2])
  return(0.5*colSums(abs(x-y)))
}

tv_true = lapply(c(1:M),mytv_dist,x=w)
tv_true = data.frame(matrix(unlist(tv_true), nrow=M, byrow=TRUE))

tv_true = data.frame(tv_true, "Mouse 1" = c("1", "2", "3", "4"))
names(tv_true)[1:M] = c("1", "2", "3", "4")
tv_true <-  pivot_longer(tv_true,
                         cols = !Mouse.1,
                         names_to = "Mouse.2", 
                         values_to = "TV"
)
ggplot(tv_true, aes(x = Mouse.1, y = Mouse.2, fill = TV)) +
  geom_tile() +
  labs(x = "Mouse",
       y = "Mouse") +
  theme_bw() +
  scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                       values = scales::rescale(c(0, 0.5, 1)),
                       limits = c(0,1)) +
  geom_text(aes(label = round(tv_true$TV,3)), color = "black", size = 4) +
  coord_fixed()
