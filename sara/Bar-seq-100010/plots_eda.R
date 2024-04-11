##############################################
## Hans paper - Exploratory data analysis
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 

  
bar_data1 <- read.csv('./data/Bar-seq-100010/BRAIN_C9.csv')
bar_data2 <- read.csv('./data/Bar-seq-100010/BRAIN_C14.csv')
bar_data3 <- read.csv('./data/Bar-seq-100010/BRAIN_C28.csv')

data_bar <- list(t(bar_data1),
                  t(bar_data2),
                  t(bar_data3))

M = 3 # number of Mice
R = dim(data_bar[[1]])[1] # number of regions
# number of neurons
n <- list(t(dim(data_bar[[1]])[2]),
                  t(dim(data_bar[[2]])[2]),
                  t(dim(data_bar[[3]])[2]))

# Histogram of empirical projection strengths
#Mouse 1
phat = data_bar[[1]]/t(matrix(colSums(data_bar[[1]]),n[[1]], R))
p11 <- ggplot() +
  geom_histogram(aes(x = phat[1,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[1])) 
p12 <- ggplot() +
  geom_histogram(aes(x = phat[2,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[2])) 
p13 = ggplot() +
  geom_histogram(aes(x = phat[3,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[3])) 
p14 = ggplot() +
  geom_histogram(aes(x = phat[4,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[4])) 
p15 = ggplot() +
  geom_histogram(aes(x = phat[5,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[5])) 
p16 = ggplot() +
  geom_histogram(aes(x = phat[6,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[6])) 
p17 = ggplot() +
  geom_histogram(aes(x = phat[7,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[7])) 
p18 = ggplot() +
  geom_histogram(aes(x = phat[8,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[8])) 
p19 = ggplot() +
  geom_histogram(aes(x = phat[9,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[9])) 
p110 = ggplot() +
  geom_histogram(aes(x = phat[10,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[10])) 
p111 = ggplot() +
  geom_histogram(aes(x = phat[11,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",1,", Region = ", row.names(phat)[11])) 

# Mouse 2
phat2 = data_bar[[2]]/t(matrix(colSums(data_bar[[2]]),n[[2]], R))
p21 <- ggplot() +
  geom_histogram(aes(x = phat2[1,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[1])) 
p22 <- ggplot() +
  geom_histogram(aes(x = phat2[2,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[2])) 
p23 = ggplot() +
  geom_histogram(aes(x = phat2[3,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[3])) 
p24 = ggplot() +
  geom_histogram(aes(x = phat2[4,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[4])) 
p25 = ggplot() +
  geom_histogram(aes(x = phat2[5,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[5])) 
p26 = ggplot() +
  geom_histogram(aes(x = phat2[6,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[6])) 
p27 = ggplot() +
  geom_histogram(aes(x = phat2[7,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[7])) 
p28 = ggplot() +
  geom_histogram(aes(x = phat2[8,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[8])) 
p29 = ggplot() +
  geom_histogram(aes(x = phat2[9,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[9])) 
p210 = ggplot() +
  geom_histogram(aes(x = phat2[10,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[10])) 
p211 = ggplot() +
  geom_histogram(aes(x = phat2[11,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",2,", Region = ", row.names(phat)[11])) 

# Mouse 3
phat3 = data_bar[[3]]/t(matrix(colSums(data_bar[[3]]),n[[3]], R))
p31 <- ggplot() +
  geom_histogram(aes(x = phat3[1,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[1])) 
p32 <- ggplot() +
  geom_histogram(aes(x = phat3[2,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[2])) 
p33 = ggplot() +
  geom_histogram(aes(x = phat3[3,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[3])) 
p34 = ggplot() +
  geom_histogram(aes(x = phat3[4,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[4])) 
p35 = ggplot() +
  geom_histogram(aes(x = phat3[5,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[5])) 
p36 = ggplot() +
  geom_histogram(aes(x = phat3[6,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[6])) 
p37 = ggplot() +
  geom_histogram(aes(x = phat3[7,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[7])) 
p38 = ggplot() +
  geom_histogram(aes(x = phat3[8,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[8])) 
p39 = ggplot() +
  geom_histogram(aes(x = phat3[9,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[9])) 
p310 = ggplot() +
  geom_histogram(aes(x = phat3[10,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[10])) 
p311 = ggplot() +
  geom_histogram(aes(x = phat3[11,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",3,", Region = ", row.names(phat)[11])) 

library(patchwork)

p11+ p21 + p31 + p12 + p22 + p32 + plot_layout(ncol = 3)

p13+ p23 + p33 + p14 + p24 + p34 + plot_layout(ncol = 3)

p15+ p25 + p35 + p16 + p26 + p36 + plot_layout(ncol = 3)

p17+ p27 + p37 + p18 + p28 + p38 + plot_layout(ncol = 3)

p19+ p29 + p39 + p110 + p210 + p310 + plot_layout(ncol = 3)

p110+ p210 + p310 + p111 + p211 + p311 + plot_layout(ncol = 3)


# Scatter plots in 2d - space

df = data.frame(x = c(phat[1,],phat2[1,],phat3[1,]), y = c(phat[2,],phat2[2,],phat3[2,]), mouse = as.factor(c(rep(1,n[[1]]), rep(2,n[[2]]), rep(3,n[[3]]))))
ggplot(df) +
  geom_point(aes(x = x, y = y, col = mouse)) +
  theme_bw() +
  labs(x = row.names(phat)[1], y = row.names(phat)[2]) 

df = data.frame(x = c(phat[3,],phat2[3,],phat3[3,]), y = c(phat[4,],phat2[4,],phat3[4,]), mouse = as.factor(c(rep(1,n[[1]]), rep(2,n[[2]]), rep(3,n[[3]]))))
ggplot(df) +
  geom_point(aes(x = x, y = y, col = mouse)) +
  theme_bw() +
  labs(x = row.names(phat)[3], y = row.names(phat)[4]) 

df = data.frame(x = c(phat[5,],phat2[5,],phat3[5,]), y = c(phat[6,],phat2[6,],phat3[6,]), mouse = as.factor(c(rep(1,n[[1]]), rep(2,n[[2]]), rep(3,n[[3]]))))
ggplot(df) +
  geom_point(aes(x = x, y = y, col = mouse)) +
  theme_bw() +
  labs(x = row.names(phat)[5], y = row.names(phat)[6]) 

df = data.frame(x = c(phat[7,],phat2[7,],phat3[7,]), y = c(phat[8,],phat2[8,],phat3[8,]), mouse = as.factor(c(rep(1,n[[1]]), rep(2,n[[2]]), rep(3,n[[3]]))))
ggplot(df) +
  geom_point(aes(x = x, y = y, col = mouse)) +
  theme_bw() +
  labs(x = row.names(phat)[7], y = row.names(phat)[8]) 

df = data.frame(x = c(phat[9,],phat2[9,],phat3[9,]), y = c(phat[10,],phat2[10,],phat3[10,]), mouse = as.factor(c(rep(1,n[[1]]), rep(2,n[[2]]), rep(3,n[[3]]))))
ggplot(df) +
  geom_point(aes(x = x, y = y, col = mouse)) +
  theme_bw() +
  labs(x = row.names(phat)[9], y = row.names(phat)[10]) 

df = data.frame(x = c(phat[10,],phat2[10,],phat3[10,]), y = c(phat[11,],phat2[11,],phat3[11,]), mouse = as.factor(c(rep(1,n[[1]]), rep(2,n[[2]]), rep(3,n[[3]]))))
ggplot(df) +
  geom_point(aes(x = x, y = y, col = mouse)) +
  theme_bw() +
  labs(x = row.names(phat)[10], y = row.names(phat)[11]) 

# Heat map

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
  theme(axis.text.x = element_text(size = 12,angle = 90),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.title.y = element_text(size = 12,angle = 90))+
  labs(title = paste0("Mouse ",1), y = "neurons")

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
  theme(axis.text.x = element_text(size = 12,angle = 90),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.title.y = element_text(size = 12,angle = 90))+
  labs(title = paste0("Mouse ",2), y = "neurons")

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
  theme(axis.text.x = element_text(size = 12,angle = 90),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size=12, hjust = 0.5),
        axis.title.y = element_text(size = 12,angle = 90))+
  labs(title = paste0("Mouse ",3), y = "neurons")

library(ggpubr)
ggarrange(gp1, gp2, gp3, ncol=3, nrow=1, common.legend = TRUE, legend="right")
  

