##############################################
## Hans paper - Exploratory data analysis
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 

#   
# Hans_data1 <- read.csv('./data/Han-data/han_brain4_gh.csv')
# Hans_data2 <- read.csv('./data/Han-data/han_brain5_gh.csv')
# Hans_data3 <- read.csv('./data/Han-data/han_brain6_gh.csv')
# 
# data_Hans <- list(t(Hans_data1),
#                   t(Hans_data2),
#                   t(Hans_data3))

M = 4 # number of Mice
R = dim(data_Hans[[1]])[1] # number of regions
# number of neurons
n <- list(t(dim(data_Hans[[1]])[2]),
          t(dim(data_Hans[[2]])[2]),
          t(dim(data_Hans[[3]])[2]),
          t(dim(data_Hans[[4]])[2]))

# Histogram of empirical projection strengths
#Mouse 1
phat = data_Hans[[1]]/t(matrix(colSums(data_Hans[[1]]),n[[1]], R))
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

# Mouse 2
phat2 = data_Hans[[2]]/t(matrix(colSums(data_Hans[[2]]),n[[2]], R))
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

# Mouse 3
phat3 = data_Hans[[3]]/t(matrix(colSums(data_Hans[[3]]),n[[3]], R))
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


# Mouse 4
phat4 = data_Hans[[4]]/t(matrix(colSums(data_Hans[[4]]),n[[4]], R))
p41 <- ggplot() +
  geom_histogram(aes(x = phat4[1,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",4,", Region = ", row.names(phat)[1])) 
p42 <- ggplot() +
  geom_histogram(aes(x = phat4[2,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",4,", Region = ", row.names(phat)[2])) 
p43 = ggplot() +
  geom_histogram(aes(x = phat4[3,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",4,", Region = ", row.names(phat)[3])) 
p44 = ggplot() +
  geom_histogram(aes(x = phat4[4,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",4,", Region = ", row.names(phat)[4])) 
p45 = ggplot() +
  geom_histogram(aes(x = phat4[5,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",4,", Region = ", row.names(phat)[5])) 
p46 = ggplot() +
  geom_histogram(aes(x = phat4[6,])) +
  theme_bw() +
  labs(x ="y/N", title = paste0("Mouse = ",4,", Region = ", row.names(phat)[6])) 

library(patchwork)

png(file = './plots/Hans/eda/hist_li_lm.png',
    width = 1000,
    height = 400)

p11+ p21 + p31 + p41 + p12 + p22 + p32 + p42 + plot_layout(ncol = 4)

dev.off()


png(file = './plots/Hans/eda/hist_al_pm.png',
    width = 1000,
    height = 400)

p13+ p23 + p33 + p43 + p14 + p24 + p34  + p44 + plot_layout(ncol = 4)

dev.off()

png(file = './plots/Hans/eda/hist_am_rl.png',
    width = 1000,
    height = 400)

p15+ p25 + p35 + p45 + p16 + p26 + p36 + p46 + plot_layout(ncol = 4)

dev.off()

# Scatter plots in 2d - space

df = data.frame(x = c(phat[1,],phat2[1,],phat3[1,],phat4[1,]), 
                y = c(phat[2,],phat2[2,],phat3[2,],phat4[2,]),
                mouse = as.factor(c(rep(1,n[[1]]), rep(2,n[[2]]), rep(3,n[[3]]), rep(4,n[[4]]))))

png(file = './plots/Hans/eda/scatter_li_lm.png',
    width = 450,
    height = 400)

ggplot(df) +
  geom_point(aes(x = x, y = y, col = mouse)) +
  theme_bw() +
  labs(x = row.names(phat)[1], y = row.names(phat)[2]) 

dev.off()



df = data.frame(x = c(phat[3,],phat2[3,],phat3[3,],phat4[3,]), 
                y = c(phat[4,],phat2[4,],phat3[4,],phat4[4,]), 
                mouse = as.factor(c(rep(1,n[[1]]), rep(2,n[[2]]), rep(3,n[[3]]), rep(4,n[[4]]))))

png(file = './plots/Hans/eda/scatter_al_pm.png',
    width = 450,
    height = 400)

ggplot(df) +
  geom_point(aes(x = x, y = y, col = mouse)) +
  theme_bw() +
  labs(x = row.names(phat)[3], y = row.names(phat)[4]) 

dev.off()



df = data.frame(x = c(phat[5,],phat2[5,],phat3[5,],phat4[5,]), 
                y = c(phat[6,],phat2[6,],phat3[6,],phat4[6,]), 
                mouse = as.factor(c(rep(1,n[[1]]), rep(2,n[[2]]), rep(3,n[[3]]), rep(4,n[[4]]))))


png(file = './plots/Hans/eda/scatter_am_rl.png',
    width = 450,
    height = 400)

ggplot(df) +
  geom_point(aes(x = x, y = y, col = mouse)) +
  theme_bw() +
  labs(x = row.names(phat)[5], y = row.names(phat)[6]) 

dev.off()

# Heat map

#sort neurons

# gel_plot_hans <- gel_plot(data_Hans)
# 
# 
# png(file = './plots/Hans/eda/gelplot_hans.png',
#     width = 1000,
#     height = 1000/3)
# 
# ggarrange(gel_plot_hans[[1]],
#           gel_plot_hans[[2]],
#           gel_plot_hans[[3]],
#           nrow = 1,
#           widths = c(1,1,1.3))
# 
# dev.off()

# sort_neurons = function(p, R){
#   max_strength = apply(p,2,max)
#   max_region = apply(p,2,which.max)
#   s_ind = sort(max_region, index.return=T)$ix
#   for (r in c(1:R)){
#     if (sum(max_region==r)>1){
#       sr_ind = sort(max_strength[max_region==r],decreasing = TRUE,index.return=T)$ix
#       s_ind[max_region[s_ind]==r] = s_ind[max_region[s_ind]==r][sr_ind]
#     }
#   }
#   return(p[,s_ind])
# }
# 
# sphat = sort_neurons(phat,R)
# df1 <- data.frame(region = rep(row.names(phat), n[[1]]), cell = rep(1:n[[1]], each = R), ps = as.vector(sphat))
# gp1 = ggplot(df1, mapping = aes(x = factor(region, levels = row.names(phat)), y = cell, fill = ps))+
#   geom_tile()+
#   theme_void()+
#   scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
#                        values = scales::rescale(c(0, 0.5, 1)),
#                        limits = c(0,1))+
#   guides(fill=guide_legend(title="Projection\nstrength"))+
#   xlab('region')+
#   theme(axis.text.x = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         plot.title = element_text(size=12, hjust = 0.5),
#         axis.title.y = element_text(size = 12,angle = 90))+
#   labs(title = paste0("Mouse ",1), y = "neurons")
# 
# sphat2 = sort_neurons(phat2,R)
# df2 <- data.frame(region = rep(row.names(phat), n[[2]]), cell = rep(1:n[[2]], each = R), ps = as.vector(sphat2))
# gp2 = ggplot(df2, mapping = aes(x = factor(region, levels = row.names(phat)), y = cell, fill = ps))+
#   geom_tile()+
#   theme_void()+
#   scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
#                        values = scales::rescale(c(0, 0.5, 1)),
#                        limits = c(0,1))+
#   guides(fill=guide_legend(title="Projection\nstrength"))+
#   xlab('region')+
#   theme(axis.text.x = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         plot.title = element_text(size=12, hjust = 0.5),
#         axis.title.y = element_text(size = 12,angle = 90))+
#   labs(title = paste0("Mouse ",2), y = "neurons")
# 
# sphat3 = sort_neurons(phat3,R)
# df3 <- data.frame(region = rep(row.names(phat3), n[[3]]), cell = rep(1:n[[3]], each = R), ps = as.vector(sphat3))
# gp3 = ggplot(df3, mapping = aes(x = factor(region, levels = row.names(phat)), y = cell, fill = ps))+
#   geom_tile()+
#   theme_void()+
#   scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
#                        values = scales::rescale(c(0, 0.5, 1)),
#                        limits = c(0,1))+
#   guides(fill=guide_legend(title="Projection\nstrength"))+
#   xlab('region')+
#   theme(axis.text.x = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         plot.title = element_text(size=12, hjust = 0.5),
#         axis.title.y = element_text(size = 12,angle = 90))+
#   labs(title = paste0("Mouse ",3), y = "neurons")
# 
# library(ggpubr)
# ggarrange(gp1, gp2, gp3, ncol=3, nrow=1, common.legend = TRUE, legend="right")
  

