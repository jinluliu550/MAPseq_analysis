# Demonstration of clusters split

figure1_bc_p1 <- matrix(0,
                        nrow = 4,
                        ncol = 300)

# For gene 1
figure1_bc_p1[2,1:150] <- runif(n = 150, min = 0.1, max = 0.3)
figure1_bc_p1[3,1:150] <- runif(n = 150, min = 0.3, max = 0.4)
figure1_bc_p1[4,1:150] <- 1-figure1_bc_p1[2,1:150]-figure1_bc_p1[3,1:150]

# For gene 2
figure1_bc_p1[2,151:300] <- runif(n = 150, min = 0.3, max = 0.5)
figure1_bc_p1[3,151:300] <- runif(n = 150, min = 0.1, max = 0.3)
figure1_bc_p1[4,151:300] <- 1-figure1_bc_p1[2,151:300]-figure1_bc_p1[3,151:300]


figure1_bc_p1_df <- data.frame(projection_strength = as.vector(figure1_bc_p1),
                               neuron = factor(rep(1:300, each = 4), levels = 1:300),
                               region = paste('region',1:4),
                               cluster = factor(c(rep(1, 600),
                                                  rep(2, 600)),
                                                levels = c(1:2)))

ggplot(figure1_bc_p1_df)+
  geom_line(mapping = aes(x = region,
                          y = projection_strength,
                          colour = cluster,
                          group = interaction(neuron, cluster)))+
  theme_bw()+
  xlab('region')+
  ylab('projection strengths')