# Suppose there are four mice
M <- 4

# And 100 neurons in mouse 1 and 4, 200 neurons in mouse 2 and 3
C <- c(100,200,200,100)

# and 3 brain regions
R <- 3

# There are three clusters
J <- 3

# with projection strengths
q <- matrix(c(0,0.4,0.6,
              0.2,0.5,0.3,
              0.2,0,0.8),
            
            byrow = TRUE,
            nrow = 3,
            ncol = 3)

# Mouse 1: LEC only
# Mouse 2: neuron 1:100 LEC, neuron 101:200 MEC
# Mouse 3: neuron 1:100 LEC, neuron 101:200 MEC
# Mouse 4: MEC only

delta <- matrix(c(0.2,0.5,0.3,
                  0.5,0,0.5,
                  0.5,0.2,0.3,
                  0.2,0.3,0.5),
                
                nrow = J,
                ncol = M)


# Gamma
gamma <- rep(10,J)

# Effect of changing from LEC to MEC
beta <- c(-2, 2, 0)

# w_jm(xi = 0)
w_jm_LEC <- matrix(c(0.2,0.5,0.3,
                     0.5,0,0.5,
                     0.5,0.2,0.3,
                     NA,NA,NA),
                   
                   nrow = J,
                   ncol = M)

# w_jm(xi = 1)
w_j1_MEC <- rep(NA,3)

w_j2_MEC <- (delta[,2]*exp(beta))/sum(delta[,2]*exp(beta))

w_j3_MEC <- (delta[,3]*exp(beta))/sum(delta[,3]*exp(beta))

w_j4_MEC <- (delta[,4]*exp(beta))/sum(delta[,4]*exp(beta))

w_jm_MEC <- matrix(c(w_j1_MEC,
                     w_j2_MEC,
                     w_j3_MEC,
                     w_j4_MEC),
                   
                   nrow = J,
                   ncol = M)

# Allocations
Z_1 <- sample(1:J,
              size = 100,
              replace = TRUE,
              prob = w_jm_LEC[,1])


