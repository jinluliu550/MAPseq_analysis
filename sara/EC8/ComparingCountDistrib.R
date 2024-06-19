# packages
library(ggplot2)
library(survival)
library(pammtools)

## load data
EC8_new <- read.csv('./data/EC8_new/all_brains_setE_JL.csv')
region.name = names(EC8_new)[1:8]

# Convert to list
digit.list <- list(one = 1,
                   two = 2,
                   three = 3,
                   four = 4,
                   five = 5,
                   six = 6)

EC8_new_list <- lapply(1:6, 
                       function(m) EC8_new %>%
                         filter(brain == names(digit.list)[m]))

EC8_new <- lapply(1:6,
                  function(m) t(round(EC8_new_list[[m]][,1:8], 0)))

EC8_EC_label <- lapply(1:6,
                       function(m) EC8_new_list[[m]]$EC)


C <- sapply(1:6, function(m) ncol(EC8_new[[m]]))
R <- 8

#######################
## Compare Empirical CDF of counts

EC8_u = matrix(unlist(EC8_new), nrow= R)
EC8_EC_label_u = unlist(EC8_EC_label)


# Compute KS test
r =3 # region
ksl = ks.test(x= EC8_u[r,EC8_EC_label_u == "MEC"], y = EC8_u[r,EC8_EC_label_u == "LEC"], alternative = "less")
ksg = ks.test(x= EC8_u[r,EC8_EC_label_u == "MEC"], y = EC8_u[r,EC8_EC_label_u == "LEC"], alternative = "greater")

# Plot CDF
ggplot() +
  stat_ecdf(aes(x= EC8_u[r,] , colour = EC8_EC_label_u), pad = FALSE) +
  theme_bw() +
  xlim(0,100) +
  labs(title = paste0("KS- p-value: ",round(ksl$p.value,5),", KS+ p-value:",round(ksg$p.value,5)),
       x=region.name[r], y ="ECDF", color = "" )

# Plot CDF with confidence intervals
x <- EC8_u[r,EC8_EC_label_u == "MEC"]
ev <- rep(1, length(x))
sf_m <- survfit(Surv(x,ev)~1)
x <- EC8_u[r,EC8_EC_label_u == "LEC"]
ev <- rep(1, length(x))
sf_l <- survfit(Surv(x,ev)~1)

df = data.frame(x = c(sf_m$time,sf_l$time), y = c(1-sf_m$surv,1-sf_l$surv), 
                l = c(1-sf_m$upper,1-sf_l$upper), u = c(1-sf_m$lower,1-sf_l$lower),
                ind = c(rep("MEC", length(sf_m$time)),rep("LEC", length(sf_l$time))))
ggplot(df) +
  geom_step(aes(x =x, y = y, colour = ind)) +
  geom_stepribbon(aes(x = x, ymin = l, ymax = u, fill = ind), alpha = 0.5) +
  theme_bw() +
  xlim(0,100) +
  labs(title = paste0("KS- p-value: ",round(ksl$p.value,5),", KS+ p-value:",round(ksg$p.value,5)),
       x=region.name[r], y ="ECDF", color = "", fill = "" )


#######################
## Compare total neuron counts between MEC and LEC
EC8_totalcounts = lapply(EC8_new, colSums)
EC8_totalcounts_u = unlist(EC8_totalcounts)

# ALL
ggplot() +
  geom_bar(aes(x=EC8_totalcounts_u, fill = EC8_EC_label_u)) +
  theme_bw() +
  labs(x = "total neuron counts", fill = "")

# MOUSE 
m = 5
ggplot() +
  geom_bar(aes(x=EC8_totalcounts[[m]], fill = EC8_EC_label[[m]])) +
  theme_bw() +
  labs(x = "total neuron counts", fill = "")
