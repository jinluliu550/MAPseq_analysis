Rbp4_data <- read.csv('./data/Rbp4/Rbp4_data.csv')

Rbp4_data1 <- Rbp4_data %>%
  filter(atlas_no == 1) %>%
  select(X, Y)

# Scatter plot
ggplot(Rbp4_data1, aes(X, Y)) +
  geom_point() +
  scale_y_reverse() +
  scale_x_reverse()

# Breaks
min.x <- 250*floor(min(Rbp4_data1$X)/250)
max.x <- 250*ceiling(max(Rbp4_data1$X)/250)


df <- data.frame(x.group = cut(Rbp4_data1$X, breaks=seq(min.x,
                                                        max.x,
                                                        by = 250), include.lowest=TRUE),
                 y = Rbp4_data1$Y)

df$x.group <- factor(df$x.group,
                     levels = rev(unique(df$x.group)))

df2 <- df %>%
  group_by(x.group) %>%
  summarise(med = median(y))

# Draw segment
df_segment <- data.frame(x1 = df2$x.group[1:(nrow(df2)-1)],
                         x2 = df2$x.group[2:nrow(df2)],
                         y1 = df2$med[1:(nrow(df2)-1)],
                         y2 = df2$med[2:nrow(df2)])

df_segment <- data.frame(x.group = c(df2$x.group[1],
                                     rep(df2$x.group[2:(nrow(df2)-1)], 2),
                                     df2$x.group[nrow(df2)]),
                         
                         y = c(df2$med[1],
                               rep(df2$med[2:(nrow(df2)-1)], 2),
                               df2$med[nrow(df2)]))

mid_point <- rev(seq(min.x,
                     max.x,
                     by = 250))

mid_point <- sapply(1:6,
                    function(i) mean(c(mid_point[i], mid_point[i+1])))


df_segment <- df_segment %>%
  mutate(mid.point.x1 = mid_point[1:5],
         mid.point.x2 = mid_point[2:6])

df_distance <- sqrt((df_segment$mid.point.x2 - df_segment$mid.point.x1)^2 + (df_segment$y1 - df_segment$y2)^2)



df %>%
  ggplot()+
  geom_boxplot(mapping = aes(x = x.group,
                             y = y))+
  geom_point(data = df2, mapping = aes(x = x.group, y = med),
             shape = 4, color = 'red', size = 4)+
  geom_segment(data = df_segment,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               color = 'red')+
  scale_y_reverse()+
  theme_bw()+
  xlab('x-coordinate')+
  ylab('y-coordinate')+
  annotate("text", x = df2$x.group, y = df2$med-50, label = round(c(0, cumsum(df_distance)), 0))
