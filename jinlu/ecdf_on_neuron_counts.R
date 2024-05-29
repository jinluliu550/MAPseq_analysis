ecdf_counts <- data.frame(barcode_count = as.vector(do.call(cbind, EC8_new)),
                          region = rownames(EC8_new[[1]]),
                          EC = rep(unlist(EC8_EC_label), each = R))

ecdf_counts$region <- factor(ecdf_counts$region,
                             levels = rownames(EC8_new[[1]]))

ggplot(ecdf_counts, aes(barcode_count))+
  stat_ecdf(geom = 'step', mapping = aes(color = EC))+
  facet_wrap(~region, nrow = 2, scales = 'free')+
  scale_color_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))+
  theme_bw()+
  ylab('Cumulative Proportion')+
  xlab('Barcode count')+
  xlim(c(0,NA))





all_brains_set_sub <- read_excel("data/EC8_new/all_brains_set_sub.xlsx")


# EC8 new
digit.list <- list(one = 1,
                   two = 2,
                   three = 3,
                   four = 4,
                   five = 5,
                   six = 6)

EC8_subnew_list <- lapply(1:6, 
                          function(m) all_brains_set_sub %>%
                            filter(brain == names(digit.list)[m]))

EC8_sub_new <- lapply(1:6,
                      function(m) t(round(EC8_subnew_list[[m]][,1:11], 0)))

# List of brain regions
list_of_regions <- rownames(EC8_sub_new[[1]])
cortical_regions <- list_of_regions[1:8]


# Function to calculate the mean barcode counts
mean_barcode_counts <- function(brain_dataset){
  
  # Select the brain regions that you want to have included in the figure 
  # the order is important
  brain <- brain_dataset %>% select(list_of_regions)
  
  # CHANGE THE NUMBER BELOW TO MATCH THE NUMBER OF REGIONS YOU SELECTED
  number_of_regions <- 11
  # Remove neurons that might have no projections to any of the selected regions
  
  brain <- brain[apply(brain[, 1:number_of_regions], 1, function(x) any(x != 0)), ]
  
  
  # Make neurons' id a factor
  brain$id <- as.factor(brain$id)
  
  
  brain <- brain %>% gather("Region", "Count", 1:number_of_regions) %>% group_by(id) %>%
    mutate(Count_norm_area = Count / sum(Count))
  
  
  brain$Region <- factor(brain$Region, levels = list_of_regions)
  
  only_cortical <- brain %>% filter(Region %in% cortical_regions)
  
  only_cortical <- only_cortical %>% group_by(id) %>% mutate(Region_max = which.max(rank(Count, ties.method = "random")))
  
  maximum_region <- only_cortical %>% select(id, Region_max) %>% distinct()
  
  
  brain <- left_join(brain, maximum_region, by = "id")
  
  
  # If you want the function to return the dataset to calculate mean barcode count values 
  # For generating gel_plots, comment this out
  return(brain)
}



# Function to calculate the standard error of the mean and 95% confidence intervals
mean_int95_summary <- function(brain_dataset){
  brains_list <- c("one", "two", "three", "four")
  list_of_datasets <- NULL
  for (b in brains_list){
    data <- brain_dataset %>% filter(brain == b)
    data <- mean_barcode_counts(data)
    df_long <- data %>% group_by(Region_max, Region) %>% summarise(mean_barcode_count = mean(Count))
    df_long$Region_max <- with(df_long, factor(Region_max,
                                               levels = c('1', '2', '3', '4', '5', '6', '7', '8'), labels = cortical_regions))
    df_long$brain <- b
    list_of_datasets <- rbind(list_of_datasets, df_long)
  }
  
  mean_int95 <- list_of_datasets %>% group_by(Region_max, Region) %>%
    summarise(mean = mean(mean_barcode_count), int_95 = 1.96 * (sd(mean_barcode_count)/sqrt(4)), 
              sem = sd(mean_barcode_count)/sqrt(4))
  
  return(mean_int95)
}


MEC_800 <- brains %>% filter(EC == "MEC")
LEC_800 <- brains %>% filter(EC == "LEC")


MEC <- mean_int95_summary(MEC_800) %>% mutate(EC = "MEC")
LEC <- mean_int95_summary(LEC_800) %>% mutate(EC = "LEC")
MEC_LEC <- rbind(MEC, LEC)


table_of_values <- MEC_LEC %>% pivot_wider(names_from = Region, values_from = c(mean, int_95, sem))


barplot <- ggplot(MEC_LEC, aes(x = Region, y = mean, fill = EC, color = Region_max)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.90) +
  labs(x = "Region", y = "Proportion") +
  theme_classic() +
  scale_fill_manual(values = c("MEC" = "#5e9597","LEC" = "#f0bc69")) +
  theme(axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "right",
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) + 
  scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
  facet_wrap(~Region_max, scales = "free_y")

print(barplot)

# Plug in MEC_800 or LEC_800 (or any other dataset) in the function

only_cortical <- mean_barcode_counts(LEC_800)

df_long <- only_cortical %>% group_by(Region_max, Region) %>% summarise(mean_barcode_count = mean(Count))
df_long$Region_max <- with(df_long, factor(Region_max,
                                           levels = c('1', '2', '3', '4', '5', '6', '7', '8'), labels = cortical_regions))

mean_barcode_plots <- ggplot(df_long, aes(x = Region, y = mean_barcode_count)) +
  geom_line(aes(group = Region_max)) +
  labs(x = "", y = "Mean barcode count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

mean_barcode_plots_separate <- mean_barcode_plots + facet_wrap(~Region_max, scales = "free_y") + 
  theme(legend.position = "none")



