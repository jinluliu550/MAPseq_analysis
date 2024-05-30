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




# Brain dataset 
all_brains_set_sub <- read_excel("data/EC8_new/all_brains_set_sub.xlsx")
all_brains_set_sub$id <- factor(1:nrow(all_brains_set_sub),
                                levels = 1:nrow(all_brains_set_sub))

# List of brain regions
list_of_regions <- colnames(all_brains_set_sub)[1:11]
cortical_regions <- list_of_regions[1:8]

# Data frame
df00 <- all_brains_set_sub %>%
  pivot_longer(cols = 1:11,
               names_to = 'brain_region',
               values_to = 'barcode_counts') %>%
  mutate(brain_region = factor(brain_region,
                               levels = list_of_regions)) %>%
  group_by(id) %>%
  mutate(Region_max = brain_region[which.max(rank(barcode_counts[1:8], ties.method = 'random'))]) %>%
  ungroup() %>%
  group_by(Region_max,
           EC,
           brain_region) %>%
  summarise(mean_barcode_count = mean(barcode_counts),
            lower_interval = t.test(barcode_counts)$conf.int[1],
            upper_interval = t.test(barcode_counts)$conf.int[2]) %>%
  ungroup()


# Plot
df00 %>%
  ggplot(mapping = aes(x = brain_region,
                       y = mean_barcode_count,
                       fill = EC))+
  geom_bar(
           stat = 'identity',
           position = 'dodge')+
  geom_errorbar(aes(ymin = lower_interval,
                    ymax = upper_interval),
                alpha = 0.9,
                width = 0.2,
                stat = 'identity',
                position = position_dodge(width = 0.9))+
  facet_wrap(~Region_max, nrow = 3, scales = 'free')+
  theme_bw()+
  xlab('brain region')+
  ylab('barcode counts')+
  scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))


# Flow graph
df01 <- all_brains_set_sub %>%
  pivot_longer(cols = 1:11,
               names_to = 'brain_region',
               values_to = 'barcode_counts') %>%
  mutate(brain_region = factor(brain_region,
                               levels = list_of_regions)) %>%
  group_by(id) %>%
  mutate(Region_max = brain_region[which.max(rank(barcode_counts[1:8], ties.method = 'random'))]) %>%
  ungroup() %>%
  mutate(projection = ifelse(barcode_counts == 0, 'No', 'Yes')) %>%
  group_by(EC,
           brain_region,
           Region_max) %>%
  summarise(Freq = length(which(projection == 'Yes'))) %>%
  ungroup()

ggplot(data = df01,
       aes(y = Freq,
           axis1 = Region_max,
           axis2 = brain_region))+
  geom_flow(aes(fill = EC))+
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_x_discrete(limits = c("Strongest Projecting Regions", "All Regions"))+
  theme_bw()+
  scale_fill_manual(values=c(LEC = '#16697A', MEC = '#DB6400'))