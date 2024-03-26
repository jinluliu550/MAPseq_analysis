Hans_data1 <- read.csv('./data/Han-data/han_brain4_gh.csv')
Hans_data2 <- read.csv('./data/Han-data/han_brain5_gh.csv')
Hans_data3 <- read.csv('./data/Han-data/han_brain6_gh.csv')

data_Hans <- list(t(Hans_data1),
                  t(Hans_data2),
                  t(Hans_data3))
