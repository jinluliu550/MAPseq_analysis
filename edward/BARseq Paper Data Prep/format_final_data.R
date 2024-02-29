library(ggrepel)
library(tidyverse)

#load in data (barseq brain projection matrices)
brain_c9 = read.csv("brain_c9.csv",header = FALSE)
brain_c28 = read.csv("brain_c28.csv", header = FALSE)

#load MAPseq brain
brain_c14 = read.csv("brain_c14.csv", header = FALSE)

#Read in headers
column_names = read.csv("column_names.csv")
header = colnames(column_names)

#Set headers
colnames(brain_c9) = header
colnames(brain_c28) = header
colnames(brain_c14) = header

#Remove negative control column (OB)
brain_c9 = subset(brain_c9, select = -c(1))
brain_c28 = subset(brain_c28, select = -c(1))
brain_c14 = subset(brain_c14, select = -c(1))

#Save files
write.csv(brain_c9,"BRAIN_C9.csv", row.names = FALSE)
write.csv(brain_c28,"BRAIN_C28.csv",row.names = FALSE)
write.csv(brain_c14,"BRAIN_C14.csv", row.names = FALSE)