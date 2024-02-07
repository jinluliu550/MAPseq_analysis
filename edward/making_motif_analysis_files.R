library(ggrepel)
library(tidyverse)
source("motif_analysis_functions.R")

#Basis Pipeline = take the "all_brains" file for each Threshold value
#Split into MEC and LEC
#Take the MEC file
#Take only the first ten columns (neocortical regions)

#Read in data

#800
brains = read.csv("all_brains_800_maria.csv",header = TRUE)

#Double Check that all the IDs in brains are unique
length(brains$id) == length(unique(brains$id)) # if TRUE, all IDs are unique

#prepare raw data for function
#split by EC

brains$EC = as.factor(brains$EC)

split_brains = split(brains, brains$EC)

brains_MEC = split_brains$MEC

brains_LEC = split_brains$LEC

#Define regions of interest

fullSet = c("ORB", "PFC", "ACA", "SS", "PTL", "VIS", "RSC",
            "ECT","dStr","vStr","OLF")

set1 = c("ORB", "PFC", "ACA", "SS", "PTL", "VIS", "RSC")

newSet = c("PFC", "ORB", "SS", "AUD", "ACA", "PTL", "VIS", "RSC") # current focus

#Generate Motif Analysis files

#MEC Results

mec_set1_800 = generate_ma_data(brains_MEC,set1)

mec_fullSet_800 = generate_ma_data(brains_MEC,fullSet)

mec_newSet_800 = generate_ma_data(brains_MEC,newSet)

#LEC Results

lec_set1_800 = generate_ma_data(brains_LEC,set1)

lec_fullSet_800 = generate_ma_data(brains_LEC,fullSet)

lec_newSet_800 = generate_ma_data(brains_LEC,newSet)

#Save MEC Results

write.csv(mec_set1_800, "mec_set1_800.csv", row.names=FALSE)
write.csv(mec_fullSet_800, "mec_fullSet_800.csv", row.names=FALSE)
write.csv(mec_newSet_800, "mec_newSet_800.csv", row.names=FALSE)

#Save LEC Results

write.csv(lec_set1_800, "lec_set1_800.csv", row.names=FALSE)
write.csv(lec_fullSet_800, "lec_fullSet_800.csv", row.names=FALSE)
write.csv(lec_newSet_800, "lec_newSet_800.csv", row.names=FALSE)
