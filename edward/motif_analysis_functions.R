library(ggrepel)
library(tidyverse)

#Input = a raw neuron projection data file containing the regions of interest
#In context of EC paper, data-file must already be split by MEC or LEC

generate_ma_data <- function(raw_data,regions_of_interest){
  
  #Extract the Regions Of Interest
  region_counts = raw_data[c(regions_of_interest)]
  
  #Binarize the Data
  binary_counts = as.data.frame(apply(region_counts, 2, function(x) {ifelse(x>=1, 1, 0)})) #converts every value in every column to either a 1 or a zero
  
  #Add Target Numbers
  binary_counts$target_no = apply(binary_counts,1,sum) #target no = the sum of each row
  
  #Add target numbers to raw data
  region_counts$target_no = binary_counts$target_no
  
  #Remove Zero Target Barcodes
  filtered_counts = region_counts %>% filter(target_no > 0)
  
  #Return Projection Data
  final_counts = filtered_counts[c(regions_of_interest)]
  
  return(final_counts)
  
}

