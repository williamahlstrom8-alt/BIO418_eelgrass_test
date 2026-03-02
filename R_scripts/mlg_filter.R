#Load packages and data####
#packages
library("hierfstat")
library("pegas")
library("adegenet")
library("poppr")
library("tidyverse")
library("assignPOP")
library("radiator")

#data
snp_zostera <- readRDS("eelgrass_data/Manipulated_data/snp_zostera.rds")

#Create distance matrix####
eg_mat <- bitwise.dist(
  snp_zostera,
  percent = TRUE,
  mat = TRUE,
  missing_match = TRUE,
  scale_missing = FALSE,
  euclidean = FALSE,
  differences_only = FALSE,
  threads = 0L
)

#determine threshold for collapsing####
plot.phylo(upgma(eg_mat))
eg_mat
#plot nmll over genetic distance cutoff, and three threshold methods
zost_filtered <- filter_stats(snp_zostera, distance = bitwise.dist, plot = TRUE)

###different ways of calculating thresholds####
#Farthest threshold, should be used (Largest value)
print(farthest_thresh <- cutoff_predictor(zost_filtered$farthest$THRESHOLDS))

#Average threshold, should NOT be used (Smallest value)
print(average_thresh  <- cutoff_predictor(zost_filtered$average$THRESHOLDS))

#Nearest threshold, same as largest for some reason (idk why)
print(nearest_thresh  <- cutoff_predictor(zost_filtered$nearest$THRESHOLDS))

#MLG filter####
#apply the filter
mlg.filter(snp_zostera, distance = bitwise.dist, algorithm = "a") <- average_thresh

snp_zostera

#returns number of MLGs 
mlg(snp_zostera)

#plots number of individuals over each MLG for each region
zost_table <- mlg.table(snp_zostera, strata = ~Region/Site)

#Clone correction####
library(poppr)

# Remove clones (keep one per MLG)
snpclone_cc <- clonecorrect(snp_zostera)

# Check new number of individuals
nInd(snpclone_cc)

##Save no clone data####
saveRDS(snpclone_cc, file = "eelgrass_data/Manipulated_data/clone_corrected_zostera.rds")
