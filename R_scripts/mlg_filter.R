#Load packages and data####
#packages
library("hierfstat")
library("pegas")
library("adegenet")
library("poppr")
library("tidyverse")
library("assignPOP")
library("radiator")
library(readr)


#data
snp_zostera <- readRDS("eelgrass_data/Manipulated_data/snp_zostera.rds")

#Split west/east-coast####
#To determine thresholds for collapsing MLLs and clone correction
# Define population groups, 1 = West, 2 = East
group1 <- c("KOD", "STE", "GOT", "GRO", "HOG", "ALA")
group2 <- c("YST", "FUR", "KAR", "KLI", "HOR", "SLI", "KRA", "NYN", "BJO")

# Create the grous
snp_zostera_west <- snp_zostera[pop(snp_zostera) %in% group1]
snp_zostera_east <- snp_zostera[pop(snp_zostera) %in% group2]

##West-coast####
###Create distance matrises####
eg_mat_west <- bitwise.dist(
  snp_zostera_west,
  percent = TRUE,
  mat = TRUE,
  missing_match = TRUE,
  scale_missing = FALSE,
  euclidean = FALSE,
  differences_only = FALSE,
  threads = 0L
)

###determine threshold for collapsing####
plot.phylo(upgma(eg_mat_west))
eg_mat_west
#plot nmll over genetic distance cutoff, and three threshold methods
zost_filtered_west <- filter_stats(snp_zostera_west, distance = bitwise.dist, plot = TRUE)

###different ways of calculating thresholds####
#Farthest threshold, should be used (Largest value)
print(farthest_thresh <- cutoff_predictor(zost_filtered_west$farthest$THRESHOLDS))

#Average threshold, should NOT be used (Smallest value)
print(average_thresh  <- cutoff_predictor(zost_filtered_west$average$THRESHOLDS))

#Nearest threshold, same as largest for some reason (idk why)
print(nearest_thresh  <- cutoff_predictor(zost_filtered_west$nearest$THRESHOLDS))

# West coast gives a threshold of 0.01496726 from farthest_thresh

##East-coast####
###Create distance matrises####
eg_mat_east <- bitwise.dist(
  snp_zostera_east,
  percent = TRUE,
  mat = TRUE,
  missing_match = TRUE,
  scale_missing = FALSE,
  euclidean = FALSE,
  differences_only = FALSE,
  threads = 0L
)

###determine threshold for collapsing####
plot.phylo(upgma(eg_mat_east))
eg_mat_east
#plot nmll over genetic distance cutoff, and three threshold methods
zost_filtered_east <- filter_stats(snp_zostera_east, distance = bitwise.dist, plot = TRUE)

###different ways of calculating thresholds####
#These fail for east-coast for some reason
#Farthest threshold, should be used 
print(farthest_thresh <- cutoff_predictor(zost_filtered_east$farthest$THRESHOLDS))

#Average threshold, should NOT be used 
print(average_thresh  <- cutoff_predictor(zost_filtered_east$average$THRESHOLDS))

#Nearest threshold, same as largest for some reason
print(nearest_thresh  <- cutoff_predictor(zost_filtered_east$nearest$THRESHOLDS))

#Threshold calculations return wrong threshold, 
#looking at the zost_filtered_east plot reveals 0.016 as a good threshold


#MLG filter####
#Use original dataset and 0.016 threshold from East coast
#apply the filter
mlg.filter(snp_zostera, distance = bitwise.dist, algorithm = "a") <- 0.016

snp_zostera
#returns number of MLGs 
mlg(snp_zostera)

#plots number of individuals over each MLG for each region
zost_table <- mlg.table(snp_zostera)

#Save mlg filtered data####
saveRDS(snp_zostera, file = "eelgrass_data/Manipulated_data/mlg_filtered_zostera.rds")

#Clone correction####
# Remove clones (keep one per MLG)
snp_zostera@pop
snpclone_cc <- clonecorrect(snp_zostera, strata = ~Site)
snpclone_cc@pop

popNames(snpclone_cc) <- popNames(snp_zostera)
# Check new number of individuals
nInd(snpclone_cc)

##Save no clone data####
saveRDS(snpclone_cc, file = "eelgrass_data/Manipulated_data/clone_corrected_zostera.rds")

#MLG table with Clonal fraction####
##Import CF####
goodCF <- read_csv("eelgrass_data/Manipulated_data/goodCF.csv")
View(goodCF)

##add CF to pop for mlg.table####
snp_zostera@pop
#convert dataframe to vector
clonal_fraction <- setNames(goodCF$ClonalFraction, goodCF$Population)

# Extract current levels
old_levels <- levels(snp_zostera@pop)

#Create new pop levels with clonal fractions added (with 3 decimals)
new_levels <- paste0(
  old_levels,
  " (",
  round(clonal_fraction[old_levels], 3),
  ")"
)

levels(snp_zostera@pop) <- new_levels

levels(snp_zostera@pop)

#MLG table with clonal fraction
zost_table_cf <- mlg.table(snp_zostera)
