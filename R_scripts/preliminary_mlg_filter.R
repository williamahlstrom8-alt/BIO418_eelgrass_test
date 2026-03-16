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
snp_zostera@pop

#Split west/east-coast####
# Define population groups
group1 <- c("KOD", "STE", "GOT", "GRO", "HOG", "ALA")
group2 <- c("YST", "FUR", "KAR", "KLI", "HOR", "SLI", "KRA", "NYN", "BJO")

# Subset the snpclone object
snp_zostera_west <- snp_zostera[pop(snp_zostera) %in% group1]
snp_zostera_east <- snp_zostera[pop(snp_zostera) %in% group2]
#check the levels in the two new object
snp_zostera_west@pop
snp_zostera_east@pop

#MLG filtering####
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

###MLG filter Farthest####
#apply the filter
mlg.filter(snp_zostera_west, distance = bitwise.dist, algorithm = "f") <- farthest_thresh

snp_zostera_west
#returns number of MLGs 
mlg(snp_zostera_west)

#plots number of individuals over each MLG for each region
zost_table <- mlg.table(snp_zostera_west)

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
#Farthest threshold, should be used (Largest value)
print(farthest_thresh <- cutoff_predictor(zost_filtered_east$farthest$THRESHOLDS))

#Average threshold, should NOT be used (Smallest value)
print(average_thresh  <- cutoff_predictor(zost_filtered_east$average$THRESHOLDS))

#Nearest threshold, same as largest for some reason (idk why)
print(nearest_thresh  <- cutoff_predictor(zost_filtered_east$nearest$THRESHOLDS))

###MLG filter Farthest####
#apply the filter
mlg.filter(snp_zostera_east, distance = bitwise.dist, algorithm = "f") <- farthest_thresh

snp_zostera_east
#returns number of MLGs 
mlg(snp_zostera_east)

#plots number of individuals over each MLG for each region
zost_table <- mlg.table(snp_zostera_east)

#Combine west/east-coast snpsclones####
class(snp_zostera_east)
class(snp_zostera_west)


#create matrices 
#simply merging the snpobjects didn't work because they
#were missing something in their classes
mat_west <- as.matrix(snp_zostera_west)
mat_east <- as.matrix(snp_zostera_east)

#combine matrices
mat_combined <- rbind(mat_west, mat_east)

#create a genlight object which we can later convert to snpclone
zostera_genlight <- new("genlight", mat_combined)
snp_zostera_west@loc.all
#restore assignments 
pop(zostera_genlight) <- c(pop(snp_zostera_west), pop(snp_zostera_east))
ploidy(zostera_genlight) <- c(ploidy(snp_zostera_west), ploidy(snp_zostera_east))
pop(zostera_genlight) <- c(pop(snp_zostera_west), pop(snp_zostera_east))
#strata didn't work unless snp_zostera was used. Don't know why
strata(zostera_genlight) <- strata(snp_zostera)


snp_zostera_new <- as.snpclone(zostera_genlight)

#Clone correction####
library(poppr)

# Remove clones (keep one per MLG)
snp_zostera_new@pop
snpclone_cc <- clonecorrect(snp_zostera_new, strata = ~Site)
snpclone_cc@pop

mlg(snpclone_cc)
popNames(snpclone_cc) <- popNames(snp_zostera_new)
# Check new number of individuals
nInd(snpclone_cc)


##Save no clone data####
saveRDS(snpclone_cc, file = "eelgrass_data/Manipulated_data/clone_corrected_zostera.rds")

#MLG table with Clonal fraction####
##Import CF####
library(readr)
goodCF <- read_csv("eelgrass_data/Manipulated_data/goodCF.csv")
View(goodCF)

##add CF to pop for mlg.table####
snp_zostera@pop
#convert dataframe to vector
clonal_fraction <- setNames(goodCF$ClonalFraction, goodCF$Population)

# Extract current levels
old_levels <- levels(snp_zostera@pop)

# Create new labels with rounded values (3 decimals looks nice)
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

# MLG filter farthest####
#apply the filter
mlg.filter(snp_zostera, distance = bitwise.dist, algorithm = "f") <- farthest_thresh

snp_zostera
#returns number of MLGs 
mlg(snp_zostera)

#plots number of individuals over each MLG for each region
zost_table <- mlg.table(snp_zostera, strata = ~Region/Site)
