#Load packages and data####
#load packages
library('adegenet')
library('poppr')

# load the genlight data
load("eelgrass_data/genlight_zostera_data.rda")

#Deal with replicates####
#obtain individual names of samples (some are replicates)
indNames(gl_zostera_25)

##Identify and remove replicates####
# Identify individuals that are NOT replicates
keep <- !grepl("-rep$", indNames(gl_zostera_25))

# Subset the genlight object
gl_zostera_no_rep <- gl_zostera_25[keep, ]

#number of replicates
sum(grepl("-rep$", indNames(gl_zostera_25)))  

#number of sampled individuals (no replicates)
nInd(gl_zostera_no_rep)

#number of samples per region
table(pop(gl_zostera_no_rep))

#convert genlight to snpclone####
#we will use snpclone in the rest of the project
snp_zostera <- as.snpclone(gl_zostera_no_rep)

#Save snp_zostera in the folder for manipulated data 
saveRDS(snp_zostera, file = "eelgrass_data/Manipulated_data/snp_zostera.rds")

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

#### Test av möjlig PCA-plot ####

library(adegenet)
library(ade4)

# PCA
pca <- glPca(snp_zostera, nf = 3)

# Plot
s.class(pca$scores,
        fac = pop(snp_zostera),
        clab = 1,
        col = transp(funky(length(levels(pop(snp_zostera))))),
        csta = 0,
        cpoint = 4,
        cellipse = 1,
        xax = 1,
        yax = 2)

### Find clusters with the appropriate K-value ####

grp <- find.clusters(snp_zostera, n.pca=50, max.n=40,scale=FALSE)
