#Load packages and data####
#packages
library("hierfstat")
library("pegas")
library("adegenet")
library("poppr")
library("tidyverse")
library("assignPOP")
library("radiator")

library(poppr)

sessionInfo()
#data
snp_zostera <- readRDS("eelgrass_data/Manipulated_data/snp_zostera.rds")
class(snp_zostera)

#plot nmll over genetic distance cutoff, and three threshold methods
zost_filtered <- filter_stats(snp_zostera, distance = bitwise.dist, plot = TRUE)

###different ways of calculating thresholds####
#Farthest threshold, should be used (Largest value)
print(farthest_thresh <- cutoff_predictor(zost_filtered$farthest$THRESHOLDS))

#Average threshold, should NOT be used (Smallest value)
print(average_thresh  <- cutoff_predictor(zost_filtered$average$THRESHOLDS))

#Nearest threshold, same as largest for some reason (idk why)
print(nearest_thresh  <- cutoff_predictor(zost_filtered$nearest$THRESHOLDS))

#MLG filter UPGMA####
#apply the filter
mlg.filter(snp_zostera, distance = bitwise.dist, algorithm = "a") <- average_thresh

snp_zostera
#returns number of MLGs 
mlg(snp_zostera)

# Skapa MLG tabell #

mlg_table <- poppr::mlg.table(snp_zostera, plot = FALSE)
dim(mlg_table)

# Definierar funktion #

myCF <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){
    res <- rowSums(x > 0) / rowSums(x)
  } else {
    res <- sum(x > 0) / sum(x)
  }
  return(res)
}

# Räkna ut klonal fraktion #

myCF(mlg_table)

# Klassisk klonal fraktion #

goodCF <- 1 - myCF(mlg_table)

class(goodCF)
str(goodCF)

##convert goodCF to data frame####
goodCF_df <- data.frame(
  Population = names(goodCF),
  ClonalFraction = goodCF
)

write.csv(goodCF_df, "eelgrass_data/Manipulated_data/goodCF.csv", row.names = FALSE)


