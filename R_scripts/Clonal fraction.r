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
snp_zostera <- readRDS("eelgrass_data/Manipulated_data/mlg_filtered_zostera.rds")
class(snp_zostera)

snp_zostrera_mlg <- readRDS

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


