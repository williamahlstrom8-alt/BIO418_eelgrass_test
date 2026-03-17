#Packages and data####
library(readr)
goodCF <- read_csv("eelgrass_data/Manipulated_data/goodCF.csv")

#t-test, split = GOT####
#split into two groups, Saline = KOD-GOT, Brackish = GRÖ-BJO
goodCF$group <- ifelse(goodCF$Population %in% c("KOD","STE","GOT"),
                       "saline", "brackish")

#run t-test
t.test(ClonalFraction ~ group, data = goodCF, alternative = "two.sided")

#create boxplot 
boxplot(ClonalFraction ~ group, data = goodCF)

#Reload data to remove group
goodCF <- read_csv("eelgrass_data/Manipulated_data/goodCF.csv")
