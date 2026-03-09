#Packages and data####
library(readr)
goodCF <- read_csv("eelgrass_data/Manipulated_data/goodCF.csv")

#t-tests####
##split = ALA####
# exempel på t-test formel
#t.test(formula, data, subset, na.action, …)
# Create group belonging
goodCF$group <- ifelse(goodCF$Population %in% c("KOD","STE","GOT","GRO","HOG","ALA"),
                       "west", "east")

#run t-test
t.test(ClonalFraction ~ group, data = goodCF, alternative = "two.sided")

#create boxplot 
boxplot(ClonalFraction ~ group, data = goodCF)

#Reload data to remove group
goodCF <- read_csv("eelgrass_data/Manipulated_data/goodCF.csv")

##split = HOG####
goodCF$group <- ifelse(goodCF$Population %in% c("KOD","STE","GOT","GRO","HOG"),
                       "west", "east")

#run t-test
t.test(ClonalFraction ~ group, data = goodCF, alternative = "two.sided")

#create boxplot 
boxplot(ClonalFraction ~ group, data = goodCF)

#Reload data to remove group
goodCF <- read_csv("eelgrass_data/Manipulated_data/goodCF.csv")

##split = GOT####
goodCF$group <- ifelse(goodCF$Population %in% c("KOD","STE","GOT"),
                       "saline", "brackish")

#run t-test
t.test(ClonalFraction ~ group, data = goodCF, alternative = "two.sided")

#create boxplot 
boxplot(ClonalFraction ~ group, data = goodCF)

#Reload data to remove group
goodCF <- read_csv("eelgrass_data/Manipulated_data/goodCF.csv")
