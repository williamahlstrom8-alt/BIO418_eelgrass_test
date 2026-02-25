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



