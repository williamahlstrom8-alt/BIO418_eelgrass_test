#load packages
library('adegenet')
library('poppr')

# load the genlight data
load("eelgrass_data/genlight_zostera_data.rda")

#convert genlight to snpclone
snp_zostera <- as.snpclone(gl_zostera_25)

#Save snp_zostera in the folder for manipulated data 
saveRDS(snp_zostera, file = "eelgrass_data/Manipulated_data/snp_zostera.rds")



