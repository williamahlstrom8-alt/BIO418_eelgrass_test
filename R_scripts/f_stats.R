#Load packages and data####
library(hierfstat)
library(adegenet)
library(dplyr)
library(viridis)
library(adegenet)
library(reshape2)
library(ggplot2)
library(ape)


load("eelgrass_data/genlight_zostera_data.rda")

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

#Fst calculations####
##Convert to data frame####
# example: convert genlight to a data frame
zost_mat <- as.matrix(gl_zostera_no_rep)
zost_df <- data.frame(pop = pop(gl_zostera_no_rep), zost_mat)

##basic Fst####
# now suitable for wc()
wc(zost_df)

basic.stats(zost_df)
##Pairwise Fst####
# run a pairwise fst
matFst <- pairwise.WCfst(zost_df) #this might take a while

# Here we replace all the NA values with 0 for easier plotting
matFst <- matFst %>% 
  replace(., is.na(.), 0)

###Fst table with color####
df <- melt(matFst)

ggplot(df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis(name = "Fst") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = NULL, y = NULL)

####magma colors####
ggplot(df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis(name = "Fst", option = "magma", direction = -1) +
  theme_minimal() +
  theme(text = element_text(size = 15),axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = NULL, y = NULL)

###Fst boxplot####
fst_eg <- matFst
diag(fst_eg) <- NA
boxplot(fst_eg, col=funky(15), ylab="Fst", las=2)

###Fst tree####
eg.tree <- nj(matFst)
plot(eg.tree, type="unr", tip.col=funky(nPop(gl_zostera_no_rep)), font=2)
add.scale.bar()

####replace popnames with numbers####
# change population names to numbers
pop_num <- as.numeric(gl_zostera_no_rep$pop)

#make a data frame with population numbers in the first column, followed by the data
eg_dat_num <- data.frame(pop_num,gl_zostera_no_rep[,-1])

# Performs bootstrapping over loci of pairwise Fst
fst_boot <- boot.ppfst(eg_dat_num) # This might take a while
# lower CI
fst_boot$ll
# upper CI
fst_boot$ul