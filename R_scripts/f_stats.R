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

##Multivariate analsis####
x.eg <- tab(gl_zostera_no_rep, freq=TRUE, NA.method="mean")
pca.eg <- dudi.pca(x.eg, center=TRUE, scale=FALSE, scannf = FALSE, nf = 3)

pca.eg

###basic PCA####
s.class(pca.eg$li, fac=pop(gl_zostera_no_rep), col=funky(15))

###more advanced####
s.class(pca.eg$li, fac=pop(gl_zostera_no_rep), col=transp(funky(15),.6), axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.eg$eig[1:50],3,1,2, ratio=.3)

### 1st & 3rd axis####
s.class(pca.eg$li, fac=pop(gl_zostera_no_rep),
        xax=1, yax=3, col=transp(funky(15),.6),
        axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.eg$eig[1:50],3,1,3, ratio=.3)

eig.perc <- 100*pca.eg$eig/sum(pca.eg$eig)
head(eig.perc)

##diveRsity trial####
library(diveRsity)
library(radiator)
new_ind_name <- paste0(gl_zostera_no_rep$pop, 1:length(gl_zostera_no_rep$pop))
zost_renamed <- `indNames<-` (gl_zostera_no_rep, new_ind_name)

# Time to write the genepopfile
genomic_converter(data = zost_renamed, output = c("genepop"), filename = "eelgrass")

##StrataG####
#devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)
library("strataG") #OBS!!! UNCOMMENT ABOVE IF THIS DOESN'T WORK

#Overall test
ovt_zost <- overallTest(genlight2gtypes(gl_zostera_no_rep),by.locus=F,nrep=1000) #all loci at once

#Pairwise test
pwt_zost <- pairwiseTest(genlight2gtypes(gl_zostera_no_rep),by.locus=F,nrep=1000) #all loci at once

pwt_summary <- pairwiseSummary(pwt_zost)
pwt_summary

pwt_summary_fst <- pwt_summary %>% 
  dplyr::select(label,strata.1,strata.2,n.1,n.2,Fst,Fst_p.val) %>% 
  mutate(Fst_p.val_FDR=p.adjust(Fst_p.val))
pwt_summary_fst

