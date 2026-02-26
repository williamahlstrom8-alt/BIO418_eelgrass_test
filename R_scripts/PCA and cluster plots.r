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

grp <- find.clusters(snp_zostera, n.pca=50, max.n=50,scale=FALSE)

grp$grp %>% as.vector

dapc <- dapc(snp_zostera, pop=grp$grp, n.pca=50)

# Scatter plot #

scatter(dapc, xax=1, yax=2, grp=dapc$grp, 
        col=transp(c("forestgreen","dodgerblue4","deeppink","orange2")),
        pch=19, bg="white",cstar = 1, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)

# Save this plot for later #
dapc_plot <- recordPlot()

## Plotting PCA again ####

s.class(pca$scores, fac=dapc$grp, clab=1,
        col = transp(funky(length(unique(dapc$grp)))), 
        csta=1, cpoint=2, cellipse =1, xax=1, yax=2)

## Testar om det funkar att skapa en karta för geografisk distribution ####

# Gör en tabell med populationer och kluster #

data2 <- data.frame(
  Population = pop(snp_zostera),
  Group = grp$grp
)

# Räkna antal individer per kluster per lokal #

library(dplyr)

count <- data2 %>%
  count(Population, Group)

# --- 1. Skapa latlong-tabell ---
latlong <- data.frame(
  Population = c("KOD","STE","GOT","GRO","HOG","ALA","YST",
                 "FUR","KAR","KLI","HOR","SLI","KRA","NYN","BJO"),
  
  Lat = c(58.882538,58.05143,57.39626,56.64116,56.19690,
          55.940683,55.4215587861172,56.0953999996185,56.969566,57.401162,
          57.695151,57.713318,58.6906470,58.8811079,59.83484288645309),
  
  Long = c(10.991989,11.81013,12.02123,12.78014,12.55141,
           12.771034,13.8463615898543,14.7202999750773,16.912377,18.152719,
           16.730625,18.8311672,17.4658689,17.9542217,19.078679980494496),
  
  stringsAsFactors = FALSE
)

# --- 2. Säkerställ samma datatyp ---
latlong$Population <- as.character(latlong$Population)
count$Population   <- as.character(count$Population)

# --- 3. Join ---
cluster_per_site <- inner_join(latlong, count, by="Population")

# --- 4. Kontroll ---
head(cluster_per_site)
unique(latlong$Population)
unique(count$Population)
names(cluster_per_site)
# Gör objekt för pie charts #

library(mapplots)

xyz <- make.xyz(
  cluster_per_site$Long,
  cluster_per_site$Lat,
  cluster_per_site$n,
  cluster_per_site$Group
)

# Ladda shapefile över Sverige #

install.packages(c("sf","rnaturalearth","rnaturalearthdata"))

library(sf)
library(rnaturalearth)
library(mapplots)

shape <- ne_countries(country="Sweden", returnclass="sf")

plot(st_geometry(shape), col="grey85")

cols <- c("forestgreen",
          "dodgerblue4",
          "deeppink",
          "orange2",
          "gold",
          "purple")

draw.pie(xyz$x, xyz$y, xyz$z,
         radius=0.3,
         col=cols)
