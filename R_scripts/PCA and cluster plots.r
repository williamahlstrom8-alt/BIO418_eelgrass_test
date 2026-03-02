#Load packages and data####
#load packages
library('adegenet')
library('poppr')

# load the genlight data
load("eelgrass_data/genlight_zostera_data.rda")


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
snp_zostera <- readRDS("eelgrass_data/Manipulated_data/clone_corrected_zostera.rds")


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
