# load packages
library("hierfstat")
library("pegas")
library("adegenet")
library("poppr")
library("tidyverse")
library("assignPOP")
library("radiator")

# load the genlight data
load("eelgrass_data/genlight_zostera_data.rda")

View(gl_zostera_25)

# view number of individuals per location
table(pop(gl_zostera_25))

f
