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

head(gl_zostera_25)

gl_zostera_25@pop

# view number of individuals per location
table(pop(gl_zostera_25))

data(gl_zostera_25)


# create SNPclone 
snp_zostera <- as.snpclone(gl_zostera_25)
snp_zostera
view(snp_zostera)
head(mll(snp_zostera, "original"), 20)

# change mll to custom
mll(x) <- "custom"
x

head(mll(x, "custom"), 20) # Showing the definitions for the first 20 samples

# Naive thing ####
mll(x) <- "original"
x
## Grid example ####
grid_example <- matrix(c(1, 4,
                         1, 1,
                         5, 1,
                         9, 1,
                         9, 4),
                       ncol = 2,
                       byrow = TRUE)
rownames(grid_example) <- LETTERS[1:5]
colnames(grid_example) <- c("x", "y")
grid_example

x1 <- as.genclone(df2genind(grid_example, ploidy = 1))
tab(x1)  # Look at the multilocus genotype table

nmll(x1) # count the number of multilocus genotypes
mll(x1)  # show the multilocus genotype definitions

x <- as.genclone(df2genind(rbind(grid_example, new = c(5, NA)), ploidy = 1))
tab(x)  # Note the missing data at locus 2.

nmll(x)
mll(x)

grid_new <- rbind(grid_example,
                  new = c(5, NA),
                  mut = c(5, 2)
)
x <- as.genclone(df2genind(grid_new, ploidy = 1))
tab(x)

(xt <- apply(tab(x), 1, paste, collapse = ""))

rank(xt, ties.method = "first")

library("ape")
raw_dist <- function(x){
  dist(genind2df(x, usepop = FALSE))
}
(xdis <- raw_dist(x))

plot.phylo(upgma(xdis))

# collapses MLGs that are close to each other
mlg.filter(x, distance = xdis) <- 1 + .Machine$double.eps^0.5

x

mll(x) <- "original" # We'll reset to the naive definition

mll(x, "original")

mlg.filter(x, distance = xdis, threshold = 1, stats = c("mlg", "thresholds"))

(e <- .Machine$double.eps^0.5) # A very tiny number

mlg.filter(x, distance = xdis, threshold = 1 + e, stats = c("mlg", "thresholds"))

# threshold of 0 allows extremely similat MLGs to remain separate
mlg.filter(x, distance = xdis, threshold = 0, stats = c("mlg", "thresholds"))

mll(x, "original")

x

mlg.table(x) # Before: 7 MLGs

mlg.filter(x, distance = xdis) <- 1 + e
x

mlg.table(x) # After: 5 MLGs

mlg.filter(x) <- 4.51
x

mlg.table(x)

# DNAGER DANGER WILL CAUSE ERROR!
rm(xdis) # NOOOOOO!
try(mlg.filter(x) <- 1 + e)

# solution: raw_dist, we defined earlier
mlg.filter(x, distance = raw_dist) <- 1 + e
x

# Bruvo dist much safer
bruvo.dist(x, replen = c(1, 1))

mlg.filter(x, distance = bruvo.dist, replen = c(1, 1)) <- 0.44
x

# access original 
mll(x, "original")

# access the collapsed mll's
mll(x) # contracted

# replace collapsed data with original 
mll(x) <- "original"
mll(x) # original

# Load the data set "Pinf"
data(Pinf)
Pinf

pinfreps <- fix_replen(Pinf, c(2, 2, 6, 2, 2, 2, 2, 2, 3, 3, 2))
pinf_filtered <- filter_stats(Pinf, distance = bruvo.dist, replen = pinfreps, plot = TRUE)

### Determining thresholds ####
# Determine Threshold for farthest
print(farthest_thresh <- cutoff_predictor(pinf_filtered$farthest$THRESHOLDS))

# Threshold for average 
print(average_thresh  <- cutoff_predictor(pinf_filtered$average$THRESHOLDS))

# Threshold for nearest
print(nearest_thresh  <- cutoff_predictor(pinf_filtered$nearest$THRESHOLDS))

## Define MLL using following criteria: ####
# [t] Threshold = 0.113 (see largest threshold)
# [d] Distance = Bruvo's distance 
# [a] Algorithm = Farthest Neighbour (from Threshold)

mlg.filter(Pinf, distance = bruvo.dist, replen = pinfreps, algorithm = "f") <- farthest_thresh
Pinf

# Custom method, for parial clones ####
data(partial_clone)
pc <- as.genclone(partial_clone)
mll(pc)
pc

LETTERS[mll(pc)]  # The new MLGs

mll.custom(pc) <- LETTERS[mll(pc)]
# Shows 4 subplots, one for each population
mlg.table(pc)

pcpal <- colorRampPalette(c("blue", "gold"))
set.seed(9001)
pcmsn <- bruvo.msn(pc, replen = rep(1, nLoc(pc)), palette = pcpal,
                   vertex.label.color = "firebrick", vertex.label.font = 2,
                   vertex.label.cex = 1.5)
# If we have evidence that Q = M we can do this to collapse them
mll.levels(pc)[mll.levels(pc) == "Q"] <- "M"

# PLot again with collapsed Q and M
set.seed(9001)
pcmsn <- bruvo.msn(pc, replen = rep(1, nLoc(pc)), palette = pcpal,
                   vertex.label.color = "firebrick", vertex.label.font = 2,
                   vertex.label.cex = 1.5)

# Diversity analysis ####
data(monpop)
splitStrata(monpop) <- ~Tree/Year/Symptom
montab <- mlg.table(monpop, strata = ~Symptom/Year)

# returns basic diversity stats
(monstat <- diversity_stats(montab))

# confidence intervals of the stats
diversity_ci(montab, n = 100L, raw = FALSE)

?diversity_ci

# Add your own custom statistics 
myCF <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){ # if it's a matrix
    res <- rowSums(x > 0)/rowSums(x)
  } else {                 # if it's a vector
    res <- sum(x > 0)/sum(x)
  }
  return(res)
}
(monstat2 <- diversity_stats(montab, CF = myCF))


# Repeat lengths are necessary
reps <- fix_replen(monpop,
                   c(CHMFc4 = 7, CHMFc5 = 2, CHMFc12 = 4,
                     SEA = 4, SED = 4, SEE = 2, SEG = 6,
                     SEI = 3, SEL = 4, SEN = 2,
                     SEP = 4, SEQ = 2, SER = 4))

# Adding a little bit, so the threshold is included.
e <- .Machine$double.eps^0.5

# Using the default farthest neighbor algorithm to collapse genotypes
mlg.filter(monpop, distance = bruvo.dist, replen = reps) <- (0.5/13) + e
montabf <- mlg.table(monpop, strata = ~Symptom/Year)