# load packages and data ####
library("hierfstat")
library("pegas")
library("adegenet")
library("poppr")
library("tidyverse")
library("assignPOP")
library("radiator")

# load the genlight data
load("eelgrass_data/genlight_zostera_data.rda")

# Data manipulation ####
## convert to SNPclone ####
x <- as.snpclone(gl_zostera_25)
x
head(mll(x, "original"), 20)

mll(x) <- "original"
x

tab(x)

nmll(x) # count the number of multilocus genotypes

mll(x)  # show the multilocus genotype definitions

(xt <- apply(tab(x), 1, paste, collapse = ""))

# don't know what this does
# I think KOD-04 would mean individual number 4 from KOD
# And I think the number shown underneath/after is the MLG 
rank(xt, ties.method = "first")

# Clone correction ####
nmll(x, "original")

nmll(x, "contracted")

mll(x) <- "contracted"

x %>%
  clonecorrect(strata = NA) %>%  # 1. clone correct whole data set
  dist() %>%                     # 2. calculate distance
  sum()                          # 3. take the sum of the distance

set.seed(999)
x[sample(nInd(x))] %>% # 1. shuffle samples
  clonecorrect(strata = NA) %>%  # 2. clone correct whole data set
  dist() %>%                     # 3. calculate distance
  sum()                          # 4. take the sum of the distance

set.seed(1000)
x[sample(nInd(x))] %>% # 1. shuffle samples
  clonecorrect(strata = NA) %>%  # 2. clone correct whole data set
  dist() %>%                     # 3. calculate distance
  sum()                          # 4. take the sum of the distance

bruvo.dist(x, replen = c(1, 1))

# After help from Maru ####
# identify names, some are replicates
indNames(gl_zostera_25)

# create a matrix of distances
eg_mat <- bitwise.dist(
  gl_zostera_25,
  percent = TRUE,
  mat = TRUE,
  missing_match = TRUE,
  scale_missing = FALSE,
  euclidean = FALSE,
  differences_only = FALSE,
  threads = 0L
)

# returns largest distance 
max(eg_mat)

## Filter to only use duplicated samples ####
# Get individual names
ids <- indNames(gl_zostera_25)

# Remove "-rep" to get base sample names
base_ids <- sub("-rep$", "", ids)

# Identify samples that have replicates (appear more than once)
dup_ids <- base_ids[duplicated(base_ids) | duplicated(base_ids, fromLast = TRUE)]

# Keep only those samples (both original and replicate)
gl_zostera_rep <- gl_zostera_25[base_ids %in% dup_ids]

# Check result
indNames(gl_zostera_rep)

# matrix of distances in duplications
eg_mat_rep <- bitwise.dist(
  gl_zostera_rep,
  percent = TRUE,
  mat = TRUE,
  missing_match = TRUE,
  scale_missing = FALSE,
  euclidean = FALSE,
  differences_only = FALSE,
  threads = 0L
)

##Replicate dist####

###Identify replicates####
ids <- indNames(gl_zostera_25)

# Remove "-rep" suffix
base_ids <- sub("-rep$", "", ids)

# Keep only duplicated base IDs
dup_ids <- base_ids[duplicated(base_ids) | duplicated(base_ids, fromLast = TRUE)]

gl_rep <- gl_zostera_25[base_ids %in% dup_ids]

###Compute distance####
library(adegenet)

dist_mat <- dist(as.matrix(gl_rep))

###Transform distance to be between 0 and 1####
mat <- as.matrix(gl_rep)

rep_distances_scaled <- sapply(rep_pairs, function(pair) {
  ind1 <- mat[pair[1], ]
  ind2 <- mat[pair[2], ]
  
  # count loci that are not NA in both
  valid <- !is.na(ind1) & !is.na(ind2)
  
  # proportion of different genotypes
  mean(ind1[valid] != ind2[valid])
})

rep_distances_scaled

summary(rep_distances_scaled)

#create a minimum distance
#Anything smaller will be collapsed
zost_dist <- max(rep_distances_scaled)

#convert genlight to snpclone
snp_zostera <- as.snpclone(gl_zostera_25)

#mlg filter to our data
mlg.filter(snp_zostera, distance = eg_mat, threshold = zost_dist + .Machine$double.eps^0.5, stats = c("mlg", "thresholds"))
snp_zostera


## determine threshold####
plot.phylo(upgma(eg_mat))
eg_mat

zost_filtered <- filter_stats(snp_zostera, distance = bitwise.dist, plot = TRUE)

# different ways of calculating thresholds

print(farthest_thresh <- cutoff_predictor(zost_filtered$farthest$THRESHOLDS))

#shortest. Use one of the other
print(average_thresh  <- cutoff_predictor(zost_filtered$average$THRESHOLDS))

print(nearest_thresh  <- cutoff_predictor(zost_filtered$nearest$THRESHOLDS))


###MLG filter####
mlg.filter(snp_zostera, distance = bitwise.dist, algorithm = "f") <- farthest_thresh
snp_zostera
mlg(snp_zostera)

nmll(snp_zostera)

snp_zostera

zost_table <- mlg.table(snp_zostera, strata = ~Region/Site)
