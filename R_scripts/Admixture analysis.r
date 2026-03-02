#### Admixture analysis ####

# Loading file ####

zostera <- readRDS("/Users/rasmusgreen/Desktop/BIO418_eelgrass_test/eelgrass_data/Manipulated_data/clone_corrected_zostera.rds")

class(zostera)

# Korrigerad pipeline #

library(poppr)
library(adegenet)

zostera_gl <- as(zostera, "genlight")

class(zostera_gl)

# Installera och ladda LEA (om ej gjort) #

install.packages("BiocManager")
BiocManager::install("LEA")
library(LEA)

# Konvertera genlight → LEA geno-format #

geno_matrix <- as.matrix(zostera_gl)

dim(geno_matrix)

table(geno_matrix, useNA="ifany")

# Om det finns NA → ersätt dem #
geno_matrix[is.na(geno_matrix)] <- 9

write.geno(geno_matrix, "zostera.geno")

# Kör admixture-analys (snmf) #

proj <- snmf("zostera.geno",
             K = 1:10,
             entropy = TRUE,
             repetitions = 10,
             project = "new")


# Hitta optimalt K #

ce <- sapply(1:10, function(k)
  min(cross.entropy(proj, K = k)))

plot(1:10, ce, type="b", pch=19,
     xlab="K", ylab="Cross-entropy")


# Extrahera ancestry proportions #




# Plot admixture-diagram#



