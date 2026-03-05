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

best_run <- which.min(cross.entropy(proj, K = 3))

qmatrix <- Q(proj, K = 3, run = best_run)


# Plot admixture-diagram#

barplot(t(qmatrix),
        col = rainbow(ncol(qmatrix)),
        border = NA,
        xlab = "Individuals",
        ylab = "Ancestry proportion")

# Test för att sortera individer efter population #

ord <- order(pop(zostera))

barplot(t(qmatrix[ord, ]),
        col = rainbow(ncol(qmatrix)),
        border = NA,
        space = 0,
        xlab="Individuals",
        ylab="Ancestry")

# Test för att göra snyggare plot #

library(tidyverse)
library(forcats)

# Gör data frame av Q-matrix #

q_df <- as.data.frame(qmatrix)

# Lägg till metadata (ID + population) #

q_df$ID <- indNames(zostera_gl)
q_df$Population <- pop(zostera_gl)

# Gör long format #

q_long <- q_df %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Cluster",
    values_to = "Prob"
  )

# Sortera individer inom population #

q_long <- q_long %>%
  arrange(Population, ID)

# Test för att försöka sortera individer ancestry komponent #

q_long <- q_long %>%
  group_by(ID) %>%
  mutate(order = Cluster[which.max(Prob)]) %>%
  ungroup() %>%
  arrange(Population, order)

# Rita plot #

anc_plot_3K <- ggplot(q_long,
                   aes(factor(ID), Prob, fill = factor(Cluster))) +
  geom_col(width=1) +
  facet_grid(~fct_inorder(as.factor(Population)),
             switch="x",
             scales="free",
             space="free") +
  theme_minimal() + 
  labs(x="",
       y="Ancestry proportion",
       fill="Cluster") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=expansion(add=1)) +
  theme(panel.spacing.x=unit(0.05,"lines"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        strip.text=element_text(size=7))

anc_plot_3K

# Lägg till samma färger som tidigare analyser #

cols <- c("forestgreen","dodgerblue4","deeppink" )
anc_plot_3KR <- anc_plot_3K + scale_fill_manual(values = cols)

anc_plot_3KR

anc_plot_3KR <- recordPlot()

## K = 2 ####

# Extrahera ancestry proportions #

best_run <- which.min(cross.entropy(proj, K = 2))

qmatrix <- Q(proj, K = 2, run = best_run)


# Plot admixture-diagram#

barplot(t(qmatrix),
        col = rainbow(ncol(qmatrix)),
        border = NA,
        xlab = "Individuals",
        ylab = "Ancestry proportion")

# Test för att sortera individer efter population #

ord <- order(pop(zostera))

barplot(t(qmatrix[ord, ]),
        col = rainbow(ncol(qmatrix)),
        border = NA,
        space = 0,
        xlab="Individuals",
        ylab="Ancestry")

# Test för att göra snyggare plot #

library(tidyverse)
library(forcats)

# Gör data frame av Q-matrix #

q_df <- as.data.frame(qmatrix)

# Lägg till metadata (ID + population) #

q_df$ID <- indNames(zostera_gl)
q_df$Population <- pop(zostera_gl)

# Gör long format #

q_long <- q_df %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Cluster",
    values_to = "Prob"
  )

# Sortera individer inom population #

q_long <- q_long %>%
  arrange(Population, ID)

# Test för att försöka sortera individer ancestry komponent #

q_long <- q_long %>%
  group_by(ID) %>%
  mutate(order = Cluster[which.max(Prob)]) %>%
  ungroup() %>%
  arrange(Population, order)

# Rita plot #

anc_plot_2K <- ggplot(q_long,
                      aes(factor(ID), Prob, fill = factor(Cluster))) +
  geom_col(width=1) +
  facet_grid(~fct_inorder(as.factor(Population)),
             switch="x",
             scales="free",
             space="free") +
  theme_minimal() + 
  labs(x="",
       y="Ancestry proportion",
       fill="Cluster") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=expansion(add=1)) +
  theme(panel.spacing.x=unit(0.05,"lines"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        strip.text=element_text(size=7))

anc_plot_2K

# Lägg till samma färger som tidigare analyser #

cols <- c("forestgreen","deeppink" )
anc_plot_2KR <- anc_plot_2K + scale_fill_manual(values = cols)

anc_plot_2KR

anc_plot_2KR <- recordPlot()

## K = 7 ####

# Extrahera ancestry proportions #

best_run <- which.min(cross.entropy(proj, K = 7))

qmatrix <- Q(proj, K = 7, run = best_run)


# Plot admixture-diagram #

barplot(t(qmatrix),
        col = rainbow(ncol(qmatrix)),
        border = NA,
        xlab = "Individuals",
        ylab = "Ancestry proportion")

# Test för att sortera individer efter population #

ord <- order(pop(zostera))

barplot(t(qmatrix[ord, ]),
        col = rainbow(ncol(qmatrix)),
        border = NA,
        space = 0,
        xlab="Individuals",
        ylab="Ancestry")

# Test för att göra snyggare plot #

library(tidyverse)
library(forcats)

# Gör data frame av Q-matrix #

q_df <- as.data.frame(qmatrix)

# Lägg till metadata (ID + population) #

q_df$ID <- indNames(zostera_gl)
q_df$Population <- pop(zostera_gl)

# Gör long format #

q_long <- q_df %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Cluster",
    values_to = "Prob"
  )

# Sortera individer inom population #

q_long <- q_long %>%
  arrange(Population, ID)

# Test för att försöka sortera individer ancestry komponent #

q_long <- q_long %>%
  group_by(ID) %>%
  mutate(order = Cluster[which.max(Prob)]) %>%
  ungroup() %>%
  arrange(Population, order)

# Rita plot #

anc_plot_7K <- ggplot(q_long,
                      aes(factor(ID), Prob, fill = factor(Cluster))) +
  geom_col(width=1) +
  facet_grid(~fct_inorder(as.factor(Population)),
             switch="x",
             scales="free",
             space="free") +
  theme_minimal() + 
  labs(x="",
       y="Ancestry proportion",
       fill="Cluster") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=expansion(add=1)) +
  theme(panel.spacing.x=unit(0.05,"lines"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid=element_blank(),
        strip.text=element_text(size=7))

anc_plot_7K

# Lägg till samma färger som tidigare analyser #

cols <- c("forestgreen","dodgerblue4","deeppink","orange2","gold","purple", "red" )
anc_plot_7KR <- anc_plot_7K + scale_fill_manual(values = cols)

anc_plot_7KR

# Save this plot for later #

anc_plot_7KR <- recordPlot()

#Take the 3 plots we saved previously, and plot them together

# Ladda paketen #

library(gridExtra)
library(gridGraphics)
library(grid)

# Konvertera dina plots till grobs #

g1 <- grid.grabExpr(replayPlot(anc_plot_2KR))
g2 <- grid.grabExpr(replayPlot(anc_plot_3KR))
g3 <- grid.grabExpr(replayPlot(anc_plot_7KR))

# Arrangerar plotterna #

grid.arrange(
  g1,
  g2,
  g3,
  ncol = 3
)

# Ovanpå #

grid.arrange(
  g1,
  g2,
  g3,
  ncol = 1
)
