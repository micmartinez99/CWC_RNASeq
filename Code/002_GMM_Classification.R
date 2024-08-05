# The purpose of this script is to stratify patients into either low, medium, or high urolithin
# based on a Gaussian Mixture modelling model. 
# To start, a PCA analysis will be conducted to identify which urolithin metabolite is driving variation
# between samples.

# Install ggbiplot package
devtools::install_github("vqv/ggbiplot")

# Load libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(PCAtools)
library(ggfortify)
library(mclust)

# Initialize a function to create output directories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# Initialize an output directory
newDir("Outputs/002_GMM_Classification")

################################################################################

# Read in the 9 metabolite delta levels file generated in script 001_Tidy_Urolithins
deltas <- read.csv("Outputs/001_Tidy_Urolithins/All_Delta_Values.csv")
deltas <- deltas %>%
  column_to_rownames(var = "X")

# Remove the taotal urolithin column
total_uros <- deltas$Total.Urolithin
deltas$Total.Urolithin <- NULL

# Transpose the data
deltas <- as.data.frame(t(deltas))

# Run Principal component analysis
pca_res <- prcomp(deltas,
                  center = TRUE, 
                  scale. = TRUE)

# Plot a biplot of the data
biplot <- ggbiplot(pca_res, obs.scale = 1, var.scale = 1, 
        ellipse = TRUE, label.label = colnames(deltas),
        circle = TRUE) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 20),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.text.y = element_text(size = 16))
ggsave("Outputs/002_GMM_Classification/Urolithin_Delta_Biplot.png", biplot, width = 8, height = 8)

# Interpretation of the biplot:
# Direction of vectors: represent the original variables (metabolites). The direction
# indicates the direction in which that metabolite increases
# Length: The length of an arrow shows the magnitude of the variable's contribution to the principal
# components. Longer arrows mean the variable contributes more to the PCs.
# Angle between arrows: Small angles mean the variables are positively correlated.
# right angle: variables are uncorrelated
# Opposite directions: variables with arrows pointing in opposite directions are negatively correlated.

# Based on this, the main cause of variation between the cluster on the right and the cluster on the left
# seems to be urolithin A, isourolithin A, Urolithin C, Urolithin M5 and Urolithin M7.

################################################################################

# Cluster data based on full urolithin delta profile

# Set seed for reproducibilit
set.seed(03061999)

# Fi the GMM model
gmm_all <- Mclust(deltas[,c(2)], G = 2)
summary(gmm_all)

# Extract classifications
classifications <- gmm_all$classification
print(classifications)
deltas$Class <- classifications

# Plot the results
plot(gmm_all, what = "classification")

# Extract PCA scores
scores <- as.data.frame(pca_res$x)
scores$Group <- classifications
scores$Group <- ifelse(scores$Group == 1, "Low", "High")
scores$Group <- factor(scores$Group, levels = c("Low", "High"))


# Plot PCA colored by groupings
ggplot(scores, aes(x = PC1, y = PC2, color= Group)) +
  geom_point() +
  labs(y = "PC2 (18.2% explained var.)",
       x = "PC1 (38.0% explained var.)") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 20),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.text.y = element_text(size = 16))
  

# Compute the median values for each group
group1 <- apply(deltas[classifications == 1,], 2, median)
group2 <- apply(deltas[classifications == 2,], 2, median)
group3 <- apply(deltas[classifications == 3,], 2, median)




