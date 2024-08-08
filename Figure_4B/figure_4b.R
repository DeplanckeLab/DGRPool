##################################################
## Project: DGRPool
## Script purpose: Figure 4B
## Version: 1.0.0
## Date Created: 2024 Aug 05
## Date Modified: 2024 Aug 08
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

## Libraries
suppressPackageStartupMessages(library(data.table)) # For fast reading of tsv
suppressPackageStartupMessages(library(ggplot2)) # For plotting
suppressPackageStartupMessages(library(Cairo)) # For saving pdf in good resolution
suppressPackageStartupMessages(library(ramwas)) # Fast QQ-Plot & Manhattan plot

# Load object from file (~1mn30s)
data.gwas <- fread("../resources/data.gwas_24_08_02.filtered.tsv.gz", sep = "\t", data.table = F)
rownames(data.gwas) <- data.gwas$V1
data.gwas <- data.gwas[,-1]

# Retrieve all variants significant at p <= 1E-6
cutoff <- 1E-6

# Count values <= cutoff per row (Takes ~22s)
significant_counts <- rowSums(data.gwas[, 8:ncol(data.gwas)] <= cutoff, na.rm = TRUE)

## We display the results as a Manhattan plot
# Prepare
chroms <- factor(data.gwas$chr)
levels(chroms) <- c("2L", "2R", "3L", "3R", "X", "4")
m <- manPlotPrepare(pvalues = significant_counts, ismlog10 = T, chr=chroms, pos = data.gwas$pos)

# Top hits
message("Variants significant in 12 different phenotypes: ", paste0(names(significant_counts)[significant_counts == 12], collapse = ", "))
message("Variants significant in 11 different phenotypes: ", paste0(names(significant_counts)[significant_counts == 11], collapse = ", "))
message("Variants significant in 10 different phenotypes: ", paste0(names(significant_counts)[significant_counts == 10], collapse = ", "))
message("Variants significant in 9 different phenotypes: ", paste0(names(significant_counts)[significant_counts == 9], collapse = ", "))
message("Variants significant in 8 different phenotypes: ", paste0(names(significant_counts)[significant_counts == 8], collapse = ", "))
message("Variants significant in 7 different phenotypes: ", paste0(names(significant_counts)[significant_counts == 7], collapse = ", "))

# Plot in Pdf
CairoPDF("Figure_4B.pdf", height = 2.8, width = 7)
manPlotFast(man = m, lwd = 1, colorSet = c("black", "darkgrey"), cex = 1, yaxmax = 20)
dev.off()
