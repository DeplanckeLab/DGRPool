##################################################
## Project: DGRPool
## Script purpose: Figure 4A
## Version: 1.0.0
## Date Created: 2024 Aug 05
## Date Modified: 2024 Aug 05
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

## Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2)) # for plotting
suppressPackageStartupMessages(library(dplyr)) # for fast filtering / summarizing

# Load object from file (~1mn30s)
data.gwas <- fread("../resources/data.gwas_24_08_02.filtered.tsv.gz", sep = "\t", data.table = F)
rownames(data.gwas) <- data.gwas$V1
data.gwas <- data.gwas[,-1]

# Retrieve all variants significant at p <= 1E-6
cutoff <- 1E-6
significant_counts <- summarise(data.gwas[,8:ncol(data.gwas)], across(everything(), ~sum(. <= cutoff, na.rm = TRUE)))

# Plot the histogram
df <- data.frame(Phenotype = names(significant_counts), SignificantVariants = unlist(significant_counts[1,]))
df$SignificantVariants[df$SignificantVariants > 50] <- 50
p <- ggplot(df, aes(x = SignificantVariants)) + 
  geom_histogram(bins = 50, fill = "#6A73B4") +
  labs(x = bquote("Nb significant GWAS variants / phenotype at" ~ p <= 1 %*% 10^-6), y = "Frequency") +
  geom_hline(yintercept = -10, color = "black") +
  theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), panel.grid.minor = element_line(color = "lightgrey"), panel.grid.major = element_line(color = "lightgrey"), panel.background = element_blank(), plot.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line.y = element_line(color = "black"), plot.title = element_text(hjust = 0.5, size = 8)) +
  scale_x_continuous(breaks = seq(0, 50, by = 10), labels = c("0", "10", "20", "30", "40", ">50")) +
  scale_y_continuous(breaks = seq(0, 800, by = 50), expand = c(0, 0))

ggsave(plot = p, filename = "Figure_4A.pdf", width = 3, height = 2.5)
