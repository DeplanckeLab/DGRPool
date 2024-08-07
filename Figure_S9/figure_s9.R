##################################################
## Project: DGRPool
## Script purpose: Figure 4a
## Version: 1.0.0
## Date Created: 2024 Aug 05
## Date Modified: 2024 Aug 05
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

## Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2)) # for plotting
suppressPackageStartupMessages(library(ggpubr)) # for arranging plots in grid
suppressPackageStartupMessages(library(dplyr)) # for fast filtering / summarizing

# Load object from file (~1mn30s)
data.gwas <- fread("../resources/data.gwas_24_08_02.filtered.tsv.gz", sep = "\t", data.table = F)
rownames(data.gwas) <- data.gwas$V1
data.gwas <- data.gwas[,-1]

# All plots
p_list <- list()

# Retrieve all variants significant at p <= 0.05
cutoff <- 0.05
significant_counts <- summarise(data.gwas[,8:ncol(data.gwas)], across(everything(), ~sum(. <= cutoff, na.rm = TRUE)))

# Plot the histogram
df <- data.frame(Phenotype = names(significant_counts), SignificantVariants = unlist(significant_counts[1,]))
p_list[[as.character(cutoff)]] <- ggplot(df, aes(x = SignificantVariants)) + 
  geom_histogram(bins = 60, fill = "#6A73B4") +
  labs(title = bquote(p <= 0.05), x = bquote("Nb significant GWAS variants / phenotype at" ~ p <= 0.05), y = "Frequency") +
  geom_hline(yintercept = -10, color = "black") +
  theme(plot.margin = unit(c(2,2,2,2), 'lines'), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), panel.grid.minor = element_line(color = "lightgrey"), panel.grid.major = element_line(color = "lightgrey"), panel.background = element_blank(), plot.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line.y = element_line(color = "black"), plot.title = element_text(hjust = 0.5, size = 8)) +
  scale_x_continuous(breaks = seq(0, 160000, by = 20000), labels = c("0", "20k", "40k", "60k", "80k", "100k", "120k", "140k", "160k")) +
  scale_y_continuous(breaks = seq(0, 200, by = 20), expand = c(0, 0))

# Retrieve all variants significant at p <= 0.01
cutoff <- 0.01
significant_counts <- summarise(data.gwas[,8:ncol(data.gwas)], across(everything(), ~sum(. <= cutoff, na.rm = TRUE)))

# Plot the histogram
df <- data.frame(Phenotype = names(significant_counts), SignificantVariants = unlist(significant_counts[1,]))
p_list[[as.character(cutoff)]] <- ggplot(df, aes(x = SignificantVariants)) + 
  geom_histogram(bins = 60, fill = "#6A73B4") +
  labs(title = bquote(p <= 0.01), x = bquote("Nb significant GWAS variants / phenotype at" ~ p <= 0.01), y = "Frequency") +
  geom_hline(yintercept = -10, color = "black") +
  theme(plot.margin = unit(c(2,2,2,2), 'lines'), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), panel.grid.minor = element_line(color = "lightgrey"), panel.grid.major = element_line(color = "lightgrey"), panel.background = element_blank(), plot.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line.y = element_line(color = "black"), plot.title = element_text(hjust = 0.5, size = 8)) +
  scale_x_continuous(breaks = seq(0, 160000, by = 20000), labels = c("0", "20k", "40k", "60k", "80k", "100k", "120k", "140k", "160k")) +
  scale_y_continuous(breaks = seq(0, 200, by = 20), expand = c(0, 0))

# Retrieve all variants significant at p <= 1E-3
cutoff <- 1E-3
significant_counts <- summarise(data.gwas[,8:ncol(data.gwas)], across(everything(), ~sum(. <= cutoff, na.rm = TRUE)))

# Plot the histogram
df <- data.frame(Phenotype = names(significant_counts), SignificantVariants = unlist(significant_counts[1,]))
p_list[[as.character(cutoff)]] <- ggplot(df, aes(x = SignificantVariants)) + 
  geom_histogram(bins = 60, fill = "#6A73B4") +
  labs(title = bquote(p <= 1 %*% 10^-3), x = bquote("Nb significant GWAS variants / phenotype at" ~ p <= 1 %*% 10^-3), y = "Frequency") +
  geom_hline(yintercept = -10, color = "black") +
  theme(plot.margin = unit(c(2,2,2,2), 'lines'), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), panel.grid.minor = element_line(color = "lightgrey"), panel.grid.major = element_line(color = "lightgrey"), panel.background = element_blank(), plot.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line.y = element_line(color = "black"), plot.title = element_text(hjust = 0.5, size = 8)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = c("0", "5k", "10k", "15k", "20k", "25k", "30k")) +
  scale_y_continuous(breaks = seq(0, 800, by = 50), expand = c(0, 0))

# Retrieve all variants significant at p <= 1E-4
cutoff <- 1E-4
significant_counts <- summarise(data.gwas[,8:ncol(data.gwas)], across(everything(), ~sum(. <= cutoff, na.rm = TRUE)))

# Plot the histogram
df <- data.frame(Phenotype = names(significant_counts), SignificantVariants = unlist(significant_counts[1,]))
df$SignificantVariants[df$SignificantVariants > 5000] <- 5000
p_list[[as.character(cutoff)]] <- ggplot(df, aes(x = SignificantVariants)) + 
  geom_histogram(bins = 60, fill = "#6A73B4") +
  labs(title = bquote(p <= 1 %*% 10^-4), x = bquote("Nb significant GWAS variants / phenotype at" ~ p <= 1 %*% 10^-4), y = "Frequency") +
  geom_hline(yintercept = -10, color = "black") +
  theme(plot.margin = unit(c(2,2,2,2), 'lines'), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), panel.grid.minor = element_line(color = "lightgrey"), panel.grid.major = element_line(color = "lightgrey"), panel.background = element_blank(), plot.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line.y = element_line(color = "black"), plot.title = element_text(hjust = 0.5, size = 8)) +
  scale_x_continuous(breaks = seq(0, 5000, by = 1000), labels = c("0", "1k", "2k", "3k", "4k", ">5k")) +
  scale_y_continuous(breaks = seq(0, 800, by = 50), expand = c(0, 0))

# Retrieve all variants significant at p <= 1E-5
cutoff <- 1E-5
significant_counts <- summarise(data.gwas[,8:ncol(data.gwas)], across(everything(), ~sum(. <= cutoff, na.rm = TRUE)))

# Plot the histogram
df <- data.frame(Phenotype = names(significant_counts), SignificantVariants = unlist(significant_counts[1,]))
df$SignificantVariants[df$SignificantVariants > 1000] <- 1000
p_list[[as.character(cutoff)]] <- ggplot(df, aes(x = SignificantVariants)) + 
  geom_histogram(bins = 60, fill = "#6A73B4") +
  labs(title = bquote(p <= 1 %*% 10^-5), x = bquote("Nb significant GWAS variants / phenotype at" ~ p <= 1 %*% 10^-5), y = "Frequency") +
  geom_hline(yintercept = -10, color = "black") +
  theme(plot.margin = unit(c(2,2,2,2), 'lines'), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), panel.grid.minor = element_line(color = "lightgrey"), panel.grid.major = element_line(color = "lightgrey"), panel.background = element_blank(), plot.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line.y = element_line(color = "black"), plot.title = element_text(hjust = 0.5, size = 8)) +
  scale_x_continuous(breaks = seq(0, 1000, by = 200), labels = c("0", "200", "400", "600", "800", ">1000")) +
  scale_y_continuous(breaks = seq(0, 800, by = 100), expand = c(0, 0))

# Retrieve all variants significant at p <= 1E-6
cutoff <- 1E-6
significant_counts <- summarise(data.gwas[,8:ncol(data.gwas)], across(everything(), ~sum(. <= cutoff, na.rm = TRUE)))

# Plot the histogram
df <- data.frame(Phenotype = names(significant_counts), SignificantVariants = unlist(significant_counts[1,]))
df$SignificantVariants[df$SignificantVariants > 200] <- 200
p_list[[as.character(cutoff)]] <- ggplot(df, aes(x = SignificantVariants)) + 
  geom_histogram(bins = 60, fill = "#6A73B4") +
  labs(title = bquote(p <= 1 %*% 10^-6), x = bquote("Nb significant GWAS variants / phenotype at" ~ p <= 1 %*% 10^-6), y = "Frequency") +
  geom_hline(yintercept = -10, color = "black") +
  theme(plot.margin = unit(c(2,2,2,2), 'lines'), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), panel.grid.minor = element_line(color = "lightgrey"), panel.grid.major = element_line(color = "lightgrey"), panel.background = element_blank(), plot.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line.y = element_line(color = "black"), plot.title = element_text(hjust = 0.5, size = 8)) +
  scale_x_continuous(breaks = seq(0, 200, by = 25), labels = c("0", "25", "50", "75", "100", "125", "150", "175", ">200")) +
  scale_y_continuous(breaks = seq(0, 800, by = 100), expand = c(0, 0))

# Retrieve all variants significant at p <= 1E-7
cutoff <- 1E-7
significant_counts <- summarise(data.gwas[,8:ncol(data.gwas)], across(everything(), ~sum(. <= cutoff, na.rm = TRUE)))

# Plot the histogram
df <- data.frame(Phenotype = names(significant_counts), SignificantVariants = unlist(significant_counts[1,]))
df$SignificantVariants[df$SignificantVariants > 50] <- 50
p_list[[as.character(cutoff)]] <- ggplot(df, aes(x = SignificantVariants)) + 
  geom_histogram(bins = 50, fill = "#6A73B4") +
  labs(title = bquote(p <= 1 %*% 10^-7), x = bquote("Nb significant GWAS variants / phenotype at" ~ p <= 1 %*% 10^-7), y = "Frequency") +
  geom_hline(yintercept = -10, color = "black") +
  theme(plot.margin = unit(c(2,2,2,2), 'lines'), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), panel.grid.minor = element_line(color = "lightgrey"), panel.grid.major = element_line(color = "lightgrey"), panel.background = element_blank(), plot.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line.y = element_line(color = "black"), plot.title = element_text(hjust = 0.5, size = 8)) +
  scale_x_continuous(breaks = seq(0, 50, by = 10), labels = c("0", "10", "20", "30", "40", ">50")) +
  scale_y_continuous(breaks = seq(0, 800, by = 100), expand = c(0, 0))

# Retrieve all variants significant at p <= bonferroni (~2.67E-8)
cutoff <- 0.05 / nrow(data.gwas)
message("Bonferroni threshold = ", formatC(cutoff, format = "e", digits = 2))
significant_counts <- summarise(data.gwas[,8:ncol(data.gwas)], across(everything(), ~sum(. <= cutoff, na.rm = TRUE)))

# Plot the histogram
df <- data.frame(Phenotype = names(significant_counts), SignificantVariants = unlist(significant_counts[1,]))
df$SignificantVariants[df$SignificantVariants > 30] <- 30
p_list[[as.character(cutoff)]] <- ggplot(df, aes(x = SignificantVariants)) + 
  geom_histogram(bins = 30, fill = "#6A73B4") +
  labs(title = bquote(p <= 2.67 %*% 10^-8 ~ " (Bonferroni)"), x = bquote("Nb significant GWAS variants / phenotype at" ~ p <= 2.67 %*% 10^-8), y = "Frequency") +
  geom_hline(yintercept = -10, color = "black") +
  theme(plot.margin = unit(c(2,2,2,2), 'lines'), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), panel.grid.minor = element_line(color = "lightgrey"), panel.grid.major = element_line(color = "lightgrey"), panel.background = element_blank(), plot.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line.y = element_line(color = "black"), plot.title = element_text(hjust = 0.5, size = 8)) +
  scale_x_continuous(breaks = seq(0, 30, by = 5), labels = c("0", "5", "10", "15", "20", "25", ">30")) +
  scale_y_continuous(breaks = seq(0, 800, by = 100), expand = c(0, 0))

p <- ggarrange(plotlist = p_list, ncol = 2, nrow = 4)
ggsave(plot = p, filename = "Supp.Figure.S9.pdf", width = 7, height = 10)
