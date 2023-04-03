##################################################
## Project: DGRPool
## Script purpose: Pubmed barplot
## Version: 1.0.0
## Date Created: 2023 Mar 03
## Date Modified: 2023 Mar 03
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Working directory
setwd("/data/gardeux/DGRPool/")

# Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# Read Pubmed file
data.pubmed <- fread("PubMed_Timeline_Results_by_Year.csv", sep = ",", header = T, data.table = F)
	
# Plot
p <- ggplot(data.pubmed, aes(x = as.character(Year), y = Count)) + geom_col() + ylab("Number of publications per year") + xlab("") + theme(axis.line.y = element_line(linewidth = 0.3, colour = "black"), axis.line.x = element_line(linewidth = 0.3, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, size = 10))
ggsave(plot = p, filename = "Pubmed.plot.pdf", width = 3, height = 5)
