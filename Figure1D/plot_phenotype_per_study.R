##################################################
## Project: DGRPool
## Script purpose: Plotting Histogram of phenotypes per study
## Version: 1.0.0
## Date Created: 2023 Mar 14
## Date Modified: 2023 Mar 21
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

## Working directory
setwd("DGRPool/")

# Libraries
suppressPackageStartupMessages(library(data.table))

## Load phenotype data (run download_phenotypes.R script first)
data.all_pheno <- readRDS(file = "RDS/data.all_pheno_21_03_23_filtered.rds")

## Results
allpheno <- c(paste0(colnames(data.all_pheno[["F"]]), "_F")[2:ncol(data.all_pheno[["F"]])], paste0(colnames(data.all_pheno[["M"]]), "_M")[2:ncol(data.all_pheno[["M"]])], paste0(colnames(data.all_pheno[["NA"]]), "_NA")[2:ncol(data.all_pheno[["NA"]])])
data.annot_studies <- data.frame(limma::strsplit2(allpheno, "_"))
colnames(data.annot_studies) <- c("Study", "Id", "Sex")

## Draw Histogram, By sex
data.annot_studies$Study <- factor(data.annot_studies$Study, levels = paste0("S", 1:122)) 
p <- ggplot(data.annot_studies, aes(x = Study, fill = Sex)) + geom_histogram(stat="count") + ylab("Number of phenotypes") + xlab("Study id") + theme(axis.line.y = element_line(size = 0.3, colour = "black"), axis.line.x = element_line(size = 0.3, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, size = 6))
p <- p + scale_fill_manual(values = c("#f28c28", "#0e4c92", "#cccc77"))
ggsave(plot = p, filename = "NbPhenoPerStudy.plot.colored.by.sex.pdf", width = 8, height = 4)
ggsave(plot = p, filename = "NbPhenoPerStudy.plot.colored.by.sex.png", width = 8, height = 4)

## Draw Histogram
data.annot_studies <- data.annot_studies[,c(1,2)]
data.annot_studies <- data.annot_studies[!duplicated(data.annot_studies),]
p <- ggplot(data.annot_studies, aes(x = Study)) + geom_histogram(stat="count") + ylab("Number of phenotypes") + xlab("Study id") + theme(axis.line.y = element_line(size = 0.3, colour = "black"), axis.line.x = element_line(size = 0.3, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, size = 6))
ggsave(plot = p, filename = "NbPhenoPerStudy.plot.pdf", width = 8, height = 4)
ggsave(plot = p, filename = "NbPhenoPerStudy.plot.png", width = 8, height = 4)
