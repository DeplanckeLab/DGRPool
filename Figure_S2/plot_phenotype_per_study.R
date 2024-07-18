##################################################
## Project: DGRPool
## Script purpose: Plotting Histogram of phenotypes per study
## Version: 1.0.0
## Date Created: 2023 Mar 14
## Date Modified: 2024 Jul 18
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

## Load phenotype data (run download_phenotypes.R script first)
data.all_pheno <- readRDS(file = "../RDS/data.all_pheno_15_07_24_filtered.rds")

## Results
allpheno <- c(paste0(colnames(data.all_pheno[["F"]]), "_F")[2:ncol(data.all_pheno[["F"]])], paste0(colnames(data.all_pheno[["M"]]), "_M")[2:ncol(data.all_pheno[["M"]])], paste0(colnames(data.all_pheno[["NA"]]), "_NA")[2:ncol(data.all_pheno[["NA"]])])
data.annot_studies <- data.frame(limma::strsplit2(allpheno, "_"))
colnames(data.annot_studies) <- c("Study", "Id", "Sex")
data.annot_studies$Study <- factor(data.annot_studies$Study, levels = paste0("S", 1:max(as.numeric(gsub(x = data.annot_studies$Study, pattern = "S", replacement = ""))))) 

## Format phenotypes
phenotypes_information <- fread("../resources/phenotypes.tsv", data.table = F)
phenotypes_information$study_id <- paste0("S", phenotypes_information$study_id)
phenotypes_information <- phenotypes_information[,c("study_id", "study_name")]
phenotypes_information <- phenotypes_information[!duplicated(phenotypes_information),]
phenotypes_information$study_name <- paste0("(",phenotypes_information$study_name,")")
rownames(phenotypes_information) <- phenotypes_information$study_id
phenotypes_information <- phenotypes_information[levels(data.annot_studies$Study),]
phenotypes_information <- subset(phenotypes_information, !is.na(study_name))
data.annot_studies$Study_Name <- phenotypes_information[as.character(data.annot_studies$Study), "study_name"]

## Draw Histogram
data.annot_studies <- data.annot_studies[,c(1,2,4)]
data.annot_studies <- data.annot_studies[!duplicated(data.annot_studies),]
p <- ggplot(data.annot_studies, aes(x = Study)) + geom_bar() + ylim(c(0,100)) + ylab("Number of phenotypes") + xlab("Study id") + theme(axis.line.y = element_line(linewidth = 0.3, colour = "black"), axis.line.x = element_line(linewidth = 0.3, colour = "black"), panel.grid.major.y = element_line(linewidth = 0.5, colour = "lightgrey"), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, size = 6, hjust=0.95, vjust=0.2))
ggsave(plot = p, filename = "NbPhenoPerStudy.plot.pdf", width = 8, height = 5)

## Draw Histogram (x axis as REF study)
p <- ggplot(data.annot_studies, aes(x = Study)) + geom_bar() + scale_x_discrete(name = "", breaks = phenotypes_information$study_id, labels = phenotypes_information$study_name) + ylim(c(0,100)) + ylab("Number of phenotypes") + theme(axis.line.y = element_line(linewidth = 0.3, colour = "black"), axis.line.x = element_line(linewidth = 0.3, colour = "black"), panel.grid.major.y = element_line(linewidth = 0.5, colour = "lightgrey"), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, size = 6, hjust=0.95, vjust=0.2))
ggsave(plot = p, filename = "NbPhenoPerStudy.plot.ref.pdf", width = 8, height = 5)
