##################################################
## Project: DGRPool
## Script purpose: Plotting phenotype-phenotype correlation heatmap for some selected phenotypes
## Version: 1.0.0
## Date Created: 2023 Mar 22
## Date Modified: 2024 Jul 24
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggplot2))

## Load datasets
data.correlation_spearman <- fread("../resources/phenotype_correlation_spearman.tsv", header = T, sep = "\t", data.table = F)
rownames(data.correlation_spearman) <- data.correlation_spearman[,1]
data.correlation_spearman <- data.correlation_spearman[,-1]

## Order samples by study
data.annot_studies <- data.frame(limma::strsplit2(colnames(data.correlation_spearman), "_"))
colnames(data.annot_studies) <- c("Study", "Phenotype", "Sex")
rownames(data.annot_studies) <- colnames(data.correlation_spearman)

## Restrict to some phenotypes
interesting_phenotypes <- as.character(fread("listphenotypes.txt", data.table = F)$V1)
message(length(interesting_phenotypes), " phenotypes were loaded")
data.annot_studies <- subset(data.annot_studies, Phenotype %in% interesting_phenotypes)
data.correlation_spearman <- data.correlation_spearman[rownames(data.annot_studies),rownames(data.annot_studies)]

## Format phenotypes
phenotypes_information <- fread("../resources/phenotypes.tsv", data.table = F)
rownames(phenotypes_information) <- phenotypes_information$phenotype_id
phenotypes_information <- subset(phenotypes_information, phenotype_id %in% interesting_phenotypes)
message(nrow(phenotypes_information), " phenotypes found")

formatted.phenotypes <- c()
for(phenotype in rownames(data.annot_studies)){
	pheno <- strsplit(phenotype, split = "_")[[1]]
	formatted.phenotypes <- c(formatted.phenotypes, paste0("(", phenotypes_information[pheno[2],"study_name"], ")[",pheno[3],"] ", phenotypes_information[pheno[2],"phenotype_name"], " - ", phenotypes_information[pheno[2],"description"]))
}

## Draw Spearman's correlation heatmap, all, no thresholding
data.heatmap <- data.correlation_spearman
rownames(data.heatmap) <- formatted.phenotypes
colnames(data.heatmap) <- formatted.phenotypes
#png("correlation_heatmap_spearman.png", width = 2200, height = 2000)
#pheatmap(data.heatmap, legend = T, scale = "none", color = colorRampPalette(c("blue", "white", "red"))(100), na_col="white", border_color ="white")
#dev.off()

#pdf("correlation_heatmap_spearman.pdf", width = 22, height = 20)
pheatmap(data.heatmap, legend = T, scale = "none", color = colorRampPalette(c("blue", "white", "red"))(100), na_col="white", border_color ="white")
#dev.off()

