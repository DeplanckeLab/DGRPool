##################################################
## Project: DGRPool
## Script purpose: Plotting phenotype-phenotype correlation heatmap for some selected phenotypes
## Version: 1.0.0
## Date Created: 2023 Mar 22
## Date Modified: 2024 Jul 24
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

## Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggplot2))

## Load correlation matrix

# Load pre-computed data
data.correlation_spearman <- fread("../resources/phenotype_correlation_spearman.tsv", header = T, sep = "\t", data.table = F)
rownames(data.correlation_spearman) <- data.correlation_spearman[,1]
data.correlation_spearman <- data.correlation_spearman[,-1]
all.phenotypes <- colnames(data.correlation_spearman)
message("Found ", length(all.phenotypes), " phenotypes in ALL studies")

## Load p-values for finding significant overlap

# Load pre-computed data
data.p_spearman <- fread("../resources/phenotype_correlation_spearman_pvalue.tsv", header = T, sep = "\t", data.table = F)
rownames(data.p_spearman) <- data.p_spearman[,1]
data.p_spearman <- data.p_spearman[,-1]

# Phenotypes of interest: mn_longevity (Arya et al, 2010), i.e. phenotypes S1_1315_F and S1_1315_M
data.p_value_F <- unlist(data.p_spearman["S1_1315_F", ])
data.p_value_F <- data.p_value_F[!is.na(data.p_value_F)] # Remove NAs (outliers / not calculated)
#fdr_F <- p.adjust(data.p_value_F, method = "fdr")
#sig.F <- names(fdr_F)[fdr_F <= 0.05]

data.p_value_M <- unlist(data.p_spearman["S1_1315_M", ])
data.p_value_M <- data.p_value_M[!is.na(data.p_value_M)] # Remove NAs (outliers / not calculated)
#fdr_M <- p.adjust(data.p_value_M, method = "fdr")
#sig.M <- names(fdr_M)[fdr_M <= 0.05]

#significant_phenotypes <- sort(unique(c(sig.F, sig.M)))
#significant_phenotypes

fdr_both <- p.adjust(c(data.p_value_F, data.p_value_M), method = "fdr")
sig.both <- names(fdr_both)[fdr_both <= 0.05]
significant_phenotypes <- sort(unique(sig.both))
significant_phenotypes_id <- unique(limma::strsplit2(significant_phenotypes, "_")[,2])

message(length(significant_phenotypes_id), " FDR 5% significant phenotypes were found")
message(length(significant_phenotypes), " FDR 5% significant sex-specific phenotypes were found")

## Format phenotypes

# Retrieve phenotype information (description)
phenotypes_information <- fread("../resources/phenotypes.tsv", data.table = F)
rownames(phenotypes_information) <- phenotypes_information$phenotype_id
phenotypes_information <- subset(phenotypes_information, phenotype_id %in% significant_phenotypes_id)

# Format
formatted.phenotypes <- c()
for(phenotype in significant_phenotypes){
  pheno <- strsplit(phenotype, split = "_")[[1]]
  formatted.phenotypes <- c(formatted.phenotypes, paste0("(", phenotypes_information[pheno[2],"study_name"], ")[",pheno[3],"] ", phenotypes_information[pheno[2],"phenotype_name"], " - ", phenotypes_information[pheno[2],"description"]))
}

## Spearman's correlation heatmap

# Prepare data
data.heatmap <- data.correlation_spearman[significant_phenotypes, significant_phenotypes]
rownames(data.heatmap) <- substr(formatted.phenotypes, start = 1, stop = 110) #formatted.phenotypes
colnames(data.heatmap) <- substr(formatted.phenotypes, start = 1, stop = 110) #formatted.phenotypes

# Draw heatmap and generate pdf
pdf("correlation_heatmap_spearman.pdf", width = 20, height = 10)
suppressWarnings(p <- pheatmap(data.heatmap, clustering_method = "ward.D2", width = 20, height = 10, legend = T, scale = "none", color = colorRampPalette(c("blue", "white", "red"))(100), na_col="white", border_color ="white", show_colnames = F, show_rownames = T, silent = T))
suppressWarnings(gridExtra::grid.arrange(p$gtable, vp=viewport(width=0.9, height=1)))
dev.off()

