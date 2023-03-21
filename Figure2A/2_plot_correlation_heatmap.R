##################################################
## Project: DGRPool
## Script purpose: Plotting phenotype-phenotype correlation heatmap
## Version: 1.0.0
## Date Created: 2022 Dec 12
## Date Modified: 2023 Mar 20
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

## Working directory
setwd("/data/gardeux/DGRPool/")

# Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggplot2))

## Load datasets
# data.correlation_pearson <- fread("phenotype_correlation_pearson.tsv", header = T, sep = "\t", data.table = F)
# rownames(data.correlation_pearson) <- data.correlation_pearson[,1]
# data.correlation_pearson <- data.correlation_pearson[,-1]

data.correlation_spearman <- fread("phenotype_correlation_spearman.tsv", header = T, sep = "\t", data.table = F)
rownames(data.correlation_spearman) <- data.correlation_spearman[,1]
data.correlation_spearman <- data.correlation_spearman[,-1]

# data.lm_r2 <- fread("phenotype_regression_R2.tsv", header = T, sep = "\t", data.table = F)
# rownames(data.lm_r2) <- data.lm_r2[,1]
# data.lm_r2 <- data.lm_r2[,-1]
# 
# data.lm_p <- fread("phenotype_regression_pvalue.tsv", header = T, sep = "\t", data.table = F)
# rownames(data.lm_p) <- data.lm_p[,1]
# data.lm_p <- data.lm_p[,-1]

## Order samples by study
sample_studies <- limma::strsplit2(colnames(data.correlation_spearman), "_")[,1]
studies_ordered <- unique(sample_studies) # It is ordered by the previous script
samples_ordered <- colnames(data.correlation_spearman) # It is ordered by the previous script
xpos <- cumsum(table(sample_studies)[studies_ordered])
data.annot_studies <- data.frame(Study = sample_studies)
rownames(data.annot_studies) <- samples_ordered
fwrite(as.list(studies_ordered), "legend.heatmap.txt", sep = "\n")
curated.studies <- subset(data.annot_studies, Study %in% paste0("S", 1:41)) 

## Draw Spearman's correlation heatmap, all, no thresholding
data.heatmap <- data.correlation_spearman
data.heatmap[lower.tri(data.heatmap)] <- NA
png("correlation_heatmap_spearman.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, gaps_col = xpos, gaps_row=xpos, annotation_col = data.annot_studies, annotation_names_col = F, annotation_row = data.annot_studies, annotation_names_row = F, color = colorRampPalette(c("blue", "white", "red"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
dev.off()

pdf("correlation_heatmap_spearman.pdf", width = 20, height = 20)
pheatmap(data.heatmap, legend = F, annotation_legend = F, gaps_col = xpos, gaps_row=xpos, annotation_col = data.annot_studies, annotation_names_col = F, annotation_row = data.annot_studies, annotation_names_row = F, color = colorRampPalette(c("blue", "white", "red"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
dev.off()

## Draw Spearman's correlation heatmap, all, no thresholding (curated studies)
png("correlation_heatmap_spearman_curated.png", width = 2000, height = 2000)
pheatmap(data.heatmap[row.names(curated.studies), row.names(curated.studies)], legend = F, annotation_legend = F, gaps_col = xpos[unique(curated.studies$Study)], gaps_row=xpos[unique(curated.studies$Study)], annotation_col = data.annot_studies, annotation_names_col = F, annotation_row = data.annot_studies, annotation_names_row = F, color = colorRampPalette(c("blue", "white", "red"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
dev.off()

pdf("correlation_heatmap_spearman_curated.pdf", width = 20, height = 20)
pheatmap(data.heatmap[row.names(curated.studies), row.names(curated.studies)], legend = T, annotation_legend = F, gaps_col = xpos[unique(curated.studies$Study)], gaps_row=xpos[unique(curated.studies$Study)], annotation_col = data.annot_studies, annotation_names_col = F, annotation_row = data.annot_studies, annotation_names_row = F, color = colorRampPalette(c("blue", "white", "red"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
dev.off()

## Draw Spearman's correlation heatmap, thresholded (all studies)
data.heatmap[abs(data.heatmap) < 0.3] <- 0
data.heatmap[abs(data.heatmap) >= 0.3] <- 1

png("correlation_heatmap_spearman_thresholded.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, gaps_col = xpos, gaps_row=xpos, annotation_col = data.annot_studies, annotation_names_col = F, annotation_row = data.annot_studies, annotation_names_row = F, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
dev.off()

pdf("correlation_heatmap_spearman_thresholded.pdf", width = 20, height = 20)
pheatmap(data.heatmap, legend = F, annotation_legend = F, gaps_col = xpos, gaps_row=xpos, annotation_col = data.annot_studies, annotation_names_col = F, annotation_row = data.annot_studies, annotation_names_row = F, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
dev.off()

## Draw Spearman's correlation heatmap, thresholded (curated studies)
png("correlation_heatmap_spearman_thresholded_curated.png", width = 2000, height = 2000)
pheatmap(data.heatmap[row.names(curated.studies), row.names(curated.studies)], legend = F, annotation_legend = F, gaps_col = xpos[unique(curated.studies$Study)], gaps_row=xpos[unique(curated.studies$Study)], annotation_col = data.annot_studies, annotation_names_col = F, annotation_row = data.annot_studies, annotation_names_row = F, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
dev.off()

pdf("correlation_heatmap_spearman_thresholded_curated.pdf", width = 20, height = 20)
pheatmap(data.heatmap[row.names(curated.studies), row.names(curated.studies)], legend = F, annotation_legend = F, gaps_col = xpos[unique(curated.studies$Study)], gaps_row=xpos[unique(curated.studies$Study)], annotation_col = data.annot_studies, annotation_names_col = F, annotation_row = data.annot_studies, annotation_names_row = F, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
dev.off()
