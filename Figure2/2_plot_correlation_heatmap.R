##################################################
## Project: DGRPool
## Script purpose: Plotting phenotype-phenotype correlation heatmap
## Version: 1.0.0
## Date Created: 2022 Dec 12
## Date Modified: 2023 Mar 17
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
data.correlation_pearson <- fread("phenotype_correlation_pearson.tsv", header = T, sep = "\t", data.table = F)
rownames(data.correlation_pearson) <- data.correlation_pearson[,1]
data.correlation_pearson <- data.correlation_pearson[,-1]

data.correlation_spearman <- fread("phenotype_correlation_spearman.tsv", header = T, sep = "\t", data.table = F)
rownames(data.correlation_spearman) <- data.correlation_spearman[,1]
data.correlation_spearman <- data.correlation_spearman[,-1]

data.lm_r2 <- fread("phenotype_regression_R2.tsv", header = T, sep = "\t", data.table = F)
rownames(data.lm_r2) <- data.lm_r2[,1]
data.lm_r2 <- data.lm_r2[,-1]

data.lm_p <- fread("phenotype_regression_pvalue.tsv", header = T, sep = "\t", data.table = F)
rownames(data.lm_p) <- data.lm_p[,1]
data.lm_p <- data.lm_p[,-1]

## Remove rows / columns of NAs
to_keep <- rownames(data.correlation_pearson)[apply(data.correlation_pearson, 1, function(x) sum(is.na(x))) != nrow(data.correlation_pearson)]
data.correlation_pearson <- data.correlation_pearson[to_keep,to_keep]
data.correlation_spearman <- data.correlation_spearman[to_keep,to_keep]
data.lm_r2 <- data.lm_r2[to_keep,to_keep]
data.lm_p <- data.lm_p[to_keep,to_keep]

## Draw Heatmap
data.heatmap <- data.correlation_spearman
data.heatmap[lower.tri(data.heatmap)] <- NA

data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.heatmap), split = "_")[,1])
rownames(data.annot_studies) <- colnames(data.heatmap)

png("correlation_heatmap_spearman.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("blue", "white", "red"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()


## Draw Heatmap
data.heatmap <- data.correlation_pearson
data.heatmap[lower.tri(data.heatmap)] <- NA

data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.heatmap), split = "_")[,1])
rownames(data.annot_studies) <- colnames(data.heatmap)

png("correlation_heatmap_pearson.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("blue", "white", "red"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()

pdf("correlation_heatmap_pearson.pdf", width = 20, height = 20)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("blue", "white", "red"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()

## Draw Heatmap
data.heatmap <- data.lm_r2
data.heatmap[lower.tri(data.heatmap)] <- NA

data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.heatmap), split = "_")[,1])
rownames(data.annot_studies) <- colnames(data.heatmap)

png("regression_heatmap_R2.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("blue", "white", "red"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()

## Draw Heatmap
data.heatmap <- data.lm_p
data.heatmap[lower.tri(data.heatmap)] <- NA

data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.heatmap), split = "_")[,1])
rownames(data.annot_studies) <- colnames(data.heatmap)

png("regression_heatmap_pvalue.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("blue", "white", "red"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()

## Draw Heatmap
data.heatmap <- data.correlation_spearman
data.heatmap[lower.tri(data.heatmap)] <- NA

data.heatmap[abs(data.heatmap) < 0.3] <- 0
data.heatmap[abs(data.heatmap) >= 0.3] <- 1

data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.heatmap), split = "_")[,1])
rownames(data.annot_studies) <- colnames(data.heatmap)

png("correlation_heatmap_spearman_thresholded.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()

## Draw Heatmap
data.heatmap <- data.correlation_spearman
data.heatmap[lower.tri(data.heatmap)] <- NA

data.heatmap[abs(data.heatmap) < 0.3] <- 0
data.heatmap[abs(data.heatmap) >= 0.3] <- 1

data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.heatmap), split = "_")[,1])
rownames(data.annot_studies) <- colnames(data.heatmap)
data.annot_studies <- subset(data.annot_studies, Study %in% paste0("S", 1:41)) 
data.heatmap <- data.heatmap[rownames(data.annot_studies),rownames(data.annot_studies)]

png("correlation_heatmap_spearman_thresholded_41first.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()

pdf("correlation_heatmap_spearman_thresholded_41first.pdf", width = 20, height = 20)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()

## Draw Heatmap
data.heatmap <- data.correlation_pearson
data.heatmap[lower.tri(data.heatmap)] <- NA

data.heatmap[abs(data.heatmap) < 0.3] <- 0
data.heatmap[abs(data.heatmap) >= 0.3] <- 1

data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.heatmap), split = "_")[,1])
rownames(data.annot_studies) <- colnames(data.heatmap)
data.annot_studies <- subset(data.annot_studies, Study %in% paste0("S", 1:41)) 
data.heatmap <- data.heatmap[rownames(data.annot_studies),rownames(data.annot_studies)]

png("correlation_heatmap_pearson_thresholded_41first.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()

## Draw Heatmap
data.heatmap <- data.lm_p
data.heatmap[lower.tri(data.heatmap)] <- NA
data.heatmap[data.heatmap > 0.05] <- "NS"
data.heatmap[data.heatmap <= 0.05] <- "S"
data.heatmap[data.heatmap == "NS"] <- 0
data.heatmap[data.heatmap == "S"] <- 1
data.heatmap <- data.frame(apply(data.heatmap, 2, as.numeric))
rownames(data.heatmap) <- colnames(data.heatmap)

data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.heatmap), split = "_")[,1])
rownames(data.annot_studies) <- colnames(data.heatmap)

png("regression_heatmap_pvalue_0.05.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()

## Draw Heatmap
data.heatmap <- data.lm_p
data.heatmap[lower.tri(data.heatmap)] <- NA
data.heatmap[data.heatmap > 0.05] <- "NS"
data.heatmap[data.heatmap <= 0.05] <- "S"
data.heatmap[data.heatmap == "NS"] <- 0
data.heatmap[data.heatmap == "S"] <- 1

data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.heatmap), split = "_")[,1])
rownames(data.annot_studies) <- colnames(data.heatmap)
data.annot_studies <- subset(data.annot_studies, Study %in% paste0("S", 1:41)) 
data.heatmap <- data.heatmap[rownames(data.annot_studies),rownames(data.annot_studies)]
data.heatmap <- data.frame(apply(data.heatmap, 2, as.numeric))
rownames(data.heatmap) <- colnames(data.heatmap)

png("regression_heatmap_pvalue_0.05_41first.png", width = 2000, height = 2000)
pheatmap(data.heatmap, legend = F, annotation_legend = F, annotation_names_col = F, annotation_col = data.annot_studies, color = colorRampPalette(c("white", "black"))(100), show_rownames = F, show_colnames = F, cluster_rows=F, cluster_cols=F, na_col="white", border_color ="white")
t <- table(data.annot_studies$Study)[unique(data.annot_studies$Study)]
xpos <- c()
cumsumx <- 0
for(i in 1:length(t)){
	xpos <- c(xpos, t[i] / 2 + cumsumx)
	cumsumx <- cumsumx + t[i]
}
xpos <- (xpos / cumsumx) * 0.995 + 0.002 # 0.9 is scale factor
grid.text(unique(data.annot_studies$Study), x=xpos, y=0.995, gp=gpar(fontsize=6), rot = 90)
dev.off()