##################################################
## Project: DGRPool
## Script purpose: Use correlation matrix generated in first script to plot phenotype-phenotype correlation heatmap
## Version: 1.0.0
## Date Created: 2022 Dec 12
## Date Modified: 2022 Mar 10
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

## Working directory
setwd("DGRPool/") # Root folder

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
data.heatmap <- data.lm_p
data.heatmap[lower.tri(data.heatmap)] <- NA
data.heatmap[data.heatmap > 0.05] <- "NS"
data.heatmap[data.heatmap <= 0.05] <- "S"
data.heatmap[data.heatmap == "NS"] <- 0
data.heatmap[data.heatmap == "S"] <- 1
data.heatmap <- apply(data.heatmap, 2, as.numeric)

data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.heatmap), split = "_")[,1])
rownames(data.annot_studies) <- colnames(data.heatmap)

png("regression_heatmap_pvalue_0.05.png", width = 2000, height = 2000)
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
data.annot_studies <- subset(data.annot_studies, Study %in% paste0("S", 1:40)) 
data.heatmap <- data.heatmap[rownames(data.annot_studies),rownames(data.annot_studies)]

png("correlation_heatmap_spearman_thresholded_40first.png", width = 2000, height = 2000)
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

## Draw Histogram
data.annot_studies <- data.frame(Study = limma::strsplit2(colnames(data.correlation_spearman), split = "_")[,1])
data.annot_studies$Study <- factor(data.annot_studies$Study, levels = paste0("S", 1:122)) 
p <- ggplot(data.annot_studies, aes(x = Study)) + geom_histogram(stat="count") + ylab("Number of phenotypes") + xlab("Study id") + theme(axis.line.y = element_line(size = 0.3, colour = "black"), axis.line.x = element_line(size = 0.3, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, size = 6))
ggsave(plot = p, filename = "NbPhenoPerStudy.plot.pdf", width = 8, height = 4)
