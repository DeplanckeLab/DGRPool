##################################################
## Project: DGRPool
## Script purpose: Plotting phenotype-phenotype correlation boxplot (within- and cross-studies)
## Version: 1.0.0
## Date Created: 2023 Mar 20
## Date Modified: 2024 Jul 22
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))

# Load correlation data
data.correlation_spearman <- fread("../resources/phenotype_correlation_spearman.tsv", header = T, sep = "\t", data.table = F)
rownames(data.correlation_spearman) <- data.correlation_spearman[,1]
data.correlation_spearman <- data.correlation_spearman[,-1]

# Order samples by study
sample_studies <- limma::strsplit2(colnames(data.correlation_spearman), "_")[,1]
studies_ordered <- unique(sample_studies) # It is ordered by the previous script
samples_ordered <- colnames(data.correlation_spearman) # It is ordered by the previous script
data.annot_studies <- data.frame(Study = sample_studies)
rownames(data.annot_studies) <- samples_ordered

## Find curated studies
phenotypes_information <- fread("../resources/phenotypes.tsv", data.table = F)
phenotypes_information <- phenotypes_information[,c("study_id", "study_name", "study_status")]
phenotypes_information <- phenotypes_information[!duplicated(phenotypes_information),]
phenotypes_information$study_name <- paste0("(",phenotypes_information$study_name,")")
rownames(phenotypes_information) <- paste0("S", phenotypes_information$study_id)
curated_studies <- paste0("S", sort(phenotypes_information$study_id[phenotypes_information$study_status == "integrated"]))
phenotypes_information <- phenotypes_information[curated_studies,]

# Values for within- and cross-study correlations
results <- data.frame(matrix(nrow = (nrow(data.correlation_spearman)*nrow(data.correlation_spearman) - nrow(data.correlation_spearman))/2, ncol = 2))
colnames(results) <- c("Value", "Type")
results$Curated <- F
n <- 1
pb <- txtProgressBar(min = 1, max = nrow(results), style = 3, width = 50, char = "=")
for(i in 1:(nrow(data.correlation_spearman) - 1)){
	for(j in (i+1):nrow(data.correlation_spearman)){
		results[n, "Value"] <- data.correlation_spearman[i, j]
		if(sample_studies[i] == sample_studies[j]){
			results[n, "Type"] <- "Within-study"
			if(sample_studies[i] %in% curated_studies) results[n, "Curated"] <- T
		} else {
			results[n, "Type"] <- "Cross-study"
			if(sample_studies[i] %in% curated_studies) {
				if(sample_studies[j] %in% curated_studies) results[n, "Curated"] <- T
			}
		}
		n <- n + 1
		setTxtProgressBar(pb, n)
	}
}

# Clean results
results$Value <- abs(as.numeric(results$Value))
results <- subset(results, !is.na(Value))

# Barplot
data.sum <- summarise(group_by(results, Type), MD = median(Value), MN = mean(Value))
p <- ggplot(results, aes(x = Type, y = Value, fill = Type)) + geom_boxplot(show.legend = FALSE, width=0.5) + theme_classic() + scale_fill_manual(values = c("#8ba8c8", "#E69A8DFF"))
p <- p + stat_compare_means(method = "wilcox.test", paired = F, comparisons = list( c("Cross-study", "Within-study")))
p <- p + geom_text(data = data.sum, aes(Type, MD, label = round(MD,3)), position = position_dodge(width = 0.8), size = 3, vjust = -0.5, color = "black")
p <- p + xlab("") + ylab("Spearman's correlation")
ggsave(p, filename = "boxplot.within.cross.study.pdf", width = 3, height = 4)
#ggsave(p, filename = "boxplot.within.cross.study.png", width = 3, height = 4)
message("Displayed values are median. Mean values are ", round(subset(data.sum, Type == "Cross-study")$MN, 3) , " for cross-study, and ", round(subset(data.sum, Type == "Within-study")$MN, 3), " for Within-study")
