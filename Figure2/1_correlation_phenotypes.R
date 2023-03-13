##################################################
## Project: DGRPool
## Script purpose: Calculating phenotype-phenotype correlations,R2,p-value from RDS
## Version: 1.0.0
## Date Created: 2022 Dec 12
## Date Modified: 2022 Dec 22
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Working directory
setwd("DGRPool/") # Root folder

# Libraries
suppressPackageStartupMessages(library(data.table))

## Load phenotype data (run download_phenotypes.R script first)
data.all_pheno <- readRDS(file = "RDS/data.all_pheno_10_03_03.rds")

## Correlations
allpheno <- c(paste0(colnames(data.all_pheno[["F"]]), "_F")[2:ncol(data.all_pheno[["F"]])], paste0(colnames(data.all_pheno[["M"]]), "_M")[2:ncol(data.all_pheno[["M"]])], paste0(colnames(data.all_pheno[["NA"]]), "_NA")[2:ncol(data.all_pheno[["NA"]])])
# Create empty correlation matrix
data.correlation_pearson <- data.frame(matrix(nrow = length(allpheno), ncol = length(allpheno)))
colnames(data.correlation_pearson) <- allpheno
rownames(data.correlation_pearson) <- allpheno
# Create empty correlation matrix
data.correlation_spearman<- data.frame(matrix(nrow = length(allpheno), ncol = length(allpheno)))
colnames(data.correlation_spearman) <- allpheno
rownames(data.correlation_spearman) <- allpheno
# Create empty R2 matrix
data.lm_r2<- data.frame(matrix(nrow = length(allpheno), ncol = length(allpheno)))
colnames(data.lm_r2) <- allpheno
rownames(data.lm_r2) <- allpheno
# Create empty pvalue matrix
data.lm_p<- data.frame(matrix(nrow = length(allpheno), ncol = length(allpheno)))
colnames(data.lm_p) <- allpheno
rownames(data.lm_p) <- allpheno
# Progress bar
pb <- txtProgressBar(min = 0, max = length(allpheno) * length(allpheno), style = 3, width = 50, char = "=")
cnt <- 0
for(sex_1 in c("F", "M", "NA")) {
	phenotype_all_1 <- data.all_pheno[[sex_1]]
	for(i in 2:ncol(phenotype_all_1)){ # First column is DGRP
		phenotype_1 <- phenotype_all_1[,i]
		phenotype_1_name <- paste0(colnames(phenotype_all_1)[i], "_", sex_1)
		non_na_1 <- which(!is.na(phenotype_1))
		dgrp_name_1 <- rownames(phenotype_all_1)[non_na_1]
		
		for(sex_2 in c("F", "M", "NA")) {
			phenotype_all_2 <- data.all_pheno[[sex_2]]
			for(j in 2:ncol(phenotype_all_2)){ # First column is DGRP
				phenotype_2 <- phenotype_all_2[,j]
				phenotype_2_name <- paste0(colnames(phenotype_all_2)[j], "_", sex_2)
				non_na_2 <- which(!is.na(phenotype_2))
				dgrp_name_2 <- rownames(phenotype_all_2)[non_na_2]
				dgrp_name <- intersect(dgrp_name_1, dgrp_name_2)
				
				# Compute correlation
				if(length(dgrp_name) > 1)
				{
					phenotype_1 <- phenotype_all_1[dgrp_name,i]
					phenotype_2 <- phenotype_all_2[dgrp_name,j]
					if(var(phenotype_1) != 0 & var(phenotype_2) != 0)
					{
						c_pearson <- cor(phenotype_1, phenotype_2, method = "pearson")
						c_spearman <- cor(phenotype_1, phenotype_2, method = "spearman")
						suppressWarnings(data_regression <- summary(lm(phenotype_2 ~ phenotype_1)))
						r2 <- data_regression$r.squared
						if(nrow(data_regression$coefficients) == 2){ # If only 2 values, regression fails
							pval <- data_regression$coefficients[2,4]
						} else {
							pval <- NA
						}
					} else {
						c_pearson <- NA
						c_spearman <- NA
						r2 <- NA
						pval <- NA
					}
				} else {
					c_pearson <- NA
					c_spearman <- NA
					r2 <- NA
					pval <- NA
				}
				data.correlation_pearson[phenotype_1_name, phenotype_2_name] <- c_pearson
				data.correlation_spearman[phenotype_1_name, phenotype_2_name] <- c_spearman
				data.lm_r2[phenotype_1_name, phenotype_2_name] <- r2
				data.lm_p[phenotype_1_name, phenotype_2_name] <- pval
				cnt <- cnt + 1
				setTxtProgressBar(pb, cnt)
			}
		}
	}
}

## Reorganize by study
reorg_cols <- c()
reorg_sex <- c()
for(s in unique(limma::strsplit2(colnames(data.correlation_pearson), split = "_")[,1])){
	data.this_study <- colnames(data.correlation_pearson)[startsWith(colnames(data.correlation_pearson), prefix = paste0(s,"_"))]
	if(length(data.this_study) > 0){
	reorg_cols <- c(reorg_cols, data.this_study)
	se <- c()
	for(row in data.this_study){
		r <- strsplit(row, split = "_")[[1]]
		se <- c(se, r[length(r)])
	}
	reorg_sex <- c(reorg_sex, se)
	}
}
data.correlation_pearson <- data.correlation_pearson[reorg_cols, reorg_cols]
fwrite(data.correlation_pearson, file = "phenotype_correlation_pearson.tsv", row.names = T, col.names = T, quote = F, sep = "\t")
data.correlation_spearman <- data.correlation_spearman[reorg_cols, reorg_cols]
fwrite(data.correlation_spearman, file = "phenotype_correlation_spearman.tsv", row.names = T, col.names = T, quote = F, sep = "\t")
data.lm_r2 <- data.lm_r2[reorg_cols, reorg_cols]
fwrite(data.lm_r2, file = "phenotype_regression_R2.tsv", row.names = T, col.names = T, quote = F, sep = "\t")
data.lm_p <- data.lm_p[reorg_cols, reorg_cols]
fwrite(data.lm_p, file = "phenotype_regression_pvalue.tsv", row.names = T, col.names = T, quote = F, sep = "\t")
