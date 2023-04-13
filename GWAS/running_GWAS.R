##################################################
## Project: DGRPool
## Script purpose: Running GWAS
## Version: 1.0.0
## Date Created: 2022 Dec 22
## Date Modified: 2023 Mar 10
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

## Working directory
setwd("/data/gardeux/DGRPool/")
dir.create("GWAS", mode="0750", showWarnings = F)

## Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(require(ramwas)) # QQ-Plot & Manhattan plot
suppressPackageStartupMessages(require(car)) # Covariate analysis

## Functions
toReadableTime = function(s){
	s <- round(s)
	if(s < 60) return(paste0(s, "s"))
	mn <- s %/% 60
	s <- s %% 60
	if(mn < 60) return(paste0(mn, "mn ", s, "s"))
	h <- mn %/% 60
	mn <- mn %% 60
	if(h < 24) return(paste0(h, "h ",  mn, "mn ", s, "s"))
	d <- h %/% 24
	h <- h %% 24
	return(paste0(d, "d ", h, "h ", mn, "mn ", s, "s"))
}

## Load phenotype data (run download_phenotypes.R script first)
data.all_pheno <- readRDS(file = "data.all_pheno_21_03_23_filtered.rds")
# DGRP names to change
for(s in names(data.all_pheno)){
	rownames(data.all_pheno[[s]]) <- gsub(rownames(data.all_pheno[[s]]), pattern = "DGRP_", replacement = "line_")
	rownames(data.all_pheno[[s]]) <- gsub(rownames(data.all_pheno[[s]]), pattern = "line_0", replacement = "line_")
	rownames(data.all_pheno[[s]]) <- gsub(rownames(data.all_pheno[[s]]), pattern = "line_0", replacement = "line_")
}
allpheno <- c(paste0(colnames(data.all_pheno[["F"]]), "_F")[2:ncol(data.all_pheno[["F"]])], paste0(colnames(data.all_pheno[["M"]]), "_M")[2:ncol(data.all_pheno[["M"]])], paste0(colnames(data.all_pheno[["NA"]]), "_NA")[2:ncol(data.all_pheno[["NA"]])])
message(length(allpheno), " phenotypes to process")

## Load covariates
# For Wolbachia: y or n
# For big inversions: ST (0), INV (1) or INV/ST (2)
dgrp.cov <- fread ("dgrp.cov.tsv", header = T, sep = "\t")
rownames(dgrp.cov) <- dgrp.cov$IID
dgrp_lines_genotyped <- dgrp.cov$IID
message("Found ", length(dgrp_lines_genotyped), " genotyped DGRP lines")
message(ncol(dgrp.cov) - 2, " covariates to account for : [", paste0(colnames(dgrp.cov)[3:ncol(dgrp.cov)], collapse = ", "), "]")

## Load VCF
# data.vcf <- fread("dgrp2.vcf", data.table = F) # 4438427 variants
# data.vcf <- subset(data.vcf, ID %in% data.vcf2$ID) # Filtering on variants after running web tool dgrp2
# # Write filtered VCF file
# fileConn<-file("dgrp2.filtered.vcf")
# writeLines(c("##fileformat=VCFv4.1","##source=DGRP2","##reference=dm3", "##INFO=<ID=REFCOUNT,Number=1,Type=Integer,Description=\"Reference Allele Count\">","##INFO=<ID=ALTCOUNT,Number=1,Type=Integer,Description=\"Alternative Allele Count\">","##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"), fileConn)
# close(fileConn)
# fwrite(x = data.vcf, file = "dgrp2.filtered.vcf", sep = "\t", quote = F, row.names = F, col.names = T, append = T)
# # Now we need to change the chromosome names for the PLINK file
# data.vcf$`#CHROM`[data.vcf$`#CHROM` == "2L"] = 1
# data.vcf$`#CHROM`[data.vcf$`#CHROM` == "2R"] = 2
# data.vcf$`#CHROM`[data.vcf$`#CHROM` == "3L"] = 3
# data.vcf$`#CHROM`[data.vcf$`#CHROM` == "4"] = 6
# data.vcf$`#CHROM`[data.vcf$`#CHROM` == "3R"] = 4
# data.vcf$`#CHROM`[data.vcf$`#CHROM` == "X"] = 5
# # Write filtered VCF file
# fileConn<-file("dgrp2.filtered.fastlmm.vcf")
# writeLines(c("##fileformat=VCFv4.1","##source=DGRP2","##reference=dm3", "##INFO=<ID=REFCOUNT,Number=1,Type=Integer,Description=\"Reference Allele Count\">","##INFO=<ID=ALTCOUNT,Number=1,Type=Integer,Description=\"Alternative Allele Count\">","##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"), fileConn)
# close(fileConn)
# fwrite(x = data.vcf, file = "dgrp2.filtered.fastlmm.vcf", sep = "\t", quote = F, row.names = F, col.names = T, append = T)
# # Create PLINK files from VCF
# system("plink2 --vcf dgrp2.filtered.fastlmm.vcf --make-bed --out dgrp2.filtered.fastlmm")
# # Replace the .fam file which is wrongly generated!
# file.copy("dgrp2.fam", "dgrp2.filtered.fastlmm.fam")
#data.vcf <- fread("dgrp2.filtered.vcf", data.table = F) # 1882162 variants (with 2L, 2R, 3L, 3R, X, 4 as chromosome names, not 1, 2, 3, 4, 5, 6)

## Load annotation
annotation.dr <- fread ("/data/genome/drosophila_melanogaster/dm3/dgrp.fb557.annot.txt", header = FALSE, sep = "\t")
annotation.dr <- annotation.dr[,c("V1", "V3", "V4")]
colnames(annotation.dr) <- c("SNP", "gene_annotation", "regulatory_annotation")

## Running GWAS with progress bar
start_time = Sys.time()
pb <- txtProgressBar(min = 0, max = length(allpheno), style = 3, width = 50, char = "=")
cnt <- 0
for(sex in c("F", "M", "NA")) {
	phenotype_all <- data.all_pheno[[sex]]
	for(i in 2:ncol(phenotype_all)){ # First column is DGRP
		## Extracting phenotype, removing NAs
		phenotype <- phenotype_all[,i]
		phenotype_name <- paste0(colnames(phenotype_all)[i], "_", sex)
		# Fix weird names
		#phenotype_name <- gsub(phenotype_name, pattern = " ", replacement = "_")
		#phenotype_name <- gsub(phenotype_name, pattern = "\ub0", replacement = "") # degree sign
		#phenotype_name <- gsub(phenotype_name, pattern = "\\(", replacement = "")
		#phenotype_name <- gsub(phenotype_name, pattern = "\\)", replacement = "")
		#phenotype_name <- gsub(phenotype_name, pattern = "-", replacement = "_")
		#phenotype_name <- gsub(phenotype_name, pattern = "/", replacement = "_")
		#phenotype_name <- gsub(phenotype_name, pattern = ":", replacement = "")
		#phenotype_name <- gsub(phenotype_name, pattern = "%", replacement = "")
		non_na <- which(!is.na(phenotype))
		dgrp_name <- rownames(phenotype_all)[non_na]
		dgrp_name <- gsub(dgrp_name, pattern = "line_0", replacement = "line_")
		phenotype <- phenotype_all[non_na,i]
		
		## Generating phenotype file for PLINK
		plink.file_content <- data.frame(FID = dgrp_name, IID = dgrp_name)
		plink.file_content[phenotype_name] <- phenotype
		fwrite(x = format(plink.file_content, nsmall = 1), file = paste0("GWAS_full_results/", phenotype_name, ".pheno.all.tsv"), sep = "\t", quote = F, row.names = F, col.names = T, scipen=50)
		
		## Restrict to genotyped samples
		plink.file_content <- subset(plink.file_content, IID %in% dgrp_lines_genotyped)
		fwrite(x = format(plink.file_content, nsmall = 1), file = paste0("GWAS/", phenotype_name, ".pheno.plink2.tsv"), sep = "\t", quote = F, row.names = F, col.names = T, scipen=50)
		## Create a version for the DGRP2 website
		dgrp2.file_content <- plink.file_content[c(2,3)]
		dgrp2.file_content$IID <- gsub(x = dgrp2.file_content$IID, pattern = "line_", replacement = "")
		fwrite(x = format(dgrp2.file_content, nsmall = 1), file = paste0("GWAS/", phenotype_name, ".pheno.dgrp2.csv"), sep = ",", quote = F, row.names = F, col.names = F, scipen=50)
		
		## Test for categorical phenotype
		is_categorical <- F
		tryCatch({
			as.numeric(plink.file_content[[phenotype_name]])
		}, warning=function(cond) {
				if(grepl(x = cond, pattern = "NAs introduced by coercion")) {
					message("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: CATEGORICAL PHENOTYPE?")
					is_categorical <<- T
				}
		})
		
		if(!is_categorical){
			## Running GWAS
			list_values <- unique(plink.file_content[,3])
			if(length(list_values) == 0){
				message("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: No values.")
			} else if(length(list_values) == 1){
				message("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: Non-variable phenotype.")
			} else if(length(list_values) == 2 && all(sort(list_values) == c(0, 2))){
				message("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: WEIRD phenotype? Only values in [0, 2]") # Error: All samples for --glm phenotype 'S7_MB_32_M' are cases.
			}	else if(length(list_values) == 2 && all(sort(list_values) == c(0, 1))){
				message("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: WEIRD phenotype? Only values in [0, 1]") # Error: All samples for --glm phenotype 'S7_MB_32_M' are cases.
			} else {
				## Covariate analysis (Anova)
				sink(paste0("GWAS/", phenotype_name, ".cov.anova.txt"))
				data.cov <- merge(plink.file_content, dgrp.cov, id = "family")
				## Testing variance of covariates
				# Generate formula
				covars <- c()
				if(length(unique(data.cov[,4])) > 1) covars <- c(covars, "factor(wolba)")
				if(length(unique(data.cov[,5])) > 1) covars <- c(covars, "factor(In_2L_t)")
				if(length(unique(data.cov[,6])) > 1) covars <- c(covars, "factor(In_2R_NS)")
				if(length(unique(data.cov[,7])) > 1) covars <- c(covars, "factor(In_3R_P)")
				if(length(unique(data.cov[,8])) > 1) covars <- c(covars, "factor(In_3R_K)")
				if(length(unique(data.cov[,9])) > 1) covars <- c(covars, "factor(In_3R_Mo)")
				if(length(covars) > 0){
					fit <- lm(data = data.cov, formula = as.formula(paste0(phenotype_name, " ~ ", paste0(covars, collapse = " + "))))
					if(deviance(fit) >= sqrt(.Machine$double.eps)){
						if(any(is.na(fit$coefficients))){
							message("ERROR: In Anova.III.lm(mod, error, singular.ok = singular.ok, ...) : there are aliased coefficients in the model")
						} else {
							print(Anova(fit, type = 'III'))
						}
					}
					print(summary(fit))
				} else message("ERROR: No invariant covariates??")
				sink() # returns output to the console
				
				if(length(covars) + 2 >= nrow(plink.file_content)){
					message("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: # samples <= # predictor columns.")
				} else {
					# PLINK 1.x (Same results, but slower)
					#std.out <- system(paste0("plink --threads 64 --linear --covar /data/gardeux/DGRPool/dgrp.cov.fastlmm.tsv --bfile /data/gardeux/DGRPool/dgrp2.filtered.fastlmm --pheno GWAS/", phenotype_name, ".pheno.plink2.tsv --out  GWAS/plink1"), intern = T)
					#file_name <- std.out[sapply(std.out, function(x) startsWith(x = x, prefix = "Writing linear model association results to"))]
					#file_name <- gsub(file_name, pattern = ".*GWAS/", replacement = "GWAS/")
					#file_name <- gsub(file_name, pattern = " ... 0%.*", replacement = "")
					# PLINK 2.x 
					std.out <- NULL
					tryCatch(
					{
						std.out <<- system(paste0("plink2 --threads 64 --glm hide-covar --geno 0.2 --maf 0.05 --covar /data/gardeux/DGRPool/dgrp.cov.tsv --bfile /data/gardeux/DGRPool/dgrp2 --pheno GWAS/", phenotype_name, ".pheno.plink2.tsv --out  GWAS_full_results/"), intern = T, ignore.stderr = T, ignore.stdout = F)
					},
					warning=function(cond) {
						if(grepl(x = cond, pattern = "had status 7")) message("ERROR [HAD STATUS 7] in processing phenotype ", phenotype_name, " [i = ", i, "]: # samples <= # predictor columns.")
						else std.out <<- system(paste0("plink2 --threads 64 --glm hide-covar --geno 0.2 --maf 0.05 --covar /data/gardeux/DGRPool/dgrp.cov.tsv --bfile /data/gardeux/DGRPool/dgrp2 --pheno GWAS/", phenotype_name, ".pheno.plink2.tsv --out  GWAS_full_results/"), intern = T, ignore.stderr = T, ignore.stdout = F)
					})
					
					if(!is.null(std.out)) {
					
						file_name <- std.out[sapply(std.out, function(x) startsWith(x = x, prefix = "Results written to"))]
						# Rename output
						if(length(file_name) > 0)
						{
							file_name <- gsub(x = file_name, pattern = "Results written to ", replacement = "")
							file_name <- gsub(x = file_name, pattern = " .", replacement = "")
							std.out <- file.rename(file_name, gsub(x = file_name, pattern = "/.", replacement = "/"))
							file_name <- gsub(x = file_name, pattern = "/.", replacement = "/")
						}
			
						## Reading Output
						plink_results <- NA
						col6 <- "BETA"
						file_suffix <- ".glm.linear"
						if(file.exists(paste0("GWAS_full_results/", phenotype_name, ".glm.linear"))){
							# Numeric / continuous
							plink_results <<- fread(paste0("GWAS_full_results/", phenotype_name, ".glm.linear"), data.table = F, nThread = 64, showProgress = F)
						} else if(file.exists(paste0("GWAS_full_results/", phenotype_name, ".glm.logistic.hybrid"))){
							# Categorical
							plink_results <<- fread(paste0("GWAS_full_results/", phenotype_name, ".glm.logistic.hybrid"), data.table = F, nThread = 64, showProgress = F)
							col6 <<- "OR"
							file_suffix <<- ".glm.logistic.hybrid"
						} else stop("ERROR? What is the output file?")
						plink_results$P <- as.numeric(plink_results$P)
						plink_results <- subset(plink_results, TEST == "ADD")
						plink_results <- plink_results[!is.na(plink_results$P),c("#CHROM", "POS", "ID", "REF", "ALT", col6, "P")]
						plink_results$FDR_BH <- p.adjust(plink_results$P, method = "fdr")
						
						## QQPLOT before filtering
						plink_results$P[plink_results$P == 0] <- 1E-255 # For avoiding errors
						png(paste0("GWAS/", phenotype_name, ".qqplot.png"), height = 500, width = 500)
						#gap::qqunif(plink_results$P) # Slow
						#qq(plink_results$P) # Slower
						qqPlotFast(plink_results$P, ci.level = NULL, col = "black", makelegend = F, lwd = 1, newplot = T)
						dev.off()
						
						# SVG is way too big. So let's go PDF
						# I recompute the plot... I know it's not optimal but I don't know how to do it else
						pdf(paste0("GWAS/", phenotype_name, ".qqplot.pdf"), height = 7, width = 7)
						#gap::qqunif(plink_results$P) # Slow
						#qq(plink_results$P) # Slower
						qqPlotFast(plink_results$P, ci.level = NULL, col = "black", makelegend = F, lwd = 1, newplot = T)
						dev.off()
			
						## Manhattan plot
						chroms <- factor(plink_results$`#CHROM`)
						levels(chroms) <- c("2L", "2R", "3L", "3R", "X", "4")
						m <- manPlotPrepare(pvalues = plink_results$P, chr=chroms, pos = plink_results$POS)
						png(paste0("GWAS/", phenotype_name, ".manhattan.png"), height = 500, width = 500)
						manPlotFast(man = m, lwd = 1, colorSet = c("black", "darkgrey"))
						#manhattan(plink_results, chr="#CHROM", bp = "POS", snp = "ID", chrlabs = c("2L", "2R", "3L", "3R", "X", "4"), suggestiveline = F, genomewideline = F, annotateTop = T, annotatePval = T)
						dev.off()
						
						# SVG is way too big. So let's go PDF
						# I recompute the plot... I know it's not optimal but I don't know how to do it else
						pdf(paste0("GWAS/", phenotype_name, ".manhattan.pdf"), height = 7, width = 7)
						manPlotFast(man = m, lwd = 1, colorSet = c("black", "darkgrey"))
						#manhattan(plink_results, chr="#CHROM", bp = "POS", snp = "ID", chrlabs = c("2L", "2R", "3L", "3R", "X", "4"), suggestiveline = F, genomewideline = F, annotateTop = T, annotatePval = T)
						dev.off()
						
						## Filtered top results and annotate them
						plink_results <- subset(plink_results, P <= 0.01)
						plink_results$P[plink_results$P == 1E-255] <- 0 # Putting it back
						plink_results <- merge(plink_results, annotation.dr, by.x = "ID", by.y = "SNP")
						plink_results <- plink_results[,c("#CHROM", "POS", "ID", "REF", "ALT", col6, "P", "FDR_BH", "gene_annotation", "regulatory_annotation")]
						plink_results <- plink_results[with(plink_results, order(P)),]
						plink_results$`#CHROM`[plink_results$`#CHROM` == 1] <- "2L"
						plink_results$`#CHROM`[plink_results$`#CHROM` == 2] <- "2R"
						plink_results$`#CHROM`[plink_results$`#CHROM` == 3] <- "3L"
						plink_results$`#CHROM`[plink_results$`#CHROM` == 4] <- "3R"
						plink_results$`#CHROM`[plink_results$`#CHROM` == 5] <- "X"
						plink_results$`#CHROM`[plink_results$`#CHROM` == 6] <- "4"
						fwrite(x = plink_results, file = paste0("GWAS/", phenotype_name, file_suffix, ".top_0.01.annot.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
						system(paste0("gzip -f GWAS/", phenotype_name, file_suffix, ".top_0.01.annot.tsv"))
					}
				}
			}
		} # !is_categorical
		## Progression bar
		cnt <- cnt + 1
		setTxtProgressBar(pb, cnt)
	}
}

message("Total time = ", toReadableTime(as.numeric(Sys.time() - start_time, unit = "secs"))) # Total time = 2h 15mn 56s


