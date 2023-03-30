##################################################
## Project: DGRPool
## Script purpose: Running GWAS for user-submitted data
## Version: 1.0.0
## Date Created: 2023 Mar 29
## Date Modified: 2023 Mar 29
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Measure execution time
start_time = Sys.time()

# Global script options
options(echo = F)
args <- commandArgs(trailingOnly = TRUE)
#args <- c("Study_42.txt", "/data/gardeux/DGRPool/dgrp2", "dgrp.cov.tsv", "dgrp.fb557.annot.txt.gz", "/data/gardeux/DGRPool/Prout", 1, 0.2, 0.05, "Ave. Social Space")
#args <- c("/data/gardeux/DGRPool/GWAS_extreme/extreme_phenotype.tsv", "/data/gardeux/DGRPool/dgrp2", "/data/gardeux/DGRPool/dgrp.cov.tsv", "/data/gardeux/DGRPool/dgrp.fb557.annot.txt.gz", "/data/gardeux/DGRPool/GWAS_extreme/low_extreme/", 32, 0.2, 0.05, "low_extreme")

## Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ramwas)) # QQ-Plot & Manhattan plot
suppressPackageStartupMessages(library(car)) # Covariate analysis

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

stop_quietly <- function(){
	opt <- options(show.error.messages = F)
	on.exit(options(opt))
	q(save="no", status=1, runLast=FALSE)
}

error.json <- function(displayed) {
	stats <- list()
	stats$displayed_error = displayed
	print(toJSON(stats, method="C", auto_unbox=T))
	stop_quietly()
}

# Result file
results <- list()

# Parameters
if(is.null(args) | length(args) < 8){
	message("Usage: running_GWAS_user.R parameters\n1. [Required] Path to phenotype file\n2. [Required] Path to PLINK genotype file prefix (bed)\n3. [Required] Path to covariate file\n4. [Required] Path to annotation file\n5. [Required] Output folder path\n6. [Required] Number of threads\n7. [Required] Minimum genotyping ratio\n8. [Required] Minimum MAF\n9. [Optional] Name of phenotype to load (if file contains multiple phenotypes)")
	error.json("You need to input at least 8 parameters: [Phenotype file path][Genotype file prefix][Covariate file path][Annotation file path][Output folder path][Number of threads][Genotype Ratio][Min MAF]")
}
results$pheno_file <- args[1]
results$genotype_file <- args[2]
results$covariate_file <- args[3]
results$annotation_file <- args[4]
results$output_folder <- args[5]
results$output_folder <- gsub(x = results$output_folder, pattern = "\\\\", replacement = "/")
if(!endsWith(results$output_folder, "/")) results$output_folder <- paste0(results$output_folder, "/")
if(!file.exists(results$output_folder)) dir.create(results$output_folder, mode = "755", showWarnings = FALSE)
results$nb_threads <- as.numeric(args[6])
results$genotype_ratio <- as.numeric(args[7])
results$min_maf <- as.numeric(args[8])
results$pheno_name <- NULL
results$pheno_column <- 3 # Default phenotype to read
results$which_sex <- "NA" # Default, in case no information
if(length(args) == 9){
	# Check if there is a "name" parameter
	results$pheno_name <- args[9]
}
if(length(args) > 9) error.json("You entered too many parameters: ", length(args))

# Read user-submitted phenotyping data
if(!file.exists(results$pheno_file)) error.json(paste0("This file: '",results$pheno_file,"' does not exist!"))
data.dgrp_pheno <- fread(results$pheno_file, sep = "\t", header = T, data.table = F, na.strings = "")
if(colnames(data.dgrp_pheno)[1] != "DGRP") error.json(paste0("First column of ",results$pheno_file," should be the 'DGRP' column"))

# Check sex column
if("sex" %in% colnames(data.dgrp_pheno)){
	# Check if it's correctly in col 2
	where_sex <- which(colnames(data.dgrp_pheno) == "sex")
	if(length(where_sex) > 1) error.json(paste0("Phenotype 'sex' was found multiple times in your phenotyping file. It should be unique, and in column 2."))
	if(where_sex != 2) error.json(paste0("Phenotype 'sex' was found in your phenotyping file, but not in the correct column (", where_sex, "). It should be in column 2."))
	# Check content of column
	data.dgrp_pheno[2][is.na(data.dgrp_pheno[2])] <- "NA"
	results$which_sex <- names(table(data.dgrp_pheno[2]))
	if(any(!results$which_sex %in% c("F", "M", "NA"))) error.json(paste0("Column 'sex' should only contain 'F', 'M', or 'NA' values."))
} else {
	# If no sex column, I create one with only "NA" values
	col_names <- colnames(data.dgrp_pheno)[2:ncol(data.dgrp_pheno)]
	data.dgrp_pheno$sex <- rep("NA", nrow(data.dgrp_pheno))
	data.dgrp_pheno <- data.dgrp_pheno[,c("DGRP", "sex", col_names)]
}
# Check pheno_name, if parameters is provided
if(!is.null(results$pheno_name)){
	# Check if pheno_name is unique and in the file
	results$pheno_column <- which(colnames(data.dgrp_pheno) == results$pheno_name)
	if(length(results$pheno_column) == 0) error.json(paste0("Phenotype '",results$pheno_name,"' was not found in any column of ", results$pheno_file))
	if(length(results$pheno_column) > 1) error.json(paste0("Phenotype '",results$pheno_name,"' was found in many column of ", results$pheno_file))
} else results$pheno_name <- colnames(data.dgrp_pheno)[results$pheno_column]

## Load covariates
# For Wolbachia: y or n
# For big inversions: ST (0), INV (1) or INV/ST (2)
dgrp.cov <- fread (results$covariate_file, header = T, sep = "\t")
rownames(dgrp.cov) <- dgrp.cov$IID
dgrp_lines_genotyped <- dgrp.cov$IID
message("Found ", length(dgrp_lines_genotyped), " genotyped DGRP lines")
message(ncol(dgrp.cov) - 2, " covariates to account for : [", paste0(colnames(dgrp.cov)[3:ncol(dgrp.cov)], collapse = ", "), "]")

## Load annotation
annotation.dr <- fread (results$annotation_file, header = FALSE, sep = "\t")
annotation.dr <- annotation.dr[,c("V1", "V3", "V4")]
colnames(annotation.dr) <- c("SNP", "gene_annotation", "regulatory_annotation")

## Create results
results$results <- list()

## Running GWAS
for(s in results$which_sex) {
	results$results[[s]] <- list()

	## Extracting phenotype
	phenotype <- subset(data.dgrp_pheno, sex == s)[,results$pheno_column]
	phenotype_name <- paste0(results$pheno_name, "_", s)
	dgrp_name <- subset(data.dgrp_pheno, sex == s)[,"DGRP"]
	
	## Transform to numeric
	suppressWarnings(phenotype <<- as.numeric(phenotype))
	
	## Fix weird names
	phenotype_name <- gsub(phenotype_name, pattern = " ", replacement = "_")
	phenotype_name <- gsub(phenotype_name, pattern = "\ub0", replacement = "") # degree sign
	phenotype_name <- gsub(phenotype_name, pattern = "\\(", replacement = "")
	phenotype_name <- gsub(phenotype_name, pattern = "\\)", replacement = "")
	phenotype_name <- gsub(phenotype_name, pattern = "-", replacement = "_")
	phenotype_name <- gsub(phenotype_name, pattern = "/", replacement = "_")
	phenotype_name <- gsub(phenotype_name, pattern = ":", replacement = "")
	phenotype_name <- gsub(phenotype_name, pattern = "%", replacement = "")
	
	## Removing NAs
	non_na <- which(!is.na(phenotype))
	phenotype <- phenotype[non_na]
	dgrp_name <- dgrp_name[non_na]
	dgrp_name <- gsub(dgrp_name, pattern = "DGRP_", replacement = "line_")
	dgrp_name <- gsub(dgrp_name, pattern = "line_0", replacement = "line_")
	dgrp_name <- gsub(dgrp_name, pattern = "line_0", replacement = "line_")
	results$results[[s]][["nb_used_dgrps"]] <- length(dgrp_name)
	
	## Generating phenotype file for PLINK
	plink.file_content <- data.frame(FID = dgrp_name, IID = dgrp_name)
	plink.file_content[phenotype_name] <- phenotype
	plink.file_content <- subset(plink.file_content, IID %in% dgrp_lines_genotyped) ## Restrict to genotyped samples
	results$results[[s]][["plink_phenotype_file"]] <- paste0(results$output_folder, phenotype_name, ".pheno.plink2.tsv")
	fwrite(x = format(plink.file_content, nsmall = 1), file = results$results[[s]][["plink_phenotype_file"]], sep = "\t", quote = F, row.names = F, col.names = T, scipen=50)
	
	## Running GWAS
	list_values <- unique(plink.file_content[,3])
	if(length(list_values) == 0){
		results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: No values.")
		message(results$results[[s]][["error_message"]])
	} else if(length(list_values) == 1){
		results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: Non-variable phenotype.")
		message(results$results[[s]][["error_message"]])
	} else if(length(list_values) == 2 && all(sort(list_values) == c(0, 2))){
		# Error: All samples for --glm phenotype 'S7_MB_32_M' are cases.
		results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: Values are only in [0, 2]. This is a specific Case/Control coding for PLINK2, can you change these values?")
		message(results$results[[s]][["error_message"]])
	}	else if(length(list_values) == 2 && all(sort(list_values) == c(0, 1))){
		results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: Values are only in [0, 1]. This is a specific Case/Control coding for PLINK2, can you change these values?")
		message(results$results[[s]][["error_message"]])
	} else{
		## Covariate analysis (Anova)
		results$results[[s]][["covariate_anova_file"]] <- paste0(results$output_folder, phenotype_name, ".cov.anova.txt")
		sink(results$results[[s]][["covariate_anova_file"]])
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
					results$results[[s]][["error_message"]] <- "ERROR: In Anova.III.lm(mod, error, singular.ok = singular.ok, ...) : there are aliased coefficients in the model"
					message(results$results[[s]][["error_message"]])
				} else {
					print(Anova(fit, type = 'III'))
				}
			}
			print(summary(fit))
		} else {
			results$results[[s]][["error_message"]] <- "ERROR: No invariant covariates??"
			message(results$results[[s]][["error_message"]])
		}
		sink() # returns output to the console
		
		if(length(covars) + 2 >= nrow(plink.file_content)){
			results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " [i = ", i, "]: # samples <= # predictor columns.")
			message(results$results[[s]][["error_message"]])
		} else {
			# PLINK 2.x 
			std.out <- NULL
			tryCatch(
			{
				std.out <<- system(paste0("plink2 --threads ",results$nb_threads," --glm hide-covar --geno ",results$genotype_ratio," --maf ",results$min_maf," --covar ", results$covariate_file," --bfile ", results$genotype_file," --pheno ",results$results[[s]]$plink_phenotype_file," --out  ", results$output_folder, "PLINK2"), intern = T, ignore.stderr = T, ignore.stdout = F)
			},
			warning=function(cond) {
				if(grepl(x = cond, pattern = "had status 7")) {
					results$results[[s]][["error_message"]] <- paste0("ERROR [HAD STATUS 7] in processing phenotype ", phenotype_name, " [i = ", i, "]: # samples <= # predictor columns.")
					message(results$results[[s]][["error_message"]])
				} else {
					# I redo it without catching the warnings... It's ugly I know
					std.out <<- system(paste0("plink2 --threads ",results$nb_threads," --glm hide-covar --geno ",results$genotype_ratio," --maf ",results$min_maf," --covar ", results$covariate_file," --bfile ", results$genotype_file," --pheno ",results$results[[s]]$plink_phenotype_file," --out  ", results$output_folder, "PLINK2"), intern = T, ignore.stderr = T, ignore.stdout = F)
				}
			})
			
			if(!is.null(std.out)) {
				file_name <- std.out[sapply(std.out, function(x) startsWith(x = x, prefix = "Results written to"))]
				## Rename output
				if(length(file_name) > 0)
				{
					file_name <- gsub(x = file_name, pattern = "Results written to ", replacement = "")
					file_name <- gsub(x = file_name, pattern = " \\.", replacement = "")
					file.rename(file_name, gsub(x = file_name, pattern = "/PLINK2\\.", replacement = "/"))
					file_name <- gsub(x = file_name, pattern = "/PLINK2\\.", replacement = "/")
					file.rename(paste0(results$output_folder, "PLINK2.log"), paste0(file_name, ".log"))
					results$results[[s]]$plink_full_output_file <- file_name
					results$results[[s]]$plink_full_log_file <- paste0(file_name, ".log")
				} else {
					results$results[[s]]$plink_full_output_file <- NA
					results$results[[s]]$plink_full_log_file <- NA
				}
				
				## Reading Output
				plink_results <- NULL
				col6 <- "BETA"
				results$results[[s]]$gwas_file_suffix <- ".glm.linear"
				if(endsWith(results$results[[s]]$plink_full_output_file, "glm.linear")){
					# Numeric / continuous
					plink_results <<- fread(results$results[[s]]$plink_full_output_file, data.table = F, nThread = results$nb_threads, showProgress = F)
					results$results[[s]]$gwas_model <- "GLM Linear"
					results$results[[s]]$gwas_file_suffix <- ".glm.linear"
					results$results[[s]]$gwas_type <- "continuous (quantitative)"
				} else if(endsWith(results$results[[s]]$plink_full_output_file, "glm.logistic.hybrid")){
					# Categorical
					plink_results <<- fread(results$results[[s]]$plink_full_output_file, data.table = F, nThread = results$nb_threads, showProgress = F)
					col6 <<- "OR"
					results$results[[s]]$gwas_model <- "GLM Logistic Hybrid"
					results$results[[s]]$gwas_file_suffix <- ".glm.logistic.hybrid"
					results$results[[s]]$gwas_type <- "categorical (qualitative)"
				} else {
					results$results[[s]][["error_message"]] <- paste0("ERROR? I don't recognize the output file:", results$results[[s]]$plink_full_output_file)
					message(results$results[[s]][["error_message"]])
				}
				
				if(!is.null(plink_results)){
					# Change column names and add FDR
					plink_results$P <- as.numeric(plink_results$P)
					plink_results <- subset(plink_results, TEST == "ADD")
					plink_results <- plink_results[!is.na(plink_results$P),c("#CHROM", "POS", "ID", "REF", "ALT", col6, "P")]
					plink_results$FDR_BH <- p.adjust(plink_results$P, method = "fdr")
					
					## QQPLOT before filtering
					results$results[[s]][["qq_plot"]] <- paste0(results$output_folder, phenotype_name, ".qqplot.png")
					plink_results$P[plink_results$P == 0] <- 1E-255 # For avoiding errors
					png(results$results[[s]][["qq_plot"]], height = 500, width = 500)
					qqPlotFast(plink_results$P, ci.level = NULL, col = "black", makelegend = F, lwd = 1)
					dev.off()
		
					## Manhattan plot
					results$results[[s]][["manhattan_plot"]] <- paste0(results$output_folder, phenotype_name, ".manhattan.png")
					chroms <- factor(plink_results$`#CHROM`)
					levels(chroms) <- c("2L", "2R", "3L", "3R", "X", "4")
					m <- manPlotPrepare(pvalues = plink_results$P, chr=chroms, pos = plink_results$POS)
					png(results$results[[s]][["manhattan_plot"]], height = 500, width = 500)
					manPlotFast(man = m, lwd = 1, colorSet = c("black", "darkgrey"))
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
					fwrite(x = plink_results, file = paste0(results$output_folder, phenotype_name, results$results[[s]]$gwas_file_suffix, ".top_0.01.annot.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
					system(paste0("gzip -f ", results$output_folder, phenotype_name, results$results[[s]]$gwas_file_suffix, ".top_0.01.annot.tsv"))
					results$results[[s]][["plink_filtered_output_file"]] <- paste0(results$output_folder, phenotype_name, results$results[[s]]$gwas_file_suffix, ".top_0.01.annot.tsv.gz")
				}
			}
		}
	}
}
message("Total time = ", toReadableTime(as.numeric(Sys.time() - start_time, unit = "secs")))

# Output
print(toJSON(results, pretty=TRUE, auto_unbox=TRUE))
