##################################################
## Project: DGRPool
## Script purpose: Running GWAS for user-submitted data
## Version: 1.0.0
## Date Created: 2023 Mar 29
## Date Modified: 2024 Aug 06
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Measure execution time
start_time = Sys.time()

# Global script options
options(echo = F)
args <- commandArgs(trailingOnly = TRUE)

## Libraries
suppressPackageStartupMessages(library(data.table)) # Fast reading of files
suppressPackageStartupMessages(library(stringr)) # String utilities
suppressPackageStartupMessages(library(ramwas)) # QQ-Plot & Manhattan plot
suppressPackageStartupMessages(library(car)) # Covariate analysis
suppressPackageStartupMessages(library(rstatix)) # Shapiro test
suppressPackageStartupMessages(library(Cairo)) # For the - in front of -log10(P) y-axis of the PDF outputs
suppressPackageStartupMessages(library(jsonlite)) # JSON output

## Functions
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
if(is.null(args) | length(args) < 9){
	message("Usage: running_GWAS_user.R parameters\n1. [Required] Path to phenotype file\n2. [Required] Path to PLINK genotype file prefix (bed)\n3. [Required] Path to covariate file\n4. [Required] Path to annotation file\n5. [Required] Output folder path\n6. [Required] Number of threads\n7. [Required] Minimum genotyping ratio\n8. [Required] Minimum MAF\n9. [Required] Filtering threshold for p-value\n10. [Optional] Name of phenotype to load (if file contains multiple phenotypes)")
	error.json("You need to input at least 9 parameters: [Phenotype file path][Genotype file prefix][Covariate file path][Annotation file path][Output folder path][Number of threads][Genotype Ratio][Min MAF]")
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
results$gwas_options <- paste0("--glm hide-covar --quantile-normalize --variance-standardize --geno ", results$genotype_ratio, " --maf ", results$min_maf) # --geno 0.2 --maf 0.01

results$filtering_pval <- as.numeric(args[9])

results$pheno_name <- NULL
results$pheno_column <- 3 # Default phenotype to read
results$which_sex <- "NA" # Default, in case no information
if(length(args) == 10){
	# Check if there is a "name" parameter
	results$pheno_name <- args[10]
}
if(length(args) > 10) error.json("You entered too many parameters: ", length(args))

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
	phenotype_name <- gsub(phenotype_name, pattern = "[^a-zA-Z0-9]+", replacement = "_")
	
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
	
	## Test for categorical phenotype
	is_categorical <- F
	tryCatch({
	  tmp.var <- as.numeric(plink.file_content[[phenotype_name]])
	}, warning=function(cond) {
	  if(grepl(x = cond, pattern = "NAs introduced by coercion")) {
	    is_categorical <<- T
	  }
	})
	
	## Running GWAS
	list_values <- unique(plink.file_content[,3])
	if(is_categorical){
	  results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " Categorical phenotype?")
	  message(results$results[[s]][["error_message"]])
	} else if(length(list_values) == 0){
		results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " : No values.")
		message(results$results[[s]][["error_message"]])
	} else if(length(list_values) == 1){
		results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " : Non-variable phenotype.")
		message(results$results[[s]][["error_message"]])
	} else if(length(list_values) == 2 && all(sort(list_values) == c(0, 2))){
		# Error: All samples for --glm phenotype 'S7_MB_32_M' are cases.
		results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " : Values are only in [0, 2]. This is a specific Case/Control coding for PLINK2, can you change these values?")
		message(results$results[[s]][["error_message"]])
	}	else if(length(list_values) == 2 && all(sort(list_values) == c(0, 1))){
		results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " : Values are only in [0, 1]. This is a specific Case/Control coding for PLINK2, can you change these values?")
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
		
		# Compute the normality test of Shapiro-Wilk
		cat("Shapiro-Wilk test of normality:\n")
		tryCatch( { 
		  shapiro <- as.data.frame(shapiro_test(residuals(fit))) 
		  cat(paste0('{\"shapiro\"={"statistic":', shapiro$statistic, ',"pvalue":', shapiro$p.value, '}}\n'))
		  if(shapiro$p.value > 0.05) cat("Interpretation: NORMAL OK (p > 0.05)\n")
		  if(shapiro$p.value <= 0.05) cat("Interpretation: NOT NORMAL (p <= 0.05)\n")
		}, error = function(e) {cat("Shapiro test failed : all residuals values are identical\n")})
		
		# Kruskal-Wallis test
		cat("\nKruskal-Wallis test\n")
		kruskal.results <- data.frame(covariate = "TOTO", chisquared = -1,  df = -1, pvalue = -1)
		for(cova in covars){
		  res <- kruskal.test(as.formula(paste0(phenotype_name, " ~ ", cova)), data = data.cov)
		  kruskal.results <- rbind(kruskal.results, data.frame(covariate = cova, chisquared = res$statistic, df = res$parameter, pvalue = res$p.value))
		}
		kruskal.results <- kruskal.results[-1,]
		rownames(kruskal.results) <- NULL
		cat(jsonlite::toJSON(kruskal.results))
		
		# returns output to the console
		sink()
		
		# Process phenotypes
		if(length(covars) + 2 >= nrow(plink.file_content)){
			results$results[[s]][["error_message"]] <- paste0("ERROR in processing phenotype ", phenotype_name, " : # samples <= # predictor columns.")
			message(results$results[[s]][["error_message"]])
		} else {
			# Run PLINK 2.x 
			std.out <- NULL
			results$command_line_gwas <- paste0("plink2 --threads ", results$nb_threads, " ", results$gwas_options, " --covar ", results$covariate_file, " --bfile ", results$genotype_file, " --pheno \"", results$results[[s]]$plink_phenotype_file, "\" --out \"", results$output_folder, "PLINK2\"")
			tryCatch(
			{
			  std.out <<- system(results$command_line_gwas, intern = T, ignore.stderr = T, ignore.stdout = F)
			},
			warning=function(cond) {
				if(grepl(x = cond, pattern = "had status 7")) {
					results$results[[s]][["error_message"]] <- paste0("ERROR [HAD STATUS 7] in processing phenotype ", phenotype_name, " : # samples <= # predictor columns.")
					message(results$results[[s]][["error_message"]])
				} else {
					# I redo it without catching the warnings... It's ugly I know
				  std.out <<- system(results$command_line_gwas, intern = T, ignore.stderr = T, ignore.stdout = F)
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
					results$results[[s]][["qq_plot_pdf"]] <- paste0(results$output_folder, phenotype_name, ".qqplot.pdf")
					plink_results$P[plink_results$P == 0] <- 1E-255 # For avoiding errors
					png(results$results[[s]][["qq_plot"]], height = 500, width = 500)
					qqPlotFast(plink_results$P, ci.level = NULL, col = "black", makelegend = F, lwd = 1, newplot = T)
					dev.off()
		
					# SVG is way too big. So let's go PDF
					# I recompute the plot... I know it's not optimal but I don't know how to do it else
					CairoPDF(results$results[[s]][["qq_plot_pdf"]], height = 7, width = 7)
					qqPlotFast(plink_results$P, ci.level = NULL, col = "black", makelegend = F, lwd = 1, newplot = T)
					dev.off()
					
					## Manhattan plot
					results$results[[s]][["manhattan_plot"]] <- paste0(results$output_folder, phenotype_name, ".manhattan.png")
					results$results[[s]][["manhattan_plot_pdf"]] <- paste0(results$output_folder, phenotype_name, ".manhattan.pdf")
					chroms <- factor(plink_results$`#CHROM`)
					levels(chroms) <- c("2L", "2R", "3L", "3R", "X", "4")
					m <- manPlotPrepare(pvalues = plink_results$P, chr=chroms, pos = plink_results$POS)
					png(results$results[[s]][["manhattan_plot"]], height = 500, width = 500)
					manPlotFast(man = m, lwd = 1, colorSet = c("black", "darkgrey"))
					dev.off()
					
					# SVG is way too big. So let's go PDF
					# I recompute the plot... I know it's not optimal but I don't know how to do it else
					CairoPDF(results$results[[s]][["manhattan_plot_pdf"]], height = 7, width = 7)
					manPlotFast(man = m, lwd = 1, colorSet = c("black", "darkgrey"))
					dev.off()
					
					## Filtered top results and annotate them
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
					
					# Writing full annotated file
					fwrite(x = plink_results, file = paste0(results$output_folder, phenotype_name, results$results[[s]]$gwas_file_suffix, ".fdr.annot.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
					system(paste0("gzip -f \"", results$output_folder, phenotype_name, results$results[[s]]$gwas_file_suffix, ".fdr.annot.tsv\""))
					results$results[[s]][["plink_output_file"]] <- paste0(results$output_folder, phenotype_name, results$results[[s]]$gwas_file_suffix, ".fdr.annot.tsv.gz")
					
					# Filtering
					plink_results <- subset(plink_results, P <= results$filtering_pval)
					fwrite(x = plink_results, file = paste0(results$output_folder, phenotype_name, results$results[[s]]$gwas_file_suffix, ".top_", results$filtering_pval, ".annot.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
					system(paste0("gzip -f \"", results$output_folder, phenotype_name, results$results[[s]]$gwas_file_suffix, ".top_", results$filtering_pval, ".annot.tsv\""))
					results$results[[s]][["plink_filtered_output_file"]] <- paste0(results$output_folder, phenotype_name, results$results[[s]]$gwas_file_suffix, ".top_", results$filtering_pval, ".annot.tsv.gz")
					
					# Run gene enrichment script on Flybase phenotypes (PYTHON)
					#system(paste0("python ", input_gene_enrichment_script, " ", output_gwas_directory, " ", input_flybase_phenotype_file, " ", input_ncsu_annot_file, " ", fstudy, " ", fphenotype, " ", fsex, " ", input_filtering_pval, " ", output_gwas_directory, "/", phenotype_name, file_suffix, ".top_", input_filtering_pval, ".gene_enrichment.phenotypes.tsv"))
					
					# Run gene enrichment script on GO terms (PYTHON)
					#system(paste0("python ", input_gene_enrichment_script, " ", output_gwas_directory, " ", input_go_file, " ", input_ncsu_annot_file, " ", fstudy, " ", fphenotype, " ", fsex, " ", input_filtering_pval, " ", output_gwas_directory, "/", phenotype_name, file_suffix, ".top_", input_filtering_pval, ".gene_enrichment.go.tsv"))
				}
			}
		}
	}
}

# Output
print(toJSON(results, pretty=TRUE, auto_unbox=TRUE))
