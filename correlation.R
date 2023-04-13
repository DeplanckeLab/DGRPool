##################################################
## Project: DGRPool
## Script purpose: Compute correlation stats for the correlation tool
## Version: 1.0.0
## Date Created: 2023 Mar 02
## Date Modified: 2023 Mar 02
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################
# Global script options
options(echo = F)
args <- commandArgs(trailingOnly = TRUE)
# Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
# Functions
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
# Parameters
if(is.null(args)) error.json("You need to input some parameters")
if(length(args) != 1) error.json("You need to input one parameter: Phenotype file")
pheno_file <- args[1]
col1 <- 3
col2 <- 4
# Read phenotype files
if(!file.exists(pheno_file)) error.json(paste0("This file: '",pheno_file,"' does not exist!"))
data.dgrp_pheno <- fread(pheno_file, sep = "\t", header = T, data.table = F, na.strings = "")
if(colnames(data.dgrp_pheno)[1] != "DGRP") error.json(paste0("First column of ",pheno_file," should be the 'DGRP' column"))
if(colnames(data.dgrp_pheno)[2] != "sex") error.json(paste0("Second column of ",pheno_file," should be the 'sex' column"))
if(ncol(data.dgrp_pheno) != 4) error.json(paste0("There should be 4 columns in the phenotype file, there is ", ))
# Compute stats
results <- list()
for(sex1 in unique(data.dgrp_pheno$sex)){
    for(sex2 in unique(data.dgrp_pheno$sex)){
        # Results for this round
        res <- list()
        
        # Restrict dataset
        data.pheno1 <- data.dgrp_pheno[data.dgrp_pheno$sex == sex1, c(1, col1)]
        data.pheno1 <- data.pheno1[!is.na(data.pheno1[,2]),]
        data.pheno2 <- data.dgrp_pheno[data.dgrp_pheno$sex == sex2, c(1, col2)]
        data.pheno2 <- data.pheno2[!is.na(data.pheno2[,2]),]
        common.dgrp <- intersect(data.pheno1$DGRP, data.pheno2$DGRP)
        res$common.notNA.dgrp <- length(common.dgrp)
        if(res$common.notNA.dgrp < 3){
            res$pearson_cor <- NA
            res$pearson_p <- NA
            res$spearman_cor <- NA
            res$spearman_p <- NA
            res$lm_b <- NA
            res$lm_a <- NA
            res$lm_p <- NA
            res$lm_r_squared <- NA
            res$lm_adj_r_squared <- NA
            res$warning <- "Not enough overlapping non-NA samples to test"
        } else {
            # Subset
            rownames(data.pheno1) <- data.pheno1$DGRP
            rownames(data.pheno2) <- data.pheno2$DGRP
            data.pheno1 <- data.pheno1[common.dgrp, 2]
            data.pheno2 <- data.pheno2[common.dgrp, 2]
            
            # Pearson's correlation
            suppressWarnings(pearson_test <- cor.test(data.pheno1, data.pheno2, method = "pearson", use = "complete.obs"))
            res$pearson_cor <- pearson_test$estimate
            res$pearson_p <- pearson_test$p.value
            
            # Spearman's correlation
            suppressWarnings(spearman_test <- cor.test(data.pheno1, data.pheno2, method = "spearman", use = "complete.obs"))
            res$spearman_cor <- spearman_test$estimate
            res$spearman_p <- spearman_test$p.value
            
            # Linear regression
            lm_test <- lm(data.pheno2 ~ data.pheno1)
            res$lm_b <- lm_test$coefficients[1]
            res$lm_a <- lm_test$coefficients[2]
            lm_test <- summary(lm_test)
            res$lm_p <- lm_test$coefficients[2, "Pr(>|t|)"]
            res$lm_r_squared <- lm_test$r.squared
            res$lm_adj_r_squared <- lm_test$adj.r.squared
        }
        results[[paste0(sex1, "_", sex2)]] <- res
    }
}
# Output
print(toJSON(results, pretty=TRUE, auto_unbox=TRUE))


