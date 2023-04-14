##################################################
## Project: DGRPool
## Script purpose: Stats over GWAS computation
## Version: 1.0.0
## Date Created: 2022 Mar 24
## Date Modified: 2023 Mar 24
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################
start_time = Sys.time()
## Working directory
setwd("/data/gardeux/DGRPool/")

## Libraries
suppressPackageStartupMessages(library(data.table))

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

## Prepare a huge table?
## Load annotation
data.vcf <- fread("dgrp2.vcf", data.table = F) # 4438427 variants
data.vcf <- data.vcf[,c("#CHROM", "POS", "ID", "REF", "ALT")]
colnames(data.vcf) <- c("chr", "pos", "id", "ref", "alt")
annotation.dr <- fread ("/data/genome/drosophila_melanogaster/dm3/dgrp.fb557.annot.txt", header = FALSE, sep = "\t")
annotation.dr <- annotation.dr[,c("V1", "V3", "V4")]
colnames(annotation.dr) <- c("SNP", "gene_annotation", "regulatory_annotation")
data.gwas <- merge(data.vcf, annotation.dr, by.x = "id", by.y = "SNP")
rm(data.vcf)
rm(annotation.dr)
rownames(data.gwas) <- data.gwas$id

## Read through GWAS data
gwas.files <- list.files(path = "GWAS_full_results", pattern = "*glm*", full.names = T)
message(length(gwas.files), " GWAS files to process")
pb <- txtProgressBar(min = 0, max = length(gwas.files), style = 3, width = 50, char = "=")
cnt <- 0
for(gwas.file in gwas.files){
	# Parse name from file path
	gwas.name <- gsub(x = gwas.file, pattern = ".*/", replacement = "")
	gwas.type = substr(gwas.name, start = unlist(gregexpr(text = gwas.name, pattern = "\\."))[1] + 1, stop = nchar(gwas.name))
	gwas.name <- gsub(x = gwas.name, pattern = "\\..*", replacement = "")
	# Read data
	data.tmp <- fread(gwas.file, sep = "\t", data.table = F)
	rows <- data.tmp$ID
	values <- data.tmp$P
	rm(data.tmp)
	names(values) <- rows
	# Add column to main data frame
	data.gwas[[gwas.name]] <- values[data.gwas$id]
	
	## Progression bar
	cnt <- cnt + 1
	setTxtProgressBar(pb, cnt)
}
message("Total time = ", toReadableTime(as.numeric(Sys.time() - start_time, unit = "secs"))) # Total time = 38mn 12s

# Save object as compressed tsv (>10 times faster than rds objects)
fwrite(data.gwas, file = "data.gwas_24_03_23.tsv.gz", sep = "\t", row.names = T, col.names = T, quote = F)

# Filter snps with only NAs
start_time = Sys.time()
s <- rep(NA, nrow(data.gwas))
for(i in 1:nrow(data.gwas)){
	s[i] <- sum(is.na(data.gwas[i,9:ncol(data.gwas)]))
}
data.gwas <- data.gwas[s != 846,]
message("Total time = ", toReadableTime(as.numeric(Sys.time() - start_time, unit = "secs"))) # Total time = 12h 32mn 6s !!

# Save object as compressed tsv
fwrite(data.gwas, file = "data.gwas_24_03_23.filtered.tsv.gz", sep = "\t", row.names = T, col.names = T, quote = F)

#start_time = Sys.time()
#check.na <- apply(data.gwas, 1, function(x) sum(is.na(x)))
#message("Total time = ", toReadableTime(as.numeric(Sys.time() - start_time, unit = "secs"))) # Total time = 2mn 30s

#start_time = Sys.time()
#check.na <- apply(data.gwas[,9:ncol(data.gwas)], 1, function(x) sum(is.na(x)))
#message("Total time = ", toReadableTime(as.numeric(Sys.time() - start_time, unit = "secs"))) # Total time = 2mn 30s

