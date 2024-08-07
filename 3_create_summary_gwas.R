##################################################
## Project: DGRPool
## Script purpose: Stats over GWAS computation
## Version: 1.0.0
## Date Created: 2023 Mar 24
## Date Modified: 2024 Aug 02
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

## Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(require(doSNOW)) # for parallelization
suppressPackageStartupMessages(require(dplyr)) # for fast filtering

## Parameters
nb_cores <- 24
dgrp2_vcf_file <- "../../dgrp2.filtered.fastlmm.vcf"
dgrp2_annotation_file <- "../../dgrp.fb557.annot.txt.gz"
gwas_folder <- "../../latest/"
output_file_prefix <- "resources/data.gwas_24_08_02"

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
data.vcf <- fread(dgrp2_vcf_file, data.table = F) # 4438427 variants (or 1882162 if filtered VCF file)
data.vcf <- data.vcf[,c("#CHROM", "POS", "ID", "REF", "ALT")]
colnames(data.vcf) <- c("chr", "pos", "id", "ref", "alt")
annotation.dr <- fread (dgrp2_annotation_file, header = FALSE, sep = "\t")
annotation.dr <- annotation.dr[,c("V1", "V3", "V4")]
colnames(annotation.dr) <- c("SNP", "gene_annotation", "regulatory_annotation")
data.gwas <- merge(data.vcf, annotation.dr, by.x = "id", by.y = "SNP")
rm(data.vcf) # Clean up RAM
rm(annotation.dr) # Clean up RAM
rownames(data.gwas) <- data.gwas$id

## Read through GWAS data
# Get all files
gwas.files <- list.files(path = gwas_folder, pattern = "*glm*", full.names = T)
message(length(gwas.files), " GWAS files to process")

# Read all files and summarize (parallelized)
cl <- makeCluster(nb_cores) # create the cluster
registerDoSNOW(cl) #register it to be used by %dopar%
message("Computing summary GWAS on ", nb_cores, " cores in parallel...")
start_time = Sys.time()
progress.bar <- txtProgressBar(style = 3, min = 1, max = length(gwas.files))
data.sum <- foreach(i = 1:length(gwas.files), .packages = c("data.table"), .combine = c, .options.snow = list(progress = function(x) setTxtProgressBar(progress.bar, x))) %dopar% {
  # Path to read
  gwas.file <- gwas.files[i]

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
  data.gwas.tmp <- list()
  data.gwas.tmp[[gwas.name]] <- values[data.gwas$id]
  return(data.gwas.tmp)
}
close(progress.bar)
stopCluster(cl)
message("Total time = ", toReadableTime(as.numeric(Sys.time() - start_time, unit = "secs"))) # 16mn 7s with 24 cores

# Add data to the global data.gwas object
for(pheno in names(data.sum)){
  data.gwas[pheno] <- data.sum[pheno]
}
rm(data.sum) # Clean up RAM

# Save object as compressed tsv (takes ~30s with fwrite, >10 times faster than rds objects)
fwrite(data.gwas, file = paste0(output_file_prefix, ".tsv.gz"), sep = "\t", row.names = T, col.names = T, quote = F)

# Load object from file (~1mn30s)
#data.gwas <- fread(paste0(output_file_prefix, ".tsv.gz"), sep = "\t", data.table = F)
#rownames(data.gwas) <- data.gwas$V1
#data.gwas <- data.gwas[,-1]

# Filter variants with only NAs
start_time = Sys.time()
data.gwas <- filter(data.gwas, !if_all(8:ncol(data.gwas), is.na)) # 1876017 variants x 855 phenotypes
message("Total time = ", toReadableTime(as.numeric(Sys.time() - start_time, unit = "secs"))) # Total time = 16s

# Save object as compressed tsv (takes ~30s with fwrite, >10 times faster than rds objects)
fwrite(data.gwas, file = paste0(output_file_prefix, ".filtered.tsv.gz"), sep = "\t", row.names = T, col.names = T, quote = F)
