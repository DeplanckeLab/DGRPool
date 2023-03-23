##################################################
## Project: DGRPool
## Script purpose: Downloading phenotypes using JSON API
## Version: 1.0.0
## Date Created: 2022 Dec 22
## Date Modified: 2023 Mar 21
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Working directory
setwd("DGRPool/")

# Libraries
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(data.table))

# Reading all studies
json_studies <- fromJSON("https://dgrpool.epfl.ch/studies.json")
json_studies <- json_studies[with(json_studies, order(id)),]
rownames(json_studies) <- json_studies$id
message(nrow(json_studies), " studies found")

# Reading all phenotypes
json_phenotypes <- fromJSON("https://dgrpool.epfl.ch/phenotypes.json?all=1")
json_phenotypes <- json_phenotypes[with(json_phenotypes, order(id)),]
rownames(json_phenotypes) <- json_phenotypes$id
message(nrow(json_phenotypes), " phenotypes found")

# Reading all "standard" DGRP lines
data_dgrps <- fread("https://dgrpool.epfl.ch/studies/1/get_file?name=dgrp_lines.tsv&namespace=downloads", sep = "\t", data.table = F)
data.dgrp_lines <- data_dgrps$dgrp
message(length(unique(data.dgrp_lines)), " DGRP lines found")

total_phenotype_number <- 0
nb_issues <- c()
nb_obsolete <- c()

# Reading each study
tmp <- data.frame(dgrp = data.dgrp_lines)
rownames(tmp) <- data.dgrp_lines
data.all_pheno <- list("F" = tmp, "M" = tmp, "NA" = tmp)
for(sid in json_studies$id) {
	sub <- subset(json_studies, id == sid)
	message("Study id = ", sid, "\ncreated at ", sub$created_at, ", modified at ", sub$updated_at)
	json_study <- fromJSON(sub$url)
	if(json_study$created_at != sub$created_at) stop("ERROR created at")
	if(json_study$updated_at != sub$updated_at) stop("ERROR updated at")
	message(json_study$authors, ", ", json_study$title, ", ", json_study$journal_id, " ", json_study$volume, "(", json_study$issue, "), ", json_study$year)
	dgrp_lines <- names(json_study$pheno_mean)
	#dgrp_lines <- names(json_study$pheno_sum)
	message(length(dgrp_lines), " DGRP lines")
	data_pheno <- list() # List by sex
	if(length(dgrp_lines) > 0){
		# First, find all phenotypes across DGRP lines
		list_pheno_all <- c()
		for(dgrp in dgrp_lines){
			data_dgrp <- json_study$pheno_mean[[dgrp]]
			#data_dgrp <- json_study$pheno_sum[[dgrp]]
			list_pheno_all <- unique(c(list_pheno_all, names(data_dgrp)[names(data_dgrp) != "sex"]))
		}
		if(length(list_pheno_all) >0){
			# Now find values for each DGRP line
			for(dgrp in dgrp_lines){
				data_dgrp <- json_study$pheno_mean[[dgrp]]
				#data_dgrp <- json_study$pheno_sum[[dgrp]]
				list_pheno <- names(data_dgrp)[names(data_dgrp) != "sex"]
				if(!all(list_pheno %in% list_pheno_all)) stop("ERROR PHENO")
				data_dgrp$sex[is.na(data_dgrp$sex)] <- "NA" # For avoiding issues
				for(s in data_dgrp$sex){
					phe <- data_pheno[[s]]
					if(is.null(phe)){
						phe <- data.frame(matrix(nrow = length(dgrp_lines), ncol = length(list_pheno) + 1))
						colnames(phe) <- c("dgrp", list_pheno)
						phe$dgrp <- dgrp_lines
						rownames(phe) <- dgrp_lines
					}
					for(p in list_pheno){
						data.subpheno <- subset(json_phenotypes, id == p)
						if(nrow(data.subpheno) == 0){
							nb_issues <- c(nb_issues, p)
						} else {
							if(data.subpheno$obsolete){
								nb_obsolete <- c(nb_obsolete, p)
							} else {
								phe[dgrp, p] <- data_dgrp[[p]][data_dgrp$sex == s]
							}
						}
					}
					data_pheno[[s]] <- phe
				}
			}
			
			# Match our data
			#dgrp_lines <- gsub(x = dgrp_lines, pattern = "DGRP", replacement = "line")
			#dgrp_lines <- gsub(x = dgrp_lines, pattern = "line_0", replacement = "line_")
			#dgrp_lines <- gsub(x = dgrp_lines, pattern = "line_0", replacement = "line_")
			#common_lines <- intersect(data.dgrp_lines, dgrp_lines)
			#if(length(common_lines) != length(dgrp_lines)){
			#	message("[WARNING] Only ", length(common_lines), " lines intersect with the reference ", length(data.dgrp_lines), " lines. ", length(dgrp_lines) - length(common_lines), " line phenotypes are ignored.")	
			#}
			for(sex in names(data_pheno))
			{
				phenotypes <- data_pheno[[sex]]
				#phenotypes$dgrp <- gsub(x = phenotypes$dgrp, pattern = "DGRP", replacement = "line")
				rownames(phenotypes) <- phenotypes$dgrp
				data.all_pheno[[sex]][phenotypes$dgrp,paste0("S", sid, "_", colnames(phenotypes)[colnames(phenotypes) != "dgrp"])] <- phenotypes[,colnames(phenotypes) != "dgrp"]
			}
			message("Available sex (",length(names(data_pheno)), "): [", paste(names(data_pheno), collapse = ", "), "]")
			message("Available phenotypes (",ncol(data_pheno[[1]]) - 1, "): [", paste(colnames(data_pheno[[1]])[-1], collapse = ", "), "]\n")
			total_phenotype_number <- total_phenotype_number + (ncol(data_pheno[[1]]) - 1)
		} else { message("No phenotype available for this study\n") }
	} else { message("No DGRP available for this study\n") }
}

## Stats
nb_obsolete <- unique(nb_obsolete)
nb_issues <- unique(nb_issues)
message("Finished. Total number of phenotypes loaded: ", total_phenotype_number)
message("Finished. Total number of phenotypes OBSOLETE: ", length(nb_obsolete))
message("Finished. Total number of phenotypes WITH ISSUES: ", length(nb_issues))

## Filter NA poop
for(sex in c("M", "F", "NA")){
	to.keep <- apply(data.all_pheno[[sex]], 2, function(x) sum(is.na(x))) != nrow(data.all_pheno[[sex]])
	data.all_pheno[[sex]] <- data.all_pheno[[sex]][,to.keep]
	message("Sex ", sex, ": ", ncol(data.all_pheno[[sex]]) - 1, " phenotypes after removing NAs")
}

## Check again
phenotypes.all <- unique(c(colnames(data.all_pheno[["F"]])[2:ncol(data.all_pheno[["F"]])], colnames(data.all_pheno[["M"]])[2:ncol(data.all_pheno[["M"]])], colnames(data.all_pheno[["NA"]])[2:ncol(data.all_pheno[["NA"]])]))
message("Number phenotypes = ", length(phenotypes.all))
message("Number phenotypes (by sex) = ", ncol(data.all_pheno[["M"]]) + ncol(data.all_pheno[["F"]]) + ncol(data.all_pheno[["NA"]]) - 3)

message("Curated studies")
message(sum(limma::strsplit2(x = phenotypes.all, split = "_")[,1] %in% paste0("S", 1:41)), " phenotypes with associated data")
message("- ", sum(limma::strsplit2(x = colnames(data.all_pheno[["M"]]), split = "_")[,1] %in% paste0("S", 1:41)), " phenotypes with male data")
message("- ", sum(limma::strsplit2(x = colnames(data.all_pheno[["F"]]), split = "_")[,1] %in% paste0("S", 1:41)), " phenotypes with female data")
message("- ", sum(limma::strsplit2(x = colnames(data.all_pheno[["NA"]]), split = "_")[,1] %in% paste0("S", 1:41)), " phenotypes with undefined sex data")

## Save object
saveRDS(data.all_pheno, file = "RDS/data.all_pheno_21_03_23_filtered.rds")

