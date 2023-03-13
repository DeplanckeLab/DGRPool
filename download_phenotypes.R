##################################################
## Project: DGRPool
## Script purpose: Downloading phenotypes using JSON API
## Version: 1.0.0
## Date Created: 2022 Dec 22
## Date Modified: 2023 Mar 10
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Working directory
setwd("DGRPool/") # Root folder

# Libraries
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(data.table))

# Reading all studies
json_studies <- fromJSON("https://dgrpool.epfl.ch/studies.json") # Get it online from the website => Always most up-to-date version
json_studies <- json_studies[with(json_studies, order(id)),]
rownames(json_studies) <- json_studies$id
message(nrow(json_studies), " studies found")

# Reading all "standard" DGRP lines
data.dgrp_lines <- fread("dgrp2.fam")$V1

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
	message(length(dgrp_lines), " DGRP lines")
	data_pheno <- list() # List by sex
	if(length(dgrp_lines) > 0){
		# First, find all phenotypes across DGRP lines
		list_pheno_all <- c()
		for(dgrp in dgrp_lines){
			data_dgrp <- json_study$pheno_mean[[dgrp]]
			list_pheno_all <- unique(c(list_pheno_all, names(data_dgrp)[names(data_dgrp) != "sex"]))
		}
		if(length(list_pheno_all) >0){
			# Now find values for each DGRP line
			for(dgrp in dgrp_lines){
				data_dgrp <- json_study$pheno_mean[[dgrp]]
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
						phe[dgrp, p] <- data_dgrp[[p]][data_dgrp$sex == s]
					}
					data_pheno[[s]] <- phe
				}
			}
			
			# Match our data
			dgrp_lines <- gsub(x = dgrp_lines, pattern = "DGRP", replacement = "line")
			dgrp_lines <- gsub(x = dgrp_lines, pattern = "line_0", replacement = "line_")
			dgrp_lines <- gsub(x = dgrp_lines, pattern = "line_0", replacement = "line_")
			common_lines <- intersect(data.dgrp_lines, dgrp_lines)
			if(length(common_lines) != length(dgrp_lines)){
				message("[WARNING] Only ", length(common_lines), " lines intersect with the reference ", length(data.dgrp_lines), " lines. ", length(dgrp_lines) - length(common_lines), " line phenotypes are ignored.")	
			}
			for(sex in names(data_pheno))
			{
				phenotypes <- data_pheno[[sex]]
				phenotypes$dgrp <- gsub(x = phenotypes$dgrp, pattern = "DGRP", replacement = "line")
				rownames(phenotypes) <- phenotypes$dgrp
				data.all_pheno[[sex]][phenotypes$dgrp,paste0("S", sid, "_", colnames(phenotypes)[colnames(phenotypes) != "dgrp"])] <- phenotypes[,colnames(phenotypes) != "dgrp"]
			}
			message("Available sex (",length(names(data_pheno)), "): [", paste(names(data_pheno), collapse = ", "), "]")
			message("Available phenotypes (",ncol(data_pheno[[1]]) - 1, "): [", paste(colnames(data_pheno[[1]])[-1], collapse = ", "), "]\n")
		} else { message("No phenotype available for this study\n") }
	} else { message("No DGRP available for this study\n") }
}

## Save object
saveRDS(data.all_pheno, file = "data.all_pheno_10_03_03.rds")
