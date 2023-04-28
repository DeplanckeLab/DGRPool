##################################################
## Project: DGRPool
## Script purpose: Extreme analysis phenotype
## Version: 1.0.0
## Date Created: 2023 Mar
## Date Modified: 2023 Mar
## Author: Roel Bevers (contact Deplancke Lab)
##################################################

# Set environment ----------------------------------------- R version 4.0.3 (2020-10-10)

library(tidyverse)#2.0.0
library(data.table)#1.14.2
library(tidylog)#1.0.2
library(jsonlite)#1.8.0
library(ggplot2)#3.3.5
library(patchwork)#1.1.2

# DGRPool API ---------------------------------------------
json_studies <- fromJSON("https://dgrpool.epfl.ch/studies.json") %>% 
  as_tibble() %>% 
  arrange(id) 
message(nrow(json_studies), " studies found")

json_phenotypes <- fromJSON("https://dgrpool.epfl.ch/phenotypes.json?all=1") %>% 
  as_tibble() %>% 
  arrange(id) 
message(nrow(json_phenotypes), " phenotypes found")

# Load objects --------------------------------------------

base_all_pheno <- readRDS(file="RDS/data.all_pheno_21_03_23_filtered.rds")
names(base_all_pheno)

phenotypes <- read_tsv("resources/phenotypes.tsv", col_names = T) %>% 
  mutate(phenotype_tag = paste("S",paste(study_id, phenotype_id, sep = "_"), sep =""))
dgrp_lines <- read_tsv("resources/dgrp_lines.tsv", col_names = T)
studies <- read_tsv("resources/studies.tsv", col_names = T)

# Set parameters ------------------------------------------ 

extreme_threshold <- 0.15
extreme_categories <- c("Life history traits","Immunity","Toxicity",
                        "Resistance","Fecundity","Aging")

# Filter phenotypes for categories and curated only -------

curated_studies <- studies %>% 
  filter(status == "integrated") %>% 
  filter(str_detect(categories, paste0(extreme_categories,sep ="", collapse="|")))

# Retrieve phenotype date and reformat --------------------
tmpF <- base_all_pheno$`F` %>% 
  as_tibble() %>% 
  mutate_all(as.character) %>% 
  pivot_longer(-dgrp,
               names_to = "phenotype_tag",
               values_to = "phenotype_value") %>% 
  add_column(sex = "F")
tmpM <- base_all_pheno$`M` %>% 
  as_tibble() %>% 
  mutate_all(as.character) %>% 
  pivot_longer(-dgrp,
               names_to = "phenotype_tag",
               values_to = "phenotype_value") %>% 
  add_column(sex = "M")
tmpNA <- base_all_pheno$`NA` %>% 
  as_tibble() %>% 
  mutate_all(as.character) %>% 
  pivot_longer(-dgrp,
               names_to = "phenotype_tag",
               values_to = "phenotype_value") %>% 
  add_column(sex = "NA")

phenotype_pivot <- bind_rows(tmpF, tmpM, tmpNA) %>% 
  separate(phenotype_tag, 
           into = c("study_id","phenotype_id"), 
           sep = "_", 
           remove = F) %>% 
  mutate(study_id = str_extract(study_id, "[:digit:]+")) %>% 
  mutate(study_id = as.double(study_id),
         phenotype_id = as.double(phenotype_id)) %>% 
  filter(is.na(phenotype_value)==F)

phenotype_pivot_curated <- phenotype_pivot %>% 
  filter(study_id %in% curated_studies$study_id) %>% 
  left_join(phenotypes %>% 
              select(-study_id, -phenotype_id),
            by = "phenotype_tag") 

# Calculate ranks and extremes ----------------------------

phenotype_pivot_curated_ranked <- phenotype_pivot_curated %>% 
    mutate(phenotype_tag_sex = paste(phenotype_tag, sex, sep = "_")) %>% 
    filter(is_numeric == TRUE & is_continuous == TRUE) %>% 
    group_by(phenotype_tag_sex) %>% 
      mutate(phenotype_value = as.double(phenotype_value)) %>% 
      arrange(phenotype_value) %>% 
      add_column(rank = 1:nrow(.)) %>% 
      mutate(rank_real = min_rank(phenotype_value)) %>% 
      mutate(rank_thes_base_cutoff = ceiling(max(rank_real)*extreme_threshold),
             rank_thes_prop = rank_real/max(rank_real),
             rank_thes_eval = case_when(
               rank_real <= rank_thes_base_cutoff ~ -1,
               rank_real >= max(rank_real)-rank_thes_base_cutoff ~ 1,
               TRUE ~ 0
             )) %>% 
      mutate(z_score = (phenotype_value - mean(phenotype_value))/sd(phenotype_value)) %>% 
      mutate(pnorm_score = pnorm(z_score)) %>% 
      mutate(extreme_eval = case_when(
        pnorm_score < extreme_threshold ~ -1,
        pnorm_score > 1-extreme_threshold ~ 1,
        TRUE ~ 0
      )) %>% 
      mutate(prop_scal = (phenotype_value - min(phenotype_value))/(max(phenotype_value)-min(phenotype_value))) %>% 
      mutate(extreme_eval_prop_scal = case_when(
        prop_scal < extreme_threshold ~ -1,
        prop_scal > 1-extreme_threshold ~ 1,
        TRUE ~ 0)) %>% 
      mutate(value_diversity = n_distinct(phenotype_value)) %>% 
      mutate(shapiro_pval = shapiro.test(phenotype_value)[[2]]) %>% 
      mutate(shapiro_caution = case_when(
        shapiro_pval < 0.05 ~ "CAUTION_NON_NORMAL_DISTRIBUTION",
        TRUE ~ "NORMAL_DISTRIBUTION"
      )) %>% 
      select(dgrp,phenotype_tag_sex, phenotype_tag, sex,
             phenotype_value, contains("nber_sex"), value_diversity, shapiro_pval, shapiro_caution,
             rank, 
             rank_real, rank_thes_base_cutoff, rank_thes_prop, rank_thes_eval, 
             z_score, pnorm_score, extreme_eval, 
             prop_scal, extreme_eval_prop_scal) %>% 
    ungroup()

save(phenotype_pivot_curated_ranked, file = "Extremes/phenotype_pivot_curated_ranked.RData")

# Fraction of extreme/medicore phenotypes -----------------

# Following the rank percentages
extreme_dgrp_fraction_thes <- phenotype_pivot_curated_ranked %>% 
  group_by(dgrp, sex) %>% 
  count(rank_thes_eval) %>% 
  mutate(total_nr_phenotypes = sum(n),
         fraction_of_group = n/total_nr_phenotypes) %>% 
  ungroup %>% 
  mutate(extreme_class = case_when(
    rank_thes_eval == -1 ~ "low_extreme",
    rank_thes_eval == 1 ~ "high_extreme",
    rank_thes_eval == 0 ~ "mediocre",
    TRUE ~ "xxx")) %>% 
  arrange(dgrp) 

# Reformatting
extreme_dgrp_fraction_thes_summ <- extreme_dgrp_fraction_thes %>% 
  select(-rank_thes_eval, -n, -total_nr_phenotypes) %>% 
  pivot_wider(c(dgrp,sex),
              names_from = extreme_class,
              values_from = fraction_of_group) %>% 
  left_join(extreme_dgrp_fraction_thes %>% 
              distinct(dgrp, sex, total_nr_phenotypes),
            by = c("dgrp","sex")) %>% 
  mutate(fraction_extreme = rowSums(across(c(low_extreme, high_extreme)), na.rm = T))

save(extreme_dgrp_fraction_thes, file = "Extremes/extreme_dgrp_fraction_thes.RData")
save(extreme_dgrp_fraction_thes_summ, file = "Extremes/extreme_dgrp_fraction_thes_summ.RData")
write_tsv(extreme_dgrp_fraction_thes_summ, file = "Extremes/extreme_phenotype.tsv", col_names = T)



