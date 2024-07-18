# Set environment ----------------------------------------- R version 4.0.3 (2020-10-10)

library(biomaRt)
library(tidyverse)#2.0.0
library(data.table)#1.14.2
library(tidylog)#1.0.2
library(jsonlite)#1.8.0
library(ggplot2)#3.3.5
library(patchwork)#1.1.2

setwd("/Users/roelbevers/OneDrive - Genomics England Ltd/Documents/Personal/EPFL/dgrpool/")

# Load objects --------------------------------------------

allhits_foecov_fem <- read_tsv("GWAS_longevity_extremes/longevity_f_foe_cov/Arya_Longevity_F_F.glm.linear")
tophits_foecov_fem <- read_tsv("GWAS_longevity_extremes/longevity_f_foe_cov/Arya_Longevity_F_F.glm.linear.top_0.01.annot.tsv")

allhits_stdcov_fem <- read_tsv("GWAS_longevity_extremes/longevity_f_std_cov/Arya_Longevity_F_F.glm.linear")
tophits_stdcov_fem <- read_tsv("GWAS_longevity_extremes/longevity_f_std_cov/Arya_Longevity_F_F.glm.linear.top_0.01.annot.tsv")

allhits_foecov_mal <- read_tsv("GWAS_longevity_extremes/longevity_m_foe_cov/Arya_Longevity_M_M.glm.linear")
tophits_foecov_mal <- read_tsv("GWAS_longevity_extremes/longevity_m_foe_cov/Arya_Longevity_M_M.glm.linear.top_0.01.annot.tsv")

allhits_stdcov_mal <- read_tsv("GWAS_longevity_extremes/longevity_m_std_cov/Arya_Longevity_M_M.glm.linear")
tophits_stdcov_mal <- read_tsv("GWAS_longevity_extremes/longevity_m_std_cov/Arya_Longevity_M_M.glm.linear.top_0.01.annot.tsv")

# All annotations -----------------------------------------

all_annotations <- bind_rows(tophits_foecov_fem %>% distinct(ID, gene_annotation, regulatory_annotation),
                             tophits_stdcov_fem %>% distinct(ID, gene_annotation, regulatory_annotation),
                             tophits_foecov_mal %>% distinct(ID, gene_annotation, regulatory_annotation),
                             tophits_stdcov_mal %>% distinct(ID, gene_annotation, regulatory_annotation)
) %>% distinct

# Gene lists for conversions ------------------------------

listMarts(host="useast.ensembl.org")
ensembl_us_east = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="useast.ensembl.org")
listDatasets(ensembl_us_east) %>% 
  as_tibble() %>% 
  filter(str_detect(dataset, "dmel"))
ensembl = useEnsembl(biomart="ensembl", dataset="dmelanogaster_gene_ensembl")
listAttributes(ensembl) %>% as_tibble() %>% slice(1:20)

dmel_chroms <- c('2R','2L','3R','3L','X','Y','4','mitochondrion_genome')

genes_and_ids <- NULL
for( i in 1:length(dmel_chroms)) { 
  genes_and_ids <- getBM(attributes = 
                           c('flybase_gene_id',
                             'external_gene_name',
                             'chromosome_name',
                             'flybase_transcript_id',
                             'flybase_annotation_id'
                           ),
                         filters = 'chromosome_name',
                         values = dmel_chroms[i],
                         mart = ensembl) %>% 
    as_tibble() %>% 
    mutate(chromosome_name = as.character(chromosome_name)) %>% 
    bind_rows(genes_and_ids) %>% 
    distinct
}
genes_and_ids

genes_and_ids %>% 
  count(chromosome_name, sort = T)

#save(genes_and_ids, file = "Extremes/genes_and_ids_25042023.RData")

genes_and_ids_clean <- genes_and_ids %>% 
  distinct(flybase_gene_id, flybase_annotation_id, external_gene_name) %>% 
  pivot_longer(-flybase_gene_id,
               names_to = "type",
               values_to = "id") %>% 
  group_by(flybase_gene_id) %>% 
  distinct %>% 
  filter(str_detect(id, "^CG|^CR") | type == "external_gene_name") %>% 
  ungroup() %>% 
  pivot_wider(flybase_gene_id,
              names_from = type,
              values_from = id,
              values_fill = NA) %>% 
  left_join(genes_and_ids %>% select(-flybase_transcript_id, -flybase_annotation_id), 
            by = c("flybase_gene_id","external_gene_name")) %>% 
  distinct

genes_and_ids_clean %>%   
  filter(flybase_gene_id == "FBgn0001297")


# Gene list "determination of adult lifespan -------------

listMarts(host="useast.ensembl.org")
ensembl_us_east = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="useast.ensembl.org")
listDatasets(ensembl_us_east) %>% 
  as_tibble() %>% 
  filter(str_detect(dataset, "dmel"))
ensembl = useEnsembl(biomart="ensembl", dataset="dmelanogaster_gene_ensembl")

lifespan_genes <- getBM(attributes = 
                          c('flybase_gene_id',
                            'go_id',
                            'ensembl_transcript_id'),
                        filters = 'go',
                        values = 'GO:0008340',
                        mart = ensembl) %>% 
  as_tibble() %>% 
  filter(go_id == 'GO:0008340')

lifespan_genes_annot <- lifespan_genes %>% 
  left_join(genes_and_ids_clean,
            by = "flybase_gene_id") %>% 
  add_column(lifespan_gene = "PASS")

lifespan_genes_annot %>% distinct(flybase_gene_id)
lifespan_genes_annot %>% 
  write_tsv("lifespan_genes_GO0008340.tsv", col_names = T)

lifespan_genes_annot %>% 
  distinct(flybase_gene_id) %>% 
  mutate(Xtype = paste("X",1:nrow(.),sep= "")) %>% 
  pivot_wider(names_from = Xtype,
              values_from = flybase_gene_id) %>% 
  add_column(ID = "lifespan_GO0008340",
             source = "NA") %>% 
  #unite("gene_list", contains("X"), sep = "\t") %>% 
  select(ID, source, contains("X")) %>% 
  write_tsv("lifespan_genes_GO0008340.gmt", col_names = F)

# Evaluate -----------------------------------------------

tophits_foecov_fem # 17,550 hits
tophits_stdcov_fem # 17,861 hits # Standard model has more hits

tophits_foecov_mal # 22,804 hits
tophits_stdcov_mal # 23,106 hits # Standard model has more hits

# Number of significant hits ----------------------------

p_eval_stdcovmal <- tophits_stdcov_mal %>% 
  mutate(p_category = case_when(
    P <= 0.0000000001 ~ "<= 0.0000000001",
    P <= 0.000000001 ~ "<= 0.000000001",
    P <= 0.00000001 ~ "<= 0.00000001",
    P <= 0.0000001 ~ "<= 0.0000001",
    P <= 0.000001 ~ "<= 0.000001",
    P <= 0.00001 ~ "<= 0.00001",
    P <= 0.0001 ~ "<= 0.0001",
    P <= 0.001 ~ "<= 0.001",
    P <= 0.01 ~ "<= 0.01",
    TRUE ~ "xxx"
  ))

p_eval_foecovmal <- tophits_foecov_mal %>% 
  mutate(p_category = case_when(
    P <= 0.0000000001 ~ "<= 0.0000000001",
    P <= 0.000000001 ~ "<= 0.000000001",
    P <= 0.00000001 ~ "<= 0.00000001",
    P <= 0.0000001 ~ "<= 0.0000001",
    P <= 0.000001 ~ "<= 0.000001",
    P <= 0.00001 ~ "<= 0.00001",
    P <= 0.0001 ~ "<= 0.0001",
    P <= 0.001 ~ "<= 0.001",
    P <= 0.01 ~ "<= 0.01",
    TRUE ~ "xxx"
  ))

p_eval_stdcovmal %>% 
  count(p_category) %>% 
  add_column(model = "standard model") %>% 
  bind_rows(p_eval_foecovmal %>% 
              count(p_category) %>% 
              add_column(model = "FoE model")
  ) # relatively equal

p_eval_stdcovfem <- tophits_stdcov_fem %>% 
  mutate(p_category = case_when(
    P <= 0.0000000001 ~ "<= 0.0000000001",
    P <= 0.000000001 ~ "<= 0.000000001",
    P <= 0.00000001 ~ "<= 0.00000001",
    P <= 0.0000001 ~ "<= 0.0000001",
    P <= 0.000001 ~ "<= 0.000001",
    P <= 0.00001 ~ "<= 0.00001",
    P <= 0.0001 ~ "<= 0.0001",
    P <= 0.001 ~ "<= 0.001",
    P <= 0.01 ~ "<= 0.01",
    TRUE ~ "xxx"
  ))

p_eval_foecovfem <- tophits_foecov_fem %>% 
  mutate(p_category = case_when(
    P <= 0.0000000001 ~ "<= 0.0000000001",
    P <= 0.000000001 ~ "<= 0.000000001",
    P <= 0.00000001 ~ "<= 0.00000001",
    P <= 0.0000001 ~ "<= 0.0000001",
    P <= 0.000001 ~ "<= 0.000001",
    P <= 0.00001 ~ "<= 0.00001",
    P <= 0.0001 ~ "<= 0.0001",
    P <= 0.001 ~ "<= 0.001",
    P <= 0.01 ~ "<= 0.01",
    TRUE ~ "xxx"
  ))

p_eval_stdcovfem %>% 
  count(p_category) %>% 
  add_column(model = "standard model") %>% 
  bind_rows(p_eval_foecovfem %>% 
              count(p_category) %>% 
              add_column(model = "FoE model")
  ) # relatively equal

# Overlapping hits --------------------------------------

allhits_foecov_fem %>% filter(P < 0.01)
tophits_foecov_fem
allhits_stdcov_fem %>% filter(P < 0.01)
tophits_stdcov_fem

allhits_foecov_mal %>% filter(P < 0.01)
tophits_foecov_mal
allhits_stdcov_mal %>% filter(P < 0.01)
tophits_stdcov_mal

allhits_fem <- allhits_stdcov_fem %>%
  left_join(allhits_foecov_fem %>% 
              select(-`#CHROM`,-POS,-REF,-ALT,-A1,-TEST),
            by = "ID",
            suffix = c("_std","_foe")) %>% 
  mutate(pval_compare = case_when(
    is.na(P_std) & is.na(P_foe)==F ~ "absent_P_std",
    is.na(P_foe) & is.na(P_std)==F ~ "absent_P_foe",
    is.na(P_std) & is.na(P_foe) ~ "absent_Pvals",
    (P_std <= 0.01) == F & (P_foe <= 0.01) == F ~ "not_top_in_either",
    (P_std <= 0.01) == T & (P_foe <= 0.01) == T ~ "top_in_both",
    (P_std <= 0.01) == T & (P_foe <= 0.01) == F ~ "top_in_std",
    (P_std <= 0.01) == F & (P_foe <= 0.01) == T ~ "top_in_foe",
    TRUE ~ "xxx"
  )) 
allhits_fem %>% 
  count(pval_compare) # 78.3% of the top hits for FoE-cov is overlapping. 
                      # 76.9% of overlap for the standard-cov.

allhits_mal <- allhits_stdcov_mal %>%
  left_join(allhits_foecov_mal %>% 
              select(-`#CHROM`,-POS,-REF,-ALT,-A1,-TEST),
            by = "ID",
            suffix = c("_std","_foe")) %>% 
  mutate(pval_compare = case_when(
    is.na(P_std) & is.na(P_foe)==F ~ "absent_P_std",
    is.na(P_foe) & is.na(P_std)==F ~ "absent_P_foe",
    is.na(P_std) & is.na(P_foe) ~ "absent_Pvals",
    (P_std <= 0.01) == F & (P_foe <= 0.01) == F ~ "not_top_in_either",
    (P_std <= 0.01) == T & (P_foe <= 0.01) == T ~ "top_in_both",
    (P_std <= 0.01) == T & (P_foe <= 0.01) == F ~ "top_in_std",
    (P_std <= 0.01) == F & (P_foe <= 0.01) == T ~ "top_in_foe",
    TRUE ~ "xxx"
  )) 
allhits_mal %>% 
  count(pval_compare) # 55.8% of the top hits for FoE-cov is overlapping. 
                      # 55.1% of overlap for the standard-cov.

# Top-5 hits -------------------------------------- [ FEM ]

top5hits_fem_anno <- allhits_fem %>% 
  arrange(P_foe) %>% 
  add_column(model = "FoE_model") %>% 
  slice(1:10) %>% 
  bind_rows(allhits_fem %>% 
              arrange(P_std) %>% 
              slice(1:10) %>%
              add_column(model = "Std_model")
  ) %>% 
  left_join(all_annotations, 
            by = "ID")
top5hits_fem_anno %>% View

# Top-5 hits -------------------------------------- 
# What is the most sig in females, that is unique to FoE
allhits_fem %>% 
  filter(pval_compare == "top_in_foe") %>% 
  left_join(all_annotations, 
            by = "ID") %>% 
  arrange(P_foe) %>% View # 0.000846448

allhits_mal %>% 
  filter(pval_compare == "top_in_foe") %>% 
  left_join(all_annotations, 
            by = "ID") %>% 
  arrange(P_foe) %>% View # 5.88884e-05

# Top-100 hits --------------------------------------

allhits_fem_rankeval <- allhits_fem %>% 
  left_join(all_annotations, 
            by = "ID") %>% 
  arrange(P_std) %>% 
  mutate(P_std_rank = 1:nrow(.)) %>% 
  arrange(P_foe) %>% 
  mutate(P_foe_rank = 1:nrow(.)) %>% 
  mutate(rank_delta = P_std_rank-P_foe_rank) %>% 
  mutate(rank_eval = case_when(
    P_std_rank > 100 & P_foe_rank > 100 ~ "non_top100",
    P_std_rank <= 100 & P_foe_rank <= 100 ~ "Shared_top100",
    P_std_rank <= 100 & P_foe_rank > 100 ~ "Std_top100",
    P_std_rank > 100 & P_foe_rank <= 100 ~ "FoE_top100"
  ))

allhits_fem_rankeval %>% 
  count(rank_eval)

allhits_fem_rankeval %>% 
  filter(rank_eval %in% c("Std_top100","FoE_top100")) %>% 
  summarise(mn = mean(abs(rank_delta)),
            mdn = median(abs(rank_delta)),
            min = min(abs(rank_delta)),
            max = max(abs(rank_delta)),
            n = nrow(.))

allhits_mal_rankeval <- allhits_mal %>% 
  left_join(all_annotations, 
            by = "ID") %>% 
  arrange(P_std) %>% 
  mutate(P_std_rank = 1:nrow(.)) %>% 
  arrange(P_foe) %>% 
  mutate(P_foe_rank = 1:nrow(.)) %>% 
  mutate(rank_delta = P_std_rank-P_foe_rank) %>% 
  mutate(rank_eval = case_when(
    P_std_rank > 100 & P_foe_rank > 100 ~ "non_top100",
    P_std_rank <= 100 & P_foe_rank <= 100 ~ "Shared_top100",
    P_std_rank <= 100 & P_foe_rank > 100 ~ "Std_top100",
    P_std_rank > 100 & P_foe_rank <= 100 ~ "FoE_top100"
  ))

allhits_mal_rankeval %>% 
  count(rank_eval)

allhits_mal_rankeval %>% 
  filter(rank_eval %in% c("Std_top100","FoE_top100")) %>% 
  summarise(mn = mean(abs(rank_delta)),
            mdn = median(abs(rank_delta)),
            min = min(abs(rank_delta)),
            max = max(abs(rank_delta)),
            n = nrow(.))

# Parameters of FoE Cov input -----------------------------

load("Extremes/extreme_dgrp_fraction_thes_summ.RData")

dgrp_foe <- extreme_dgrp_fraction_thes_summ %>% 
  filter(total_nr_phenotypes >= 50) %>% 
  mutate(covID = str_replace(dgrp,"DGRP","line")) %>% 
  mutate(covID = str_replace(covID,"_00","_")) %>% 
  mutate(covID = str_replace(covID,"_0","_")) %>% 
  rename(foe = "fraction_extreme")

dgrp_foe %>% 
  count(sex)

dgrp_foe %>% 
  filter(sex != "NA") %>% 
  group_by(sex) %>% 
  summarise(total = n_distinct(dgrp),
            med_nr_pheno = median(total_nr_phenotypes),
            avg_nr_pheno = mean(total_nr_phenotypes),
            mn_nr_pheno = min(total_nr_phenotypes),
            mx_nr_pheno = max(total_nr_phenotypes))

# Correlation of Beta's --------------------------------------

cor_stats <- function(x,y) {
  
  c_pearson <- cor(x, y,
                   use = "pairwise.complete.obs",
                   method = "pearson")
  
  c_spearman <- cor(x, y,
                    use = "pairwise.complete.obs",
                    method = "spearman")
  
  smlm <- summary(lm(x ~ y))
  rsq <- smlm$r.squared
  pval <- smlm$coefficients[2,4]
  message("Spearman correlation: ",round(c_spearman,4),"\n",
          "Pearson correlation: ",round(c_pearson,4),"\n",
          "Regression p-value: ",pval,", with R-squared: ",round(rsq,3))
  
}

# Across all hits
cor_stats(allhits_fem$P_std, allhits_fem$P_foe) # P-values Females
cor_stats(allhits_fem$BETA_std, allhits_fem$BETA_foe) # Beta's Females

cor_stats(allhits_mal$P_std, allhits_mal$P_foe) # P-value Males
cor_stats(allhits_mal$BETA_std, allhits_mal$BETA_foe) # Beta's Males

# Across 'top' hits
cor_stats(allhits_fem %>% 
            filter(pval_compare %in% c("top_in_both","top_in_foe","top_in_std")) %>% 
            pull(P_std), 
          allhits_fem %>% 
            filter(pval_compare %in% c("top_in_both","top_in_foe","top_in_std")) %>% 
            pull(P_foe)) # P-values Females
cor_stats(allhits_fem %>% 
            filter(pval_compare %in% c("top_in_both","top_in_foe","top_in_std")) %>% 
            pull(BETA_std), 
          allhits_fem %>% 
            filter(pval_compare %in% c("top_in_both","top_in_foe","top_in_std")) %>% 
            pull(BETA_foe)) # Beta's Females

cor_stats(allhits_mal %>% 
            filter(pval_compare %in% c("top_in_both","top_in_foe","top_in_std")) %>% 
            pull(P_std), 
          allhits_mal %>% 
            filter(pval_compare %in% c("top_in_both","top_in_foe","top_in_std")) %>% 
            pull(P_foe)) # P-values Males
cor_stats(allhits_mal %>% 
            filter(pval_compare %in% c("top_in_both","top_in_foe","top_in_std")) %>% 
            pull(BETA_std), 
          allhits_mal %>% 
            filter(pval_compare %in% c("top_in_both","top_in_foe","top_in_std")) %>% 
            pull(BETA_foe)) # Beta's Males

# Extract genes from annotation file ----------------------

regex_filt <- function(x, regexStr) {
  t1 <- str_extract_all(x,regexStr) %>% unlist %>% unique
  t2 <- paste0(t1, collapse = ",", sep = "")
  return(t2)
}

annotations_string_captured <- all_annotations %>% 
  separate(gene_annotation, 
           into = paste("X",1:50, sep=""),
           remove = F,
           sep = "\\|") %>% 
  mutate(CG_capture = unlist(map2(gene_annotation, "CG[0-9]{4,5}", regex_filt))) %>% 
  mutate(CR_capture = unlist(map2(gene_annotation, "CR[0-9]{4,5}", regex_filt))) %>% 
  mutate(FBtr_capture = unlist(map2(gene_annotation, "FBtr[0-9]{7}", regex_filt))) %>% 
  mutate(FBgn_capture = unlist(map2(gene_annotation, "FBgn[0-9]{7}", regex_filt)))

annotations_string_captured %>% 
  select(gene_annotation, X2, CG_capture, CR_capture, FBtr_capture, FBgn_capture) %>% 
  View

# Clean-up output -----------------------------------------
annotations_string_fbgn_only <- annotations_string_captured %>% 
  select(-contains("X")) %>% 
  select(ID, gene_annotation, FBgn_capture) %>% 
  filter(FBgn_capture != "") %>% 
  separate(FBgn_capture,
           into = paste("X",1:12, sep=""),
           remove = T,
           sep = ",") %>% 
  pivot_longer(-c(ID, gene_annotation),
               names_to = "type",
               values_to = "FBgn_ID") %>% 
  select(-type) %>% 
  filter(is.na(FBgn_ID)==F)

allhits_fem_ranked_with_clean_genes <- allhits_fem %>% 
  arrange(P_std) %>% 
  mutate(P_std_rank = 1:nrow(.)) %>% 
  arrange(P_foe) %>% 
  mutate(P_foe_rank = 1:nrow(.)) %>% 
  left_join(annotations_string_fbgn_only,
            by = "ID") %>% 
  left_join(genes_and_ids_clean %>% 
              distinct(flybase_gene_id, flybase_annotation_id, external_gene_name),
            by = c("FBgn_ID" = "flybase_gene_id")) # Inflates the rows a bit (multiple genes for 1 variant)

allhits_mal_ranked_with_clean_genes <- allhits_mal %>% 
  arrange(P_std) %>% 
  mutate(P_std_rank = 1:nrow(.)) %>% 
  arrange(P_foe) %>% 
  mutate(P_foe_rank = 1:nrow(.)) %>% 
  left_join(annotations_string_fbgn_only,
            by = "ID") %>% 
  left_join(genes_and_ids_clean %>% 
              distinct(flybase_gene_id, flybase_annotation_id, external_gene_name),
            by = c("FBgn_ID" = "flybase_gene_id")) # Inflates the rows a bit (multiple genes for 1 variant)


# Plot correlation of Beta's and P-values ----------------- [ Quick plots ]

ToPlot <- allhits_fem %>% 
  filter(P_std <= 0.01 | P_foe <= 0.01)

ggplot(ToPlot, aes(x = BETA_std,
                   y = BETA_foe)) +
  geom_point()

ToPlot <- allhits_mal %>% 
  filter(P_std <= 0.01 | P_foe <= 0.01)

ggplot(ToPlot, aes(x = BETA_std,
                   y = BETA_foe)) +
  geom_point()

ToPlot <- allhits_fem %>% 
  filter(P_std <= 0.01 | P_foe <= 0.01)

ggplot(ToPlot, aes(x = -log10(P_std),
                   y = -log10(P_foe))) +
  geom_point()

ToPlot <- allhits_mal %>% 
  filter(P_std <= 0.01 | P_foe <= 0.01)

ggplot(ToPlot, aes(x = -log10(P_std),
                   y = -log10(P_foe))) +
  geom_point()


# Lifespan gene focussed ---------------------------------- [ Includes Correlation plots]

# Mark all lifespan genes
allhits_mal_ranked_with_clean_genes_lfsp_marked <- allhits_mal_ranked_with_clean_genes %>% 
  select(-flybase_annotation_id) %>% 
  left_join(lifespan_genes_annot,
            by = c("FBgn_ID" = "flybase_gene_id",
                   "external_gene_name"))

allhits_fem_ranked_with_clean_genes_lfsp_marked <- allhits_fem_ranked_with_clean_genes %>% 
  select(-flybase_annotation_id) %>% 
  left_join(lifespan_genes_annot,
            by = c("FBgn_ID" = "flybase_gene_id",
                   "external_gene_name"))
# Prepare Plots
# Males
allhits_mal_ranked_with_clean_genes_lfsp_marked %>% 
  filter(lifespan_gene == "PASS") %>% 
  select(ID, BETA_std, BETA_foe, P_std, P_foe, pval_compare) %>% 
  distinct %>% 
  mutate(pval_diff = P_std-P_foe) %>% 
  mutate(sign_pval_diff = sign(pval_diff)) %>% 
  count(pval_compare, sign_pval_diff)

allhits_mal_ranked_with_clean_genes_lfsp_marked %>% 
  select(ID, BETA_std, BETA_foe, P_std, P_foe, pval_compare, lifespan_gene) %>% 
  distinct %>% 
  arrange(P_foe) %>% 
  slice(1:100) %>% 
  count(lifespan_gene)

allhits_mal_ranked_with_clean_genes_lfsp_marked %>% 
  select(ID, BETA_std, BETA_foe, P_std, P_foe, pval_compare, lifespan_gene) %>% 
  distinct %>% 
  arrange(P_std) %>% 
  slice(1:100) %>% 
  count(lifespan_gene)

# Females
allhits_fem_ranked_with_clean_genes_lfsp_marked %>% 
  filter(lifespan_gene == "PASS") %>% 
  select(ID, BETA_std, BETA_foe, P_std, P_foe, pval_compare) %>% 
  distinct %>% 
  mutate(pval_diff = P_std-P_foe) %>% 
  mutate(sign_pval_diff = sign(pval_diff)) %>% 
  count(pval_compare, sign_pval_diff) 

allhits_fem_ranked_with_clean_genes_lfsp_marked %>% 
  select(ID, BETA_std, BETA_foe, P_std, P_foe, pval_compare, lifespan_gene) %>% 
  distinct %>% 
  arrange(P_foe) %>% 
  slice(1:100) %>% 
  count(lifespan_gene)

allhits_fem_ranked_with_clean_genes_lfsp_marked %>% 
  select(ID, BETA_std, BETA_foe, P_std, P_foe, pval_compare, lifespan_gene) %>% 
  distinct %>% 
  arrange(P_std) %>% 
  slice(1:100) %>% 
  count(lifespan_gene)

# Plotting
# Fem BETA
ToPlot <- allhits_fem_ranked_with_clean_genes_lfsp_marked %>% 
  filter(P_std <= 0.01 | P_foe <= 0.01) %>% 
  arrange(lifespan_gene)

ggplot(ToPlot) +
  geom_point(aes(x = BETA_std,
                 y = BETA_foe,
                 colour = lifespan_gene)) +
  geom_point(data = subset(ToPlot, lifespan_gene == "PASS"),
             aes(x = BETA_std,
                 y = BETA_foe,
                 colour = lifespan_gene)) + 
  scale_colour_manual(values = c("red","black"))

# Fem PVAL
ggplot(ToPlot) +
  geom_point(aes(x = -log10(P_std),
                 y = -log10(P_foe),
                 colour = lifespan_gene)) +
  geom_point(data = subset(ToPlot, lifespan_gene == "PASS"),
             aes(x = -log10(P_std),
                 y = -log10(P_foe),
                 colour = lifespan_gene)) + 
  scale_colour_manual(values = c("red","black"))

# Mal BETA
ToPlot <- allhits_mal_ranked_with_clean_genes_lfsp_marked %>% 
  filter(P_std <= 0.01 | P_foe <= 0.01) %>% 
  arrange(lifespan_gene)

ggplot(ToPlot) +
  geom_point(aes(x = BETA_std,
                 y = BETA_foe,
                 colour = lifespan_gene)) +
  geom_point(data = subset(ToPlot, lifespan_gene == "PASS"),
             aes(x = BETA_std,
                 y = BETA_foe,
                 colour = lifespan_gene)) + 
  scale_colour_manual(values = c("red","black"))

# Mal PVAL
ggplot(ToPlot) +
  geom_point(aes(x = -log10(P_std),
                 y = -log10(P_foe),
                 colour = lifespan_gene)) +
  geom_point(data = subset(ToPlot, lifespan_gene == "PASS"),
             aes(x = -log10(P_std),
                 y = -log10(P_foe),
                 colour = lifespan_gene)) + 
  scale_colour_manual(values = c("red","black"))





