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

# Load objects --------------------------------------------

load(file = "Extremes/phenotype_pivot_curated_ranked.RData")   
load(file = "Extremes/extreme_dgrp_fraction_thes.RData")
load(file = "Extremes/extreme_dgrp_fraction_thes_summ.RData")

# Set parameters ------------------------------------------ 

nr_phe <- 50

# Check correlation between Females and Males -------------

extreme_dgrp_fraction_thes_summ_tmp <- extreme_dgrp_fraction_thes_summ %>% 
  filter( total_nr_phenotypes >= nr_phe ) %>% 
  filter( sex != "NA") %>%   
  group_by(sex) %>% 
  arrange(fraction_extreme) %>% 
  ungroup

extreme_dgrp_fraction_sex_diff <- extreme_dgrp_fraction_thes_summ_tmp %>% 
  select(dgrp, sex, fraction_extreme) %>% 
  pivot_wider(dgrp,
              names_from = sex,
              values_from = fraction_extreme) %>% 
  rename(extreme_F = "F",
         extreme_M = "M")

c_pearson <- cor(extreme_dgrp_fraction_sex_diff$extreme_F, 
                 extreme_dgrp_fraction_sex_diff$extreme_M,
                 use = "pairwise.complete.obs",
                 method = "pearson")
c_spearman <- cor(extreme_dgrp_fraction_sex_diff$extreme_F, 
                  extreme_dgrp_fraction_sex_diff$extreme_M,
                  use = "pairwise.complete.obs",
                  method = "spearman")

smlm <- summary(lm(extreme_dgrp_fraction_sex_diff$extreme_F ~ extreme_dgrp_fraction_sex_diff$extreme_M))
rsq <- smlm$r.squared
pval <- smlm$coefficients[2,4]
message("Spearman correlation: ",round(c_spearman,3),"\n",
        "Pearson correlation: ",round(c_pearson,3),"\n",
        "Regression p-value: ",pval,", with R-squared: ",round(rsq,3))

ggplot(extreme_dgrp_fraction_sex_diff, 
       aes(x = extreme_F, y = extreme_M,
           label = dgrp)) +
  geom_point() +
  scale_x_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1,.1)) +
  scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1,.1)) +
  geom_text(hjust = 0)


# Individual lines and extremes ---------------------------

extreme_lines <- extreme_dgrp_fraction_thes_summ_tmp %>% 
  group_by(sex) %>% 
    filter(fraction_extreme == min(fraction_extreme) | 
             fraction_extreme == max(fraction_extreme)) %>% 
  ungroup() %>% 
  arrange(sex) %>% 
  mutate(dgrp_sex_tag = paste(dgrp, sex, sep = "_"))

extreme_lines_perPheno_fraction <- phenotype_pivot_curated_ranked %>% 
  mutate(dgrp_sex_tag = paste(dgrp, sex, sep = "_")) %>% 
  filter(dgrp_sex_tag %in% extreme_lines$dgrp_sex_tag) %>% 
  filter(sex != "NA") %>% 
  select(dgrp, dgrp_sex_tag, phenotype_tag_sex, phenotype_tag, sex, 
         nber_sex_female, nber_sex_male,
         rank_real, rank_thes_base_cutoff, rank_thes_prop, rank_thes_eval) %>% 
  mutate(rank_thes_prop_adj = case_when(
    rank_thes_prop > 0.5 ~ 1 - rank_thes_prop,
    TRUE ~ rank_thes_prop))
 
save(extreme_lines_perPheno_fraction, file = "Extremes/extreme_lines_perPheno_fraction.RData")

# Lines where the sex are within 0.05 ---------------------
# Lowest and highest mean, while difference is < 0.05

extreme_lines_both_sex <- extreme_dgrp_fraction_sex_diff %>% 
  mutate(ext_diff = abs(extreme_F-extreme_M)) %>% 
  rowwise() %>% 
    mutate(ext_mn = mean(c(extreme_F,extreme_M), na.rm = T)) %>% 
  ungroup() %>% 
  filter(is.na(extreme_F)==F & is.na(extreme_M)==F) %>% 
  arrange(ext_mn) %>% 
  filter(ext_diff <= 0.05) %>% 
  slice(1:2, nrow(.):(nrow(.)-1)) 

extreme_lines_both_sex_perPheno_fraction <- phenotype_pivot_curated_ranked %>% 
  mutate(dgrp_sex_tag = paste(dgrp, sex, sep = "_")) %>% 
  filter(dgrp %in% extreme_lines_both_sex$dgrp) %>% 
  filter(sex != "NA") %>% 
  select(dgrp, dgrp_sex_tag, phenotype_tag_sex, phenotype_tag, sex, 
         nber_sex_female, nber_sex_male,
         rank_real, rank_thes_base_cutoff, rank_thes_prop, rank_thes_eval) %>% 
  mutate(rank_thes_prop_adj = case_when(
    rank_thes_prop > 0.5 ~ 1 - rank_thes_prop,
    TRUE ~ rank_thes_prop))

save(extreme_lines_both_sex_perPheno_fraction, file = "Extremes/extreme_lines_both_sex_perPheno_fraction.RData")













