# Files -----------
fig2D_phenotypes <- read_tsv("Figure 2D/listphenotypes.txt", col_names = 'phenotype_id')
phenotypes <- read_tsv("Extremes/phenotypes.tsv", col_names = T) %>% 
  mutate(phenotype_tag = paste("S",paste(study_id, phenotype_id, sep = "_"), sep =""))

load(file = "Extremes/extreme_dgrp_fraction_thes_summ.RData")
load(file = "Extremes/phenotype_pivot_curated_ranked.RData")

nr_phe <- 50

# Prep data ------
global_extreme_values_per_DGRP <- extreme_dgrp_fraction_thes_summ %>% 
  filter(total_nr_phenotypes >= nr_phe) %>% 
  filter(sex != "NA") %>% 
  select(-low_extreme, -mediocre, -high_extreme, -total_nr_phenotypes) %>% 
  rename(DGRP = "dgrp",
         FoE = "fraction_extreme")

subset_for_extremes_phenotype_values <- phenotype_pivot_curated_ranked %>% 
  separate(phenotype_tag, into = c("study_id","phenotype_id"), remove = F) %>% 
  filter(phenotype_id %in% fig2D_phenotypes$phenotype_id)

lifespan_pheno_info <- phenotypes %>% 
  filter(phenotype_id %in% fig2D_phenotypes$phenotype_id) %>% 
  filter(str_detect(tolower(phenotype_name), "longevity|lifespan")) %>% 
  distinct(study_id, phenotype_id, phenotype_tag, phenotype_name)

# Select top 5 FoE lines per sex
line_subset_extremes <- global_extreme_values_per_DGRP %>% 
  group_by(sex) %>% 
    arrange(desc(FoE)) %>% 
  slice(1:5) %>% 
  ungroup 

# Quick check number of phenotypes for lines (both sexes)
extreme_dgrp_fraction_thes_summ %>% 
  filter(dgrp %in% line_subset_extremes$DGRP) %>% 
  distinct(total_nr_phenotypes) 

# FoE for all lines
all_lines_extreme <- extreme_dgrp_fraction_thes_summ %>% 
  select(dgrp, sex, fraction_extreme) %>% 
  rename(FoE = "fraction_extreme")

# Generic mean across all lines
FoE_mean <- global_extreme_values_per_DGRP %>% 
  group_by(sex) %>% 
    summarise(mn_study = mean(FoE, na.rm = T)) %>% 
  ungroup() 

# Retrieve study means
study_means <- subset_for_extremes_phenotype_values %>% 
  filter(phenotype_tag %in% lifespan_pheno_info$phenotype_tag) %>% 
  group_by(phenotype_tag_sex) %>% 
  mutate(take_reps = case_when(
    str_detect(phenotype_tag_sex, "_F") ~ nber_sex_female,
    str_detect(phenotype_tag_sex, "_M") ~ nber_sex_male,
    str_detect(phenotype_tag_sex, "_NA") ~ nber_sex_na
  )) %>% 
  summarise(mn_study = mean(phenotype_value, na.rm = T),
            se_study = mean(phenotype_value, na.rm = T)/sqrt(take_reps)) %>% 
  distinct %>% 
  ungroup() 

# Plot ------------
ToPlot <- subset_for_extremes_phenotype_values %>% 
  filter(dgrp %in% c(line_subset_extremes$DGRP,"DGRP_852")) %>% 
  filter(phenotype_tag %in% lifespan_pheno_info$phenotype_tag) %>% 
  select(dgrp, phenotype_tag_sex, phenotype_value) %>% 
  bind_rows(study_means %>% 
              add_column(dgrp = "dgrp_mean") %>% 
              rename(phenotype_value = "mn_study")) %>% 
  mutate(line_indicator = case_when(
    str_detect(dgrp, "mean") ~ "DGRP Mean",
    dgrp == "DGRP_757" ~ "DGRP_757",
    dgrp == "DGRP_765" ~ "DGRP_765",
    dgrp == "DGRP_320" ~ "DGRP_320",
    dgrp == "DGRP_042" ~ "DGRP_042",
    dgrp == "DGRP_852" ~ "DGRP_852",
    TRUE ~ "additional_extreme_dgrp_lines"
  )) %>% 
  mutate(sex = case_when(
    str_detect(phenotype_tag_sex, "_F") ~ "F",
    str_detect(phenotype_tag_sex, "_M") ~ "M",
    str_detect(phenotype_tag_sex, "_NA") ~ "NA",
    TRUE ~ "xxx"
  )) %>% 
  mutate(dgrp_sex = paste(dgrp, sex, sep = "_")) %>% 
  left_join(all_lines_extreme,
            by = c("dgrp", "sex")) %>% 
  left_join(FoE_mean %>% 
              add_column(dgrp = "dgrp_mean") %>% 
              rename(FoE = "mn_study"),
            by = c("dgrp", "sex")) %>% 
  mutate(FoEc = coalesce(FoE.x, FoE.y)) %>% 
  select(-FoE.x, -FoE.y) %>% 
  rename(FoE = "FoEc")

ToPlot %>% 
  filter(dgrp == "dgrp_mean")

cols <- c("grey", "#000000","#800000", "#3cb44b","#000075", "#425e79", "#c86464")

lifespan_extremes_plot <- ggplot(ToPlot, aes(x = phenotype_tag_sex,
                                             y = phenotype_value,
                                             group = dgrp,
                                             fill = line_indicator)) + 
  geom_bar(stat = "identity",
           position = "dodge") +
  scale_fill_manual(values = cols) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
lifespan_extremes_plot

ggsave("Figure4E.pdf", lifespan_extremes_plot,  units = "cm", 
       width = 24, height = 10, device = "pdf")
ggsave("Figure4E.png", lifespan_extremes_plot,  units = "cm", 
       width = 24, height = 10, device = "png")



