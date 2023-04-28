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
library(wordcloud)#2.6
library(gridExtra)

# Load objects --------------------------------------------

load(file = "Extremes/extreme_dgrp_fraction_thes.RData")
load(file = "Extremes/extreme_lines_perPheno_fraction.RData")
load(file = "Extremes/extreme_lines_both_sex_perPheno_fraction.RData")

# Set parameters ------------------------------------------ 

nr_phe <- 50
col_fem <- "#F28C28"
col_mal <- "#0E4C92"


# Plotting distributions ---------------------------------- [ Females, 50 phenotypes minimal ]

input_object <- extreme_dgrp_fraction_thes

dgrp_order <- input_object %>% 
  filter(extreme_class == "mediocre") %>% 
  filter(sex == "F") %>% 
  filter(total_nr_phenotypes >= nr_phe) %>% 
  arrange(desc(fraction_of_group)) %>% 
  pull(dgrp)

ToPlot <- input_object %>% 
  filter(sex == "F") %>% 
  filter(total_nr_phenotypes >= nr_phe)

ToPlot$dgrp <- factor(ToPlot$dgrp, levels = unique(dgrp_order))

ggplot(ToPlot, aes(x = dgrp,
                   y = fraction_of_group,
                   fill = extreme_class)) + 
  geom_bar(stat = "identity") + #, fill = "gold", colour = "black") +
  scale_fill_manual(values = c("gold","darkgreen","grey")) +
  theme_bw() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(size=10, color="black"),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1,.1))

# Plot by sex ---------------------------------------------

input_object <- extreme_dgrp_fraction_thes

ToPlotF <- input_object %>% 
  filter(extreme_class == "mediocre") %>% 
  filter(sex == "F") %>% 
  filter(total_nr_phenotypes >= nr_phe) %>% 
  arrange(fraction_of_group)
ToPlotF$dgrp <- factor(ToPlotF$dgrp, levels = ToPlotF$dgrp)

ToPlotM <- input_object %>% 
  filter(extreme_class == "mediocre") %>% 
  filter(sex == "M") %>% 
  filter(total_nr_phenotypes >= nr_phe) %>% 
  arrange(fraction_of_group)
ToPlotM$dgrp <- factor(ToPlotM$dgrp, levels = ToPlotM$dgrp)

pF <- ggplot(ToPlotF, 
             aes(x = dgrp,
                 y = 1-fraction_of_group)) + 
  geom_bar(stat = "identity", fill = col_fem) +
  labs(title = "Female", x = "", y = "Fraction of extremeness", fill = "") +
  theme_bw() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(size=10, color="black"),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1,.1))

pM <- ggplot(ToPlotM, 
             aes(x = dgrp,
                 y = 1-fraction_of_group)) + 
  geom_bar(stat = "identity", fill = col_mal) +
  labs(title = "Male", x = "", y = "", fill = "") +
  theme_bw() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1,.1))

pComb <- grid.arrange(patchworkGrob(( pF + pM )), bottom = "DGRP lines")

ggsave("Extremes/4A_extremes.pdf", pComb,  units = "cm", 
       width = 20, height = 10, device = "pdf")
ggsave("Extremes/4A_extremes.png", pComb,  units = "cm", 
       width = 20, height = 10, device = "png")

# Scatter plot extremes ----------------------------------- 

load(file = "Extremes/extreme_dgrp_fraction_sex_diff.RData")

input_object <- extreme_dgrp_fraction_sex_diff %>% 
  filter(is.na(extreme_F)==F & is.na(extreme_M)==F )

c_pearson <- cor(input_object$extreme_F, 
                 input_object$extreme_M,
                 use = "pairwise.complete.obs",
                 method = "pearson")
c_spearman <- cor(input_object$extreme_F, 
                  input_object$extreme_M,
                  use = "pairwise.complete.obs",
                  method = "spearman")

smlm <- summary(lm(input_object$extreme_F ~ input_object$extreme_M))
rsq <- smlm$r.squared
pval <- smlm$coefficients[2,4]
message("Spearman correlation: ",round(c_spearman,3),"\n",
        "Pearson correlation: ",round(c_pearson,3),"\n",
        "Regression p-value: ",pval,", with R-squared: ",round(rsq,3))

pExtScatt <- ggplot(input_object, 
                    aes(x = extreme_F, 
                        y = extreme_M)) +
  geom_point(colour = "#5a0e2d") +
  labs(title = "Correlation between males and females", 
       x = "Fraction of extremeness (Females)", 
       y = "Fraction of extremeness (Males)", 
       fill = "") +
  scale_x_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1,.1)) +
  scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1,.1)) + 
  theme_bw() +
  geom_abline(intercept = smlm$coefficients[1],
              slope = smlm$coefficients[2],
              colour = "red") + 
  ggplot2::annotate("text", x = 0.1, y = 0.9, hjust = 0,
                    label = paste0("Pearson's r = ",round(c_pearson,4),"\n",
                                   "Spearman's ρ= ",round(c_spearman,4)))

pExtScatt 

ggsave("Extremes/4B_scatter_extremes.pdf", pExtScatt,  units = "cm", 
       width = 20, height = 20, device = "pdf")
ggsave("Extremes/4B_scatter_extremes.png", pExtScatt,  units = "cm", 
       width = 20, height = 20, device = "png")


# Individual extremes -------------------------------------

ToPlot <- extreme_lines_perPheno_fraction 
ToPlot$dgrp_sex_tag <- factor(ToPlot$dgrp_sex_tag, 
                              levels = c("DGRP_879_F","DGRP_757_F","DGRP_783_M","DGRP_352_M"))
ToPlot$dgrp <- factor(ToPlot$dgrp, 
                              levels = c("DGRP_879","DGRP_757","DGRP_783","DGRP_352"))

pExtLines <- ggplot(ToPlot, aes(x = dgrp,
                   y = rank_thes_prop_adj,
                   fill = sex)) + 
  geom_boxplot() +
  scale_fill_manual(values = c(col_fem, col_mal)) + 
  geom_jitter(position=position_jitter(0.2)) +
  labs(title = "Most extreme and medicore lines", 
       x = "", 
       y = "Adjusted fraction of extremeness (1-x if x > 0.5)", 
       fill = "") +
  theme_bw() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.text.y = element_text(size=10, color="black"),
        axis.text.x = element_text(size=10, color="black"),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  scale_y_continuous(expand=c(0,0), limits = c(-.01,.51), breaks = seq(0,1,.1)) 
pExtLines

ggsave("Extremes/4C_most_extreme_lines.pdf", pExtLines,  units = "cm", 
       width = 20, height = 20, device = "pdf")
ggsave("Extremes/4C_most_extreme_lines.png", pExtLines,  units = "cm", 
       width = 20, height = 20, device = "png")

# Individual extremes ------------------------------------- [ True extremes (for both sexes) ]

ToPlot <- extreme_lines_both_sex_perPheno_fraction
#ToPlot$dgrp_sex_tag <- factor(ToPlot$dgrp_sex_tag, 
#                      levels = c("DGRP_383_F","DGRP_280_F",
#                                 "DGRP_852_F","DGRP_042_F",
#                                 "DGRP_383_M","DGRP_280_M",
#                                 "DGRP_852_M","DGRP_042_M")) # Fem -> Mal
ToPlot$dgrp_sex_tag <- factor(ToPlot$dgrp_sex_tag, 
                              levels = c("DGRP_383_F","DGRP_383_M",
                                         "DGRP_280_F","DGRP_280_M",
                                         "DGRP_852_F","DGRP_852_M",
                                         "DGRP_042_F","DGRP_042_M")) # DGRPx, DGRPy, ..z, ..a
ToPlot$dgrp <- factor(ToPlot$dgrp,
                     levels = c("DGRP_383","DGRP_280","DGRP_852","DGRP_042"))

pExtPairs <- ggplot(ToPlot, aes(x = dgrp,
                   y = rank_thes_prop_adj,
                   fill = sex)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c(col_fem, col_mal)) + 
  geom_point(position = position_jitterdodge(0.2)) +
  labs(title = "Most extreme and medicore lines", 
       x = "", 
       y = "Adjusted fraction of extremeness (1-x if x > 0.5)", 
       fill = "") +
  theme_bw() +
  theme(axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  scale_y_continuous(expand=c(0,0), limits = c(-.01,.51), breaks = seq(0,1,.1)) 

ggsave("Extremes/4D_most_extreme_pairs.pdf", pExtPairs,  units = "cm", 
       width = 35, height = 20, device = "pdf")
ggsave("Extremes/4D_most_extreme_pairs.png", pExtPairs,  units = "cm", 
       width = 35, height = 20, device = "png")









