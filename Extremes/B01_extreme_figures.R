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

load(file = "Extremes/extreme_dgrp_fraction_thes.RData")
load(file = "Extremes/extreme_lines_perPheno_fraction.RData")
load(file = "Extremes/extreme_lines_both_sex_perPheno_fraction.RData")

# Set parameters ------------------------------------------ 

nr_phe <- 50

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
  geom_bar(stat = "identity") +
  scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1,.1)) +
  labs(title = "Female", x = "", y= "", fill="")

pM <- ggplot(ToPlotM, 
             aes(x = dgrp,
                 y = 1-fraction_of_group)) + 
  geom_bar(stat = "identity") +
  scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1,.1)) +
  labs(title = "Male", x = "", y= "", fill="")

pF + pM


# Individual extremes --------------------------------------------- [ True extremes ]

ToPlot <- extreme_lines_perPheno_fraction 
ToPlot$dgrp_sex_tag <- factor(ToPlot$dgrp_sex_tag, 
                              levels = c("DGRP_783_M","DGRP_352_M","DGRP_879_F","DGRP_757_F"))

ggplot(ToPlot, aes(x = dgrp_sex_tag,
                   y = rank_thes_prop_adj)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2))


# Individual extremes --------------------------------------------- [ True extremes (for both sexes) ]

ToPlot <- extreme_lines_both_sex_perPheno_fraction
ToPlot$dgrp_sex_tag <- factor(ToPlot$dgrp_sex_tag, 
                      levels = c("DGRP_383_F","DGRP_280_F",
                                 "DGRP_852_F","DGRP_042_F",
                                 "DGRP_383_M","DGRP_280_M",
                                 "DGRP_852_M","DGRP_042_M"))

ggplot(ToPlot, aes(x = dgrp_sex_tag,
                   y = rank_thes_prop_adj,
                   fill = sex)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2))










