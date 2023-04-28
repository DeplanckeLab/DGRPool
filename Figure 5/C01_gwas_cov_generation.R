dgrp_cov <- read_tsv("dgrp.cov.tsv")

dgrp_foe <- extreme_dgrp_fraction_thes_summ %>% 
  filter(total_nr_phenotypes >= 50) %>% 
  mutate(covID = str_replace(dgrp,"DGRP","line")) %>% 
  mutate(covID = str_replace(covID,"_00","_")) %>% 
  mutate(covID = str_replace(covID,"_0","_")) %>% 
  rename(foe = "fraction_extreme")

dgrp_foe %>% 
  count(sex)

dgrp_cov_F <- dgrp_cov %>% 
  left_join(dgrp_foe %>% 
              filter(sex == "F") %>% 
              select(covID, foe),
            by = c("FID" = "covID")) %>% 
  write_tsv("dgrp.cov.f.foe.tsv")

dgrp_cov_M <- dgrp_cov %>% 
  left_join(dgrp_foe %>% 
              filter(sex == "M") %>% 
              select(covID, foe),
            by = c("FID" = "covID")) %>% 
  write_tsv("dgrp.cov.m.foe.tsv")

