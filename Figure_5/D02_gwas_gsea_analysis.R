# GSEA Function -------------------------------------------
# As provided by https://bioinformaticsbreakdown.com/how-to-gsea/
GSEA <- function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=15, ## minimum gene set size
                        maxSize=400, ## maximum gene set size
                        nperm=10000) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header) 
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}

# GMT File location
GMT_file <- "Extremes/lifespan_genes_GO0008340.gmt"
inputGMT <- gmtPathways(GMT_file) 

# Set gene lists ------------------------------------------ [ Run models ]
# Female Std 
input_list_fem_std <- allhits_fem_ranked_with_clean_genes_lfsp_marked %>% 
  select(ID, BETA_std, BETA_foe, P_std, P_foe, pval_compare, FBgn_ID) %>% 
  arrange(P_std) %>% 
  add_column(rank = 1:nrow(.)) %>% 
  filter(is.na(FBgn_ID)==F) %>% 
  distinct(FBgn_ID,P_std) %>% 
  group_by(FBgn_ID) %>% 
  filter(P_std == max(P_std)) %>% 
  ungroup() %>% 
  mutate(P_std = -log10(P_std)) %>% 
  pivot_wider(names_from = FBgn_ID,
              values_from = P_std)
input_list_fem_std <- sort(input_list_fem_std, decreasing = T)

set.seed(54321)
gsea_fem_std <- fgsea(pathways = inputGMT,
                      stats = input_list_fem_std,
                      minSize=15, ## minimum gene set size
                      maxSize=400, ## maximum gene set size
                      nperm=10000) %>% 
  as_tibble()

# Female FoE 
input_list_fem_foe <- allhits_fem_ranked_with_clean_genes_lfsp_marked %>% 
  select(ID, BETA_std, BETA_foe, P_std, P_foe, pval_compare, FBgn_ID) %>% 
  arrange(P_foe) %>% 
  add_column(rank = 1:nrow(.)) %>% 
  filter(is.na(FBgn_ID)==F) %>% 
  distinct(FBgn_ID,P_foe) %>% 
  group_by(FBgn_ID) %>% 
  filter(P_foe == max(P_foe)) %>% 
  ungroup() %>% 
  mutate(P_foe = -log10(P_foe)) %>% 
  pivot_wider(names_from = FBgn_ID,
              values_from = P_foe)
input_list_fem_foe <- sort(input_list_fem_foe, decreasing = T)

set.seed(54321)
gsea_fem_foe <- fgsea(pathways = inputGMT,
                      stats = input_list_fem_foe,
                      minSize=15, ## minimum gene set size
                      maxSize=400, ## maximum gene set size
                      nperm=10000) %>% 
  as_tibble()

# Male Std 
input_list_mal_std <- allhits_mal_ranked_with_clean_genes_lfsp_marked %>% 
  select(ID, BETA_std, BETA_foe, P_std, P_foe, pval_compare, FBgn_ID) %>% 
  arrange(P_std) %>% 
  add_column(rank = 1:nrow(.)) %>% 
  filter(is.na(FBgn_ID)==F) %>% 
  distinct(FBgn_ID,P_std) %>% 
  group_by(FBgn_ID) %>% 
  filter(P_std == max(P_std)) %>% 
  ungroup() %>% 
  mutate(P_std = -log10(P_std)) %>% 
  pivot_wider(names_from = FBgn_ID,
              values_from = P_std)
input_list_mal_std <- sort(input_list_mal_std, decreasing = T)

set.seed(54321)
gsea_mal_std <- fgsea(pathways = inputGMT,
                      stats = input_list_mal_std,
                      minSize=15, ## minimum gene set size
                      maxSize=400, ## maximum gene set size
                      nperm=10000) %>% 
  as_tibble()

# Male FoE 
input_list_mal_foe <- allhits_mal_ranked_with_clean_genes_lfsp_marked %>% 
  select(ID, BETA_std, BETA_foe, P_std, P_foe, pval_compare, FBgn_ID) %>% 
  arrange(P_foe) %>% 
  add_column(rank = 1:nrow(.)) %>% 
  filter(is.na(FBgn_ID)==F) %>% 
  distinct(FBgn_ID,P_foe) %>% 
  group_by(FBgn_ID) %>% 
  filter(P_foe == max(P_foe)) %>% 
  ungroup() %>% 
  mutate(P_foe = -log10(P_foe)) %>% 
  pivot_wider(names_from = FBgn_ID,
              values_from = P_foe)
input_list_mal_foe <- sort(input_list_mal_foe, decreasing = T)

set.seed(54321)
gsea_mal_foe <- fgsea(pathways = inputGMT,
                      stats = input_list_mal_foe,
                      minSize=15, ## minimum gene set size
                      maxSize=400, ## maximum gene set size
                      nperm=10000) %>% 
  as_tibble()

# Summarise GSEA results ----------------------------------
gsea_summary <- bind_rows(gsea_fem_std %>% 
                            add_column(sex = "Female",
                                       model = "std",
                                       .before = 1),
                          gsea_fem_foe %>% 
                            add_column(sex = "Female",
                                       model = "foe",
                                       .before = 1),
                          gsea_mal_std %>% 
                            add_column(sex = "Male",
                                       model = "std",
                                       .before = 1),
                          gsea_mal_foe %>% 
                            add_column(sex = "Male",
                                       model = "foe",
                                       .before = 1)) %>% 
  mutate(leadingEdge = unlist(leadingEdge))

# Plotting ------------------------------------------------ [ Single pathway, meaningless plot ]
gsea_output <- fgsea(pathways = inputGMT,
              stats = input_list_fem_std,
              minSize=15, ## minimum gene set size
              maxSize=400, ## maximum gene set size
              nperm=10000) %>% 
  as_tibble() 

ggplot(gsea_output, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(size = size), shape=21) +
  #scale_fill_manual(values = c("firebrick2", "dodgerblue2") ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= "gsea_output")


