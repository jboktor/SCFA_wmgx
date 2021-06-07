# Feature Correlation Analysis

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")
base::load("files/Phyloseq_objects_Woltka.RData")

#_______________________________________________________________________________
#####                      Correlation Analysis                           ##### 
#_______________________________________________________________________________

# Metadata Selection
corr.meta <-
  dat.species %>%
  process_meta() %>%
  dplyr::select(contains(c(motortest_vars)))


selected_dats <-
  c(
    "Species",
    "Genus",
    "Phylum",
    "MCProteins",
    "Enzrxn",
    "Superpathways",
    "Pathways",
    "KOs",
    "GOs",
    "eggNOGs"
  )

correlation_df <- tibble()
for (obj.name in selected_dats){

  print_line(); cat(obj.name, "\n\n")
  
  tested_feats <- 
    read_tsv(paste0("data/MaAsLin2_Analysis/", obj.name, 
                    "_treatment_group_MaAsLin2/all_results.tsv"), 
             col_names = T) %>% select(feature) %>% unique() %>% 
    mutate(feature = )
  cat("Testing a total of:", length(tested_feats$feature), obj.name, 
      "\n", sep = " ")
  
  corr.abund <- 
    phylo_dats[[obj.name]] %>% 
    abundances() %>% as.data.frame() %>% 
    rownames_to_column() %>% 
    mutate(rowname = gsub("\\-", ".", rowname)) %>% 
    mutate(rowname = gsub("\\:", ".", rowname)) %>% 
    mutate(rowname = gsub("\\+", ".", rowname)) %>% 
    column_to_rownames() %>% 
    t() %>% as.data.frame() %>% 
    dplyr::select(one_of(tested_feats$feature))
  
  corr.data <-
    corr_loop_parallel(metadata = corr.meta,
                       abundance = corr.abund,
                       obj.name = obj.name)
  
  # Scatter plots
  top_n_scatterplots(dat = phylo_dats[[obj.name]], obj.name = obj.name,
                     df.cors = corr.data)
  # Heatmaps
  heatmap <- corr_heatmap(corr.data)
  
  ggsave(heatmap, filename = 
           paste0("data/Correlations/", obj.name, "_motor_score_correlations_heatmap.svg"),
         height = 10, width = 13)
  
  correlation_df <- bind_rows(correlation_df, corr.data)
}

openxlsx::write.xlsx(correlation_df, file = 'data/Correlations/correlation_data.xlsx')
save(correlation_df, file = 'data/Correlations/correlation_data.RData')


