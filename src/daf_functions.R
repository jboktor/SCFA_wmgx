# Differentially Abundant Features (DAF) Functions


#_______________________________________________________________________________


fit_models <- function(obj.name, 
                       df_input,
                       df_metadata,
                       cores = 10, 
                       plot_scatter = F){
  
  #' function to run MaAsLin2 mixed models
  #' Input variance and qc trimmed abundance data tables
  #' Outputs MaAsLin2 analysis in a specific feature folder
  
  fit_data = Maaslin2(
    input_data = df_input,
    input_metadata = df_metadata,
    output = paste0(wkd, "/data/MaAsLin2_Analysis/", obj.name, "_treatment_group_MaAsLin2"), 
    fixed_effects = c("genotype", "diet", "treatment_group"),
    reference = c("treatment_group,WT_Control"),
    min_prevalence = 0,
    analysis_method = "LM",
    normalization = "NONE",
    transform = "NONE",
    cores = cores,
    plot_scatter = plot_scatter)

}


#_______________________________________________________________________________

VennPlot <- function(Maaslin_PDvPC, Maaslin_PDvHC, qval_threshold = 0.25){
  
  ###' Function reads in maaslin2 significance results and creates Venn-Diagrams
  ###' for shared and unique features between Household & Population controls in 
  ###' reference to PD patients
  
  ## Note qval_threshold must be 0.25 or less
  require(eulerr)
  
  # Initalize Return variables
  venn_depleted <- NULL
  venn_enriched <- NULL
  
  ## Split both inputs into enriched and depleted features
  Maaslin_PDvPC <- Maaslin_PDvPC %>% filter(qval <= qval_threshold)
  Maaslin_PDvHC <- Maaslin_PDvHC %>% filter(qval <= qval_threshold)
  
  ## If either of the data frames have 0 features after qval-filtering, Stop and print error:
  if (nrow(Maaslin_PDvPC) < 1 | nrow(Maaslin_PDvHC) < 1){
    cat("One or both of comparisons doesn't meet the qval threshold \n")
    return()
  } else {
    cat("### Sufficient features for plotting \n")
  }
  
  PC.model <- mutate(Maaslin_PDvPC, direction = if_else(coef > 0, "Depleted","Enriched"))
  HC.model <- mutate(Maaslin_PDvHC, direction = if_else(coef > 0, "Depleted","Enriched")) 
  a.pc <- dplyr::filter(PC.model, direction == "Depleted")
  b.pc <- dplyr::filter(PC.model, direction == "Enriched")
  a.hc <- dplyr::filter(HC.model, direction == "Depleted")
  b.hc <- dplyr::filter(HC.model, direction == "Enriched")
  
  
  if (nrow(a.pc) > 0 & nrow(a.hc) > 0 ){
    ## Create joined matrix of features 
    down.vs.pc.df <- data.frame("features"=a.pc$feature, "Down_vs_PC" = TRUE)
    down.vs.hc.df <- data.frame("features"=a.hc$feature, "Down_vs_HC" = TRUE)
    depleted_matrix <-full_join(down.vs.pc.df, down.vs.hc.df,  by="features")
    depleted_matrix[is.na(depleted_matrix)] <-  F
    
    venn_depleted <- plot(euler(depleted_matrix[-1], shape = "ellipse"), quantities = TRUE)
    cat("1) PD DEPLETED FEATURES: Plotting shared PD depleted features \n")
    # print(venn_depleted)
    
  } else {
    cat("1) PD DEPLETED FEATURES: At least one comparison without features: no plot \n")
  }
  
  
  if (nrow(b.pc) > 0 & nrow(b.hc) > 0 ){
    
    ## Create joined matrix of features 
    up.vs.pc.df <- data.frame("features"=b.pc$feature, "Up_vs_PC" = T)
    up.vs.hc.df <- data.frame("features"=b.hc$feature, "Up_vs_HC" = T)
    enriched_matrix <-full_join(up.vs.pc.df, up.vs.hc.df,  by="features")
    enriched_matrix[is.na(enriched_matrix)] <-  F
    
    venn_enriched <- plot(euler(enriched_matrix[-1], shape = "ellipse"), quantities = TRUE)
    cat("2) PD ENRICHED FEATURES: Plotting shared PD enriched features \n")
    # print(venn_enriched)
    
  } else {
    cat("2) PD ENRICHED FEATURES: At least one comparison without features: no plot \n")
  }
  
  plots2return <- list( "venn_depleted" = venn_depleted, "venn_enriched" = venn_enriched)

  cat("Venn Diagram Plotting Complete: \n\n")
  return(plots2return)
  
}

#_______________________________________________________________________________

taxa_genus_phlyum_annotation = function(physeq, selectedTaxa){
  
  ###'  Returns Phylum level of a selected list 
  ###'  of species from phyloseq Object
  
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% selectedTaxa)]
  physeq <- prune_taxa(allTaxa, physeq)
  tax.df <- as.data.frame(tax_table(physeq)[,2])
  
  return(tax.df)
}


#_______________________________________________________________________________
# Inspired by MicrobiomeAnalystR: https://github.com/xia-lab/MicrobiomeAnalystR
# Plot (0.1 - 0.9) features by rank : to help decide on percentage cutoff 

variance_plot <- function(dat) {
  
  int.mat <- abundances(dat) %>% as.data.frame()
  filter.val <- apply(int.mat, 1, function (x) {
    diff(quantile(as.numeric(x), c(0.1, 0.9), na.rm = TRUE, names = FALSE, 
                  type = 7)) })
  
  var.df <- as.data.frame(filter.val) %>% 
    rownames_to_column(var = "features")
  
  rk <- rank(-filter.val, ties.method='random')
  rws <-  nrow(var.df)
  
  p <- ggplot(var.df, aes(x= reorder(features, -filter.val), y= log10(filter.val + .000000001))) +
    geom_point(color="#1170aa") + 
    ggthemes::theme_clean() +
    labs(x = "Ranked Features", y = " log10([0.1 - 0.9] Quantile Range)") +
    geom_vline(xintercept = c(rk[rk == round(rws*.9)], rk[rk == round(rws*.8)], rk[rk == round(rws*.7)], 
                              rk[rk == round(rws*.6)], rk[rk == round(rws*.5)], rk[rk == round(rws*.4)],
                              rk[rk == round(rws*.3)], rk[rk == round(rws*.2)]), linetype = "dashed", alpha = 0.7 ) +
    theme(axis.text.x = element_blank())
  return(p)
}


#_______________________________________________________________________________

variance_filter <- function(dat, filter.percent = 0.1) {
  
  #' This function filters features by their
  #' Inter-quartile range - larger values indicate larger spread
  #' features are ranked by IQR and a specified percentage is trimmed
  
  int.mat <- abundances(dat) %>% as.data.frame()
  filter.val <- apply(int.mat, 1, function (x) {
    diff(quantile(as.numeric(x), c(0.1, 0.9), na.rm = TRUE, names = FALSE, 
                  type = 7)) })
  
  rk <- rank(-filter.val, ties.method='random')
  var.num <- nrow(int.mat);
  remain <- rk <= var.num*(1-filter.percent);
  int.mat <- int.mat[remain,];
  
  cat("A total of", sum(!remain), "low variance features were removed based on the Quantile Range between [0.1 - 0.9]. \n")
  cat("The number of features remaining after filtering is:", nrow(int.mat), "\n")
  
  return(int.mat)
  
}

#_______________________________________________________________________________
#                         Genotype by Diet Summary
#_______________________________________________________________________________

dafs_genotype_by_diet <- function(obj.name, n_lab_diet = 3, nlab_geno = 3){
  
  ### Read-in MaAsLin2 output
  maaslin_df <-
    read_tsv(paste0("data/MaAsLin2_Analysis/",
                    obj.name,
                    "_treatment_group_MaAsLin2/all_results.tsv"),
             col_names = T)
  
  treatments_stats <- 
    maaslin_df %>% 
    select(feature, value, maaslin_stats) %>%
    pivot_wider(names_from = "value", values_from = maaslin_stats) %>% 
    mutate(diet_treatment_sig = if_else(qval_Prebiotic < 0.1 | 
                                          qval_ASO < 0.1, T, F)) %>% 
    mutate(feat_label = if_else(grepl("s__", feature), 
                                sub(".*\\s__", "", feature), feature)) %>% 
    mutate(feat_label = if_else(grepl("g__", feature), 
                                sub(".*\\.g__", "", feature), feature))
  
  
  treatments_stats_dietsig <- treatments_stats %>% 
    filter(diet_treatment_sig) %>% 
    slice_max(order_by = abs(coef_Prebiotic), n = n_lab_diet)
  treatments_stats_genosig <- treatments_stats %>% 
    filter(diet_treatment_sig) %>% 
    slice_max(order_by = abs(coef_ASO), n = nlab_geno)
  treatments_stats_labs <- 
    bind_rows(treatments_stats_dietsig, treatments_stats_genosig) %>% 
    distinct()
  
  
  treatments_stats %>% 
    ggplot(aes(x = coef_Prebiotic, y = coef_ASO)) +
    geom_point(aes(color = diet_treatment_sig), alpha = 0.6) +
    ggplot2::annotate("text", x=0.33, y=0.33, label= "ASO\nPrebiotic", color = "#45403f") +
    ggplot2::annotate("text", x=0.33, y=-0.33, label= "ASO\nControl", color = "#45403f") +
    ggplot2::annotate("text", x=-0.33, y=0.33, label= "WT\nPrebiotic", color = "#45403f") +
    ggplot2::annotate("text", x=-0.33, y=-0.33, label= "WT\nControl", color = "#45403f") +
    geom_abline(intercept = 0, slope = 1, linetype = 3, color = "darkgrey") +
    geom_abline(intercept = 0, slope = -1, linetype = 3, color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = 1, color = "grey") +
    geom_hline(yintercept = 0, linetype = 1, color = "grey") +
    lims(x = c(-0.4, 0.4), y = c(-0.4, 0.4)) +
    labs(y = expression(paste("GLM ", beta, " coefficient Diet [Prebiotic/Control]")),
         x = expression(paste("GLM ", beta, " coefficient Genotype [ASO/WT]"))) +
    scale_color_manual(values = c('#45403f', '#CC071E')) +
    geom_text_repel(
      data = treatments_stats_labs,
      aes(x = coef_Prebiotic, y = coef_ASO, label = feat_label),
      segment.color	 =  "#929292",
      segment.curvature = -0.5,
      # angle = 90,
      # nudge_y = 0.2,
      # direction = "x", 
      # vjust = 0.1,
      force = 5, 
      size = 3) +
    my_clean_theme() +
    theme(legend.position = "none")
}





