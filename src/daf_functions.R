# Differentially Abundant Features (DAF) Functions


#-------------------------------------------------------------------------------


fit_models <- function(dat, dat_pdpc, dat_pdhc, obj.name = "Species", 
                       df_input_data_pdhc, df_input_data_pdpc, cores, plot_scatter = F, 
                       cohort = "Merged"){
  
  #' function to run MaAsLin2 mixed models
  #' Input variance and qc trimmed abundance data tables
  #' Outputs MaAsLin2 analysis in a specific feature folder
  
  # Run Metadata pre-processing function
  # df_input_metadata <- 
  #   process_meta(dat, cohort = cohort) %>% 
  #   column_to_rownames(var = "donor_id")
  
  df_input_metadata_pdpc <- 
    process_meta(dat_pdpc, cohort = cohort) %>% 
    dplyr::mutate(description = factor(description, levels = c("PD Patient", "Population Control"))) %>% 
    column_to_rownames(var = "donor_id")
  
  df_input_metadata_pdhc <- 
    process_meta(dat_pdhc, cohort = cohort) %>% 
    dplyr::mutate(description = factor(description, levels = c("PD Patient", "Household Control"))) %>% 
    column_to_rownames(var = "donor_id")
  
  if (cohort == "Merged"){
    random_effects_pdpc <- "cohort"
    random_effects_pdhc <- c("paired", "cohort")
  } else {
    random_effects_pdpc <- ""
    random_effects_pdhc <- "paired"
  }
  
  # PD v PC  
  fit_data = Maaslin2(
    input_data = df_input_data_pdpc,
    input_metadata = df_input_metadata_pdpc,
    output = paste0(wkd, "/data/MaAsLin2_Analysis/", cohort, "/", obj.name, "_PDvPC_maaslin2_output"), 
    random_effects = random_effects_pdpc,
    fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
    min_prevalence = 0,
    analysis_method = "LM",
    normalization = "NONE",
    transform = "AST",
    cores = cores,
    plot_scatter = plot_scatter
  )
  
  # PD v HC  
  fit_data = Maaslin2(
    input_data = df_input_data_pdhc, 
    input_metadata = df_input_metadata_pdhc, 
    output = paste0(wkd, "/data/MaAsLin2_Analysis/", cohort, "/", obj.name, "_PDvHC_maaslin2_output"), 
    min_prevalence = 0,
    random_effects = random_effects_pdhc,
    fixed_effects = c("description"),
    analysis_method = "LM",
    normalization = "NONE",
    transform = "AST",
    cores = cores,
    plot_scatter = plot_scatter
  )
}


#-------------------------------------------------------------------------------

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

#---------------------------------------------------------------------------------------------------------------------

taxa_genus_phlyum_annotation = function(physeq, selectedTaxa){
  
  ###'  Returns Phylum level of a selected list 
  ###'  of species from phyloseq Object
  
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% selectedTaxa)]
  physeq <- prune_taxa(allTaxa, physeq)
  tax.df <- as.data.frame(tax_table(physeq)[,2])
  
  return(tax.df)
}


#---------------------------------------------------------------------------------------------------------------------

pull.phylo.sig <- function(comparison, phylo.names){
  
  ###' Function takes in a comparion (either HC or PC) and table
  ###' containing Genus and Phylum annotations 
  ###' Returns long format table with features and their qvals from maaslin
  ###' at the Genus and Phylum levels
  
  output <- NULL
  
  if (comparison == "PC"){
    groop <- "Population Control"
  } 
  else if (comparison == "HC"){
    groop <- "Household Control"
  } 
  else {
    return("ERROR: Please enter either PC or HC for comparison: ")
  }
  
  ### Read-in Maaslin Files - PHYLUM LEVEL
  Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/Phylum_PDv", comparison, "_maaslin2_output/all_results.tsv"), col_names = T) %>% 
    filter(value == groop)
  Maas.pd.pc$feature <- gsub("s__", "", Maas.pd.pc$feature)
  Phylum.qvals <- Maas.pd.pc %>% filter(feature %in% unique(phylo.names$Phylum))
  vals <- Phylum.qvals %>% dplyr::select(feature, qval)
  
  colnames(vals)[1] <- "Phylum"
  return(vals)
  
}

#---------------------------------------------------------------------------------------------------------------------

boxplot_phylobars <- function(inpt.phylo, sigvals, tile.cols){
  
  ##' Function takes in features selected for plotting 
  ##' Returns color bars, legends, and y-axis order separately, 
  ##' mapping each feature's Phylum and Genus
  ##' Adds significance symbol if below threshold in respective
  ##' MaAsLin2 model
  
  full_plot <- NULL
  Legends <- NULL
  products <- NULL
  
  phylo.plot <- inpt.phylo %>% rownames_to_column(var = "Species")
  phylo.plot$Species <- gsub("s__", "", phylo.plot$Species)
  phylo.plot <-  left_join(phylo.plot, sigvals, by = "Species")
  
  ## Sorting - Phylum then qval
  phylo.plot <- phylo.plot[order(phylo.plot$Phylum, phylo.plot$qval), ]
  phylo.plot$order_col <- 1:nrow(phylo.plot)
  phylo.plot$xaxis <- "Phylum"
  
  ## Select colors from pre-made list
  cut_colors <- tile.cols[names(tile.cols) %in% phylo.plot$Phylum]

  full_plot <- 
    ggplot(phylo.plot, aes(x= xaxis, y=reorder(Species, -order_col), fill= Phylum)) + 
    geom_tile() +
    theme_minimal() +
    scale_fill_manual(values=cut_colors) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none")
  
  Phylum_legend <-
    phylo.plot %>%
    ggplot(aes(x=xaxis, y=reorder(Species, -order_col), fill= Phylum)) +
    geom_tile() +
    scale_fill_manual(values = cut_colors, name = "Phylum")
  
  Legends <- cowplot::plot_grid(get_legend(Phylum_legend))
  
  products <- list("Bars" = full_plot, "Legends" = Legends, "Axis.order" = unique(phylo.plot$Species))
  
  return(products)
  
}

#---------------------------------------------------------------------------------------------------------------------


boxplot_hl1.bars <- function(inpt.hl1, tile.cols){
  
  ##' Function takes in features selected for plotting 
  ##' Returns color bars, legends, and y-axis order separately, 
  ##' mapping each feature's Metabolic Module and it's higher level hierarchy
  ##' Adds significance symbol if below threshold in respective
  ##' MaAsLin2 model
  
  full_plot <- NULL
  Legends <- NULL
  products <- NULL
  
  ## Sorting - HL1 then qval
  hl1.plot <- inpt.hl1[order(inpt.hl1$HL1, inpt.hl1$qval), ]
  hl1.plot$order_col <- 1:nrow(hl1.plot)
  hl1.plot$xaxis <- "level_1"
  
  ## Select colors from pre-made list
  cut_colors <- tile.cols[names(tile.cols) %in% hl1.plot$HL1]
  
  full_plot <- 
    ggplot(hl1.plot, aes(x= xaxis, y=reorder(Name, -order_col), fill= HL1)) + 
    geom_tile() +
    theme_minimal() +
    scale_fill_manual(values=cut_colors) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none")
  
  Phylum_legend <-
    hl1.plot %>%
    ggplot(aes(x=xaxis, y=reorder(Name, -order_col), fill= HL1)) +
    geom_tile() +
    scale_fill_manual(values = cut_colors, name = "HL1")
  
  Legends <- cowplot::plot_grid(get_legend(Phylum_legend))
  
  products <- list("Bars" = full_plot, "Legends" = Legends, "Axis.order" = unique(hl1.plot$Name))
  
  return(products)
  
}


#---------------------------------------------------------------------------------------------------------------------
## Adapted from: https://github.com/zellerlab/crc_meta


generalized_fold_change <- function(pd_abundance, ctrl_abundances) {
  # Initialize counts
  i = 1
  probs.fc <- seq(.1, .9, .05)
  gfc_data <- c()
  
  for (feature in 1:nrow(ctrl_abundances)) {
    # Loops through each feature row and calculates
    # the Generalized fold change for each feature
    cat("Feature number: ", feature, "\n")
    cat(
      "Testing PD: ",
      rownames(pd_abundance)[feature],
      "feature vs ",
      rownames(ctrl_abundances)[feature],
      "\n"
    )
    q.pd <- quantile(pd_abundance[feature, ], probs = probs.fc)
    q.ctl <- quantile(ctrl_abundances[feature, ], probs = probs.fc)
    
    gfc <- sum(q.pd - q.ctl) / length(q.ctl)
    print(gfc)
    
    gfc_data <- c(gfc_data, gfc)
  }
  return(gfc_data)
}



#---------------------------------------------------------------------------------------------------------------------

gfc_plot <- function(df, manual_colors, alfa = 0.5){
  ######' Generalized Fold Change (gFC) BarPlot ######
  p <- ggplot(data=df, aes(x=gFC, y= feature, fill = direction)) +
    geom_point(aes(fill=direction), shape=21, size=3, alpha = alfa) +
    geom_segment(aes(x=0, xend=gFC, y=feature, yend=feature, color=direction)) +
    theme_minimal() +
    xlim(-max(abs(df$gFC)), max(abs(df$gFC))) +
    labs(x="Average Difference") +
    ggtitle("Generalized Fold Change") +
    scale_fill_manual(values = manual_colors, name ="Group") +
    scale_color_manual(values = manual_colors, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank())
  
  return(p)
}

#---------------------------------------------------------------------------------------------------------------------

significance_barplot <- function(df){
  ######' Significance BarPlot ######
  p <- ggplot() +
    geom_col(data=df, aes(x=-log10(value), y= feature,  fill = variable), position = "dodge", width = 0.8) +
    theme_minimal() +
    labs(x=expression('-log'[10]*'(value)')) +
    ggtitle("Significance") +
    scale_fill_manual(values = c("pval" = "#d3d3d3", "qval" = "#676767")) +
    geom_vline(xintercept = -log10(0.1), linetype = "dotted", color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank())
  return(p)
}

#---------------------------------------------------------------------------------------------------------------------

daf_boxplot_sigvalues <- function(sigplot.df, abund.input){
  
  sigplot.df.QVAL <- filter(sigplot.df, variable == "qval")
  sigplot.df.QVAL$value <- round(sigplot.df.QVAL$value, digits = 3)
  
  sig.labs <- c()
  for (i in sigplot.df.QVAL$value){
    sig.labs <- c(sig.labs, sig_mapper(i, porq = "q", symbols = F))
  }
  sigplot.df.QVAL$sig.labels <- sig.labs
  sigplot.df.QVAL <- sigplot.df.QVAL %>% dplyr::select(feature, sig.labels) %>% 
    dplyr::rename(Var2 = feature) 
  abund.input2 <- left_join(abund.input, sigplot.df.QVAL, by = "Var2")
  
  return(abund.input2)
  
}

#---------------------------------------------------------------------------------------------------------------------

daf_boxplots <- function(df, fill_cols, rim_cols, alfa = 0.5, obj.name){
  
  set.seed(123)
  p <- ggplot(data=df, aes(x=value, y= Var2)) +
    geom_boxplot(aes(fill = group), alpha = alfa, outlier.alpha = 0, width = 0.8) +
    geom_point(aes(fill=group, color=group), position = position_jitterdodge(jitter.width = .2), shape=21, size=1, alpha = 0.9) +
    theme_minimal() +
    ggtitle(paste0("Differential Abundance: ", obj.name)) +
    geom_text(aes(x= max(value) + max(value)*0.15, 
                  y=Var2, label = sig.labels), size = 4, check_overlap = TRUE) +
    labs(x = "Abundance") +
    scale_fill_manual(values = fill_cols, name ="Group") +
    scale_color_manual(values = rim_cols, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank())
  
  return(p)
}

#---------------------------------------------------------------------------------------------------------------------

prevalence_barplot <- function(df, manual_colors, alfa = 0.5){
  p <- ggplot(data=df, aes(x=value*100, y= feature, fill=variable, color=variable)) +
    geom_bar(position = position_dodge(),  alpha = alfa, width = 0.8,
             stat = "identity") +
    theme_minimal() +
    labs(x=expression('Prevalence [%]')) +
    ggtitle("Prevalence") +
    scale_fill_manual(values = manual_colors, name ="Group") +
    scale_color_manual(values = manual_colors, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank())

  return(p)
}


#-----------------------------------------------------------------------------------------------------------
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


#-----------------------------------------------------------------------------------------------------------

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


#---------------------------------------------------------------------------------------------------------------------


## FUNCTION TO ADD GROUP COLUMNS FOR D3 FlashWeave NETWORK

flashweave_group_colors <- function(features, metadata_added, Maas.pd.pc.sig, Maas.pd.hc.sig){
  
  ##' Loop for identifying significant features and stratifying them
  ##' by enrichment or depletion in PC, HC, or both groups in reference 
  ##' to PD groups 
  ##' Note: a negative coefficient in Maaslin model indicates higher levels in PD
  
  features <- features %>% 
    mutate(group = if_else(Node %in% metadata_added, "Donor_Group",
                           if_else(Node %in% Maas.pd.pc.sig$feature & Node %in% Maas.pd.hc.sig$feature, "BOTH",
                                   if_else(Node %in% Maas.pd.pc.sig$feature, "PC",
                                           if_else(Node %in% Maas.pd.hc.sig$feature, "HC",
                                                   "None")))))
  # Initiate vars for loop
  n = 1
  group2 <- c()
  for (node in features$Node){
    if (features$group[n] == "BOTH") {
      if (filter(Maas.pd.pc.sig, feature == node)$coef  < 0 & 
          filter(Maas.pd.hc.sig, feature == node)$coef  < 0) {
        group2  = c(group2, "Up_PDvBoth")
      }  else if (filter(Maas.pd.pc.sig, feature == node)$coef > 0 & 
                  filter(Maas.pd.hc.sig, feature == node)$coef  > 0) {
        group2  = c(group2, "Down_PDvBoth")
      }
    } else if (features$group[n] == "PC") {
      if (filter(Maas.pd.pc.sig, feature == node)$coef < 0) {
        group2  = c(group2, "Up_PDvPC")
      } else if (filter(Maas.pd.pc.sig, feature == node)$coef > 0) {
        group2  = c(group2, "Down_PDvPC")
      }
    } else if (features$group[n] == "HC") {
      if (filter(Maas.pd.hc.sig, feature == node)$coef < 0) {
        group2  = c(group2, "Up_PDvHC")
      }  else if (filter(Maas.pd.hc.sig, feature == node)$coef > 0) {
        group2  = c(group2, "Down_PDvHC")
      }
    } else if (features$group[n] == "Donor_Group") {
      group2  = c(group2, "Donor_Group")
    } else if (features$group[n] == "None") {
      group2  = c(group2, "None")
    } else {
      group2 = c(group2, "error")
    }
    n = n + 1
  }
  features$group <- group2
  
  ## Complex color codes - DAF PD/PC (up & down) and PD/HC feature (up and down) coloring -
  features2 <- features %>%
    mutate(group_color = if_else(group == "Donor_Group", "#d62728", 
                                 if_else(group == "Up_PDvBoth", "#98df8a",
                                         if_else(group == "Down_PDvBoth", "#2ca02c",
                                                 if_else(group == "Up_PDvPC", "#ffbb78",
                                                         if_else(group == "Down_PDvPC", "#ff7f0e",
                                                                 if_else(group == "Up_PDvHC", "#aec7e8",
                                                                         if_else(group == "Down_PDvHC", "#1f77b4",
                                                                                 if_else(group == "None", "#c7c7c7",
                                                                                         "Error")))))))))
  
  return(features2)
  
}


#---------------------------------------------------------------------------------------------------------------------
## FUNCTION TO ADD GROUP COLUMNS FOR DAF Fold Change Summary Plot

daf_group_colors <- function(features, Maas.pd.pc.sig, Maas.pd.hc.sig){
  
  ##' Loop for identifying significant features and stratifying them
  ##' by enrichment or depletion in PC, HC, or both groups in reference 
  ##' to PD groups 
  ##' Note: a negative coefficient in Maaslin model indicates higher levels in PD
  
  
  # ## TROUBLE
  # features <- taxa_names(dat.GOs.slim)
  # features <- tibble("Node" = features)
  # 

  features <- features %>% 
    mutate(group = if_else(Node %in% Maas.pd.pc.sig$feature & 
                             Node %in% Maas.pd.hc.sig$feature, "BOTH",
                           if_else(Node %in% Maas.pd.pc.sig$feature, "PC", 
                                   if_else(Node %in% Maas.pd.hc.sig$feature, "HC",
                                                   "None"))))
  # Initiate vars for loop
  n = 1
  group2 <- c()
  for (node in features$Node){
    if (features$group[n] == "BOTH") {
      if (filter(Maas.pd.pc.sig, feature == node)$coef  < 0 & 
          filter(Maas.pd.hc.sig, feature == node)$coef  < 0) {
        group2  = c(group2, "PD Enriched /PC & HC")
      }  else if (filter(Maas.pd.pc.sig, feature == node)$coef > 0 & 
                  filter(Maas.pd.hc.sig, feature == node)$coef  > 0) {
        group2  = c(group2, "PD Depleted /PC & HC")
      }
    } else if (features$group[n] == "PC") {
      if (filter(Maas.pd.pc.sig, feature == node)$coef < 0) {
        group2  = c(group2, "PD Enriched /PC")
      } else if (filter(Maas.pd.pc.sig, feature == node)$coef > 0) {
        group2  = c(group2, "PD Depleted /PC")
      }
    } else if (features$group[n] == "HC") {
      if (filter(Maas.pd.hc.sig, feature == node)$coef < 0) {
        group2  = c(group2, "PD Enriched /HC")
      }  else if (filter(Maas.pd.hc.sig, feature == node)$coef > 0) {
        group2  = c(group2, "PD Depleted /HC")
      }
    } else if (features$group[n] == "Donor_Group") {
      group2  = c(group2, "Donor_Group")
    } else if (features$group[n] == "None") {
      group2  = c(group2, "None")
    } else {
      group2 = c(group2, "error")
    }
    n = n + 1
  }
  features$group_directional <- group2
  return(features)
  
}

#-------------------------------------------------------------------------------
#                   DAF Summary Panel
#-------------------------------------------------------------------------------

plot_dafs <- function(obj.name, obj, cohort = "TBC", tag = ""){
  
  load("files/low_quality_samples.RData")

  # # TROUBLE
  # obj.name = "GOs.slim"
  # obj = dat.object
  # cohort = "Merged"
  
  obj <- obj %>% 
    subset_samples(donor_id %ni% low_qc[[1]])
  
  abund_rename <-
    obj %>% 
    abundances() %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    decode_rfriendly_rows(passed_column = "rowname") %>% 
    column_to_rownames(var = "fullnames") %>%
    select(-rowname) %>% 
    otu_table(taxa_are_rows=T)
  my_sample_data <- meta(obj) %>% sample_data()
  obj <- phyloseq(abund_rename, my_sample_data)
  
  
  ## Color Schemes
  cols.pdpc <- c("PD"= "#bfbfbf", "PC" = "#ed7d31")
  cols.pdhc <- c("PD"= "#bfbfbf", "HC" = "#5b9bd5")
  # Rims
  cols.pdpc.rim <- c("PD"= "#494949", "PC" = "#c15811")
  cols.pdhc.rim <- c("PD"= "#494949", "HC" = "#2e75b5")

  # Normalization and Transformation
  dat_obj <- microbiome::transform(obj, "compositional")
  otu_table(dat_obj) <- asin(sqrt(otu_table(dat_obj)))
  
  # PD v PC
  dat_pdpc = subset_samples(dat_obj, donor_group !="HC")
  abun.pdpc <- as.data.frame.matrix(abundances(dat_pdpc))
  # PD v HC 
  dat_pdhc = subset_samples(dat_obj, paired !="No")
  abun.pdhc <- as.data.frame.matrix(abundances(dat_pdhc))
  
  ### Read-in MaAsLin2 output
  Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/", cohort, "/", obj.name, "_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
    filter(value == "Population Control")
  Maas.pd.pc.sig <- Maas.pd.pc %>% filter(qval < 0.25) %>% 
    decode_rfriendly_rows(passed_column = "feature") %>% 
    select(-feature) %>%
    dplyr::rename("feature"="fullnames")
  
  Maas.pd.hc <- read_tsv(paste0("data/MaAsLin2_Analysis/",  cohort, "/", obj.name, "_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
    filter(value == "Household Control")
  Maas.pd.hc.sig <- Maas.pd.hc %>% filter(qval < 0.25) %>% 
    decode_rfriendly_rows(passed_column = "feature") %>% 
    select(-feature) %>%
    dplyr::rename("feature"="fullnames")
  
  #--------------------------------------
  #           Venn Diagrams
  #--------------------------------------

  v <- VennPlot(Maas.pd.pc.sig, Maas.pd.hc.sig, qval_threshold = 0.25)
  
  # Save Venn Diagrams
  if (!is.null(v$venn_depleted)){
    pdf(file = paste0("data/DAF_Analysis/", cohort, "/DAF_", obj.name, "_VennDiagram_PD_depleted.pdf"),
        width = 7, 
        height = 5,
        pointsize = 12)
    plot(v$venn_depleted)
    dev.off()
  }
  
  if (!is.null(v$venn_enriched)){
    pdf(file = paste0("data/DAF_Analysis/", cohort, "/DAF_", obj.name, "_VennDiagram_PD_enriched.pdf"),
        width = 7, 
        height = 5,
        pointsize = 12)
    plot(v$venn_enriched)
    dev.off()
  }
  
  #--------------------------------------
  
  # Select top 50 features by q-value
  if (nrow(Maas.pd.pc.sig) > 25){
    Maas.pd.pc.sig <- Maas.pd.pc.sig[1:25,]
  }
  if (nrow(Maas.pd.hc.sig) > 25){
    Maas.pd.hc.sig <- Maas.pd.hc.sig[1:25,]
  }
  
  # Pull significant features from abundance tables 
  # PD v PC
  abun.pdpc.inpt <-
    abun.pdpc %>% 
    rownames_to_column() %>% 
    filter(rowname %in% Maas.pd.pc.sig$feature) %>%
    column_to_rownames(var = "rowname") %>%
    t() %>% melt() %>% 
    mutate(group = if_else(grepl("PC", Var1), "PC", "PD"))
  # PD v HC PAIRED
  abun.pdhc.inpt <- 
    abun.pdhc %>% 
    rownames_to_column() %>% 
    filter(rowname %in% Maas.pd.hc.sig$feature) %>%
    column_to_rownames(var = "rowname") %>%
    t() %>% melt() %>% 
    mutate(group = if_else(grepl("HC", Var1), "HC", "PD"))
  
  #--------------------------------------------------------------------------
  #                              PD v PC Plots 
  #--------------------------------------------------------------------------
  
  # Subset of phyloseq obj subset to get samples of interest
  dat_pdpc.PD = subset_samples(dat_pdpc, donor_group =="PD")
  dat_pdpc.PDabun <- as.data.frame.matrix(abundances(dat_pdpc.PD)) %>% 
    rownames_to_column() %>%  filter(rowname %in% Maas.pd.pc.sig$feature) %>% column_to_rownames()
  
  dat_pdpc.PC = subset_samples(dat_pdpc, donor_group =="PC")
  dat_pdpc.PCabun <- as.data.frame.matrix(abundances(dat_pdpc.PC)) %>% 
    rownames_to_column() %>%  filter(rowname %in% Maas.pd.pc.sig$feature) %>% column_to_rownames()
  
  
  ######  Generalized or pseudo-fold calculation  ######
  gfc_data <- generalized_fold_change(dat_pdpc.PDabun, dat_pdpc.PCabun)
  PDovrPC <- tibble("feature" = rownames(dat_pdpc.PCabun), "gFC" = gfc_data)
  ###### To get order for Significant feauture plots
  ### Stop chunk after next line to order by gFC
  PDovrPC.order <- PDovrPC[order(-PDovrPC$gFC),]
  ### Use for ordering by qval
  PDovrPC.order <- left_join(PDovrPC.order, Maas.pd.pc.sig, by = "feature")
  PDovrPC.order <- PDovrPC.order[order(-PDovrPC.order$qval),]
  
  
  ###### Generalized Fold Change (gFC) BarPlot ######
  PDovrPC.BP <- PDovrPC.order
  PDovrPC.BP <- mutate(PDovrPC.BP, direction = if_else(PDovrPC.BP$gFC > 0, "PD",
                                                       if_else(PDovrPC.BP$gFC < 0, "PC",  "error")))
  PDovrPC.BP$feature <- factor(PDovrPC.BP$feature, levels = PDovrPC.order$feature)
  g0 <- gfc_plot(PDovrPC.BP, cols.pdpc, alfa = 0.8)
  
  
  ###### Significance Plot ######
  sigplot.df.pdpc <- dplyr::select(Maas.pd.pc.sig, c("feature", "pval", "qval")) %>%  melt()
  sigplot.df.pdpc$feature <- factor(sigplot.df.pdpc$feature, levels = PDovrPC.order$feature)
  g2 <- significance_barplot(sigplot.df.pdpc)
  
  
  ###### BoxPlots ######
  ## Prepping Significance labels
  abun.pdpc.inpt <- daf_boxplot_sigvalues(sigplot.df.pdpc, abun.pdpc.inpt)
  abun.pdpc.inpt$Var2 <- factor(abun.pdpc.inpt$Var2, levels = PDovrPC.order$feature) 
  g1 <- daf_boxplots(abun.pdpc.inpt, fill_cols = cols.pdpc, rim_cols = cols.pdpc.rim, alfa = 0.2, obj.name = obj.name)
  
  
  ###### Prevalence Plot ######
  # Subset of phyloseq obj subset to get samples of interest
  dat_pdpc.PD = subset_samples(dat_pdpc, donor_group =="PD")
  dat_pdpc.PDprev <- tibble::enframe(prevalence(dat_pdpc.PD)) %>% filter(name %in% Maas.pd.pc.sig$feature) 
  colnames(dat_pdpc.PDprev) <- c("feature", "PD")
  dat_pdpc.PC = subset_samples(dat_pdpc, donor_group =="PC")
  dat_pdpc.PCprev <- tibble::enframe(prevalence(dat_pdpc.PC)) %>% filter(name %in% Maas.pd.pc.sig$feature)
  colnames(dat_pdpc.PCprev) <- c("feature", "PC")
  dat_pdpc.PREV <- left_join(dat_pdpc.PDprev, dat_pdpc.PCprev, by = "feature") %>% melt()
  dat_pdpc.PREV$feature <- factor(dat_pdpc.PREV$feature, levels = PDovrPC.order$feature)
  dat_pdpc.PREV$variable <- factor(dat_pdpc.PREV$variable, levels = c("PC", "PD"))
  g3 <- prevalence_barplot(dat_pdpc.PREV, cols.pdpc, alfa = 0.2)
  
  
  #--------------------------------------------------------------------------
  #                              PD v HC Plots 
  #--------------------------------------------------------------------------
  
  # Subset of phyloseq obj subset to get samples of interest
  dat_pdhc.PD = subset_samples(dat_pdhc, donor_group =="PD")
  dat_pdhc.PDabun <- as.data.frame.matrix(abundances(dat_pdhc.PD))  %>% 
    rownames_to_column() %>%  filter(rowname %in% Maas.pd.hc.sig$feature) %>% column_to_rownames()
  
  dat_pdhc.HC = subset_samples(dat_pdhc, donor_group =="HC")
  dat_pdhc.HCabun <- as.data.frame.matrix(abundances(dat_pdhc.HC)) %>% 
    rownames_to_column() %>%  filter(rowname %in% Maas.pd.hc.sig$feature) %>% column_to_rownames()
  
  
  ######  Generalized or pseudo-fold change 
  gfc_data <- generalized_fold_change(pd_abundance=dat_pdhc.PDabun, 
                                      ctrl_abundances=dat_pdhc.HCabun)
  PDovrHC <- tibble("feature" = rownames(dat_pdhc.PDabun), "gFC" = gfc_data)
  
  
  ###### To get order for Significant feauture plots
  ### Stop chunk after next line to order by gFC
  PDovrHC.order <- PDovrHC[order(-PDovrHC$gFC),]
  ### Use for ordering by qval
  PDovrHC.order <- left_join(PDovrHC.order, Maas.pd.hc.sig, by = "feature")
  PDovrHC.order <- PDovrHC.order[order(-PDovrHC.order$qval),]
  
  
  ###### Generalized Fold Change (gFC) BarPlot ######
  PDovrHC.BP <- PDovrHC.order
  PDovrHC.BP <- mutate(PDovrHC.BP, direction = if_else(PDovrHC.BP$gFC > 0, "PD",
                                                       if_else(PDovrHC.BP$gFC < 0, "HC",  "error")))
  PDovrHC.BP$feature <- factor(PDovrHC.BP$feature, levels = PDovrHC.order$feature) 
  h0 <- gfc_plot(PDovrHC.BP, cols.pdhc, alfa = 0.8)
  
  
  ###### Significance Plot ###### 
  sigplot.df.pdhc <- dplyr::select(Maas.pd.hc.sig, c("feature", "pval", "qval")) %>% melt()
  sigplot.df.pdhc$feature <- factor(sigplot.df.pdhc$feature, levels = PDovrHC.order$feature) 
  h2 <- significance_barplot(sigplot.df.pdhc)
  
  
  ###### BoxPlots ######
  ## Prepping Significance labels
  abun.pdhc.inpt <- daf_boxplot_sigvalues(sigplot.df.pdhc, abun.pdhc.inpt)
  abun.pdhc.inpt$Var2 <- factor(abun.pdhc.inpt$Var2, levels = PDovrHC.order$feature)
  h1 <- daf_boxplots(abun.pdhc.inpt, fill_cols = cols.pdhc, rim_cols = cols.pdhc.rim, alfa = 0.2, obj.name = obj.name)
  
  
  ###### Prevalence Plot ######
  # Subset of phyloseq obj subset to get samples of interest
  dat_pdhc.PDprev <- tibble::enframe(prevalence(dat_pdhc.PD)) %>% filter(name %in% Maas.pd.hc.sig$feature) 
  colnames(dat_pdhc.PDprev) <- c("feature", "PD")
  dat_pdhc.HCprev <- tibble::enframe(prevalence(dat_pdhc.HC)) %>% filter(name %in% Maas.pd.hc.sig$feature)
  colnames(dat_pdhc.HCprev) <- c("feature", "HC")
  dat_pdhc.PREV <- left_join(dat_pdhc.PDprev, dat_pdhc.HCprev, by = "feature") %>% melt()
  dat_pdhc.PREV$feature <- factor(dat_pdhc.PREV$feature, levels = PDovrHC.order$feature)
  dat_pdhc.PREV$variable <- factor(dat_pdhc.PREV$variable, levels = c("HC", "PD"))
  h3 <- prevalence_barplot(dat_pdhc.PREV, cols.pdhc, alfa = 0.2)
  
  #--------------------------------------------------------------------------
  #                            Merge figures
  #--------------------------------------------------------------------------
  
  h0a <- h0 + theme(axis.text.y = element_blank())
  h1a <- h1 + theme(legend.position = "none", axis.text.y = element_text(face = "italic"))
  h2a <- h2 + theme(axis.text.y = element_blank(), legend.position = c(.90, .85))
  h3a <- h3 + theme(axis.text.y = element_blank(), legend.position = "none")
  
  g0a <- g0 + theme(axis.text.y = element_blank())
  g1a <- g1 + theme(legend.position = "none", axis.text.y = element_text(face = "italic"))
  g2a <- g2 + theme(axis.text.y = element_blank(), legend.position = c(.90, .35))
  g3a <- g3 + theme(axis.text.y = element_blank(), legend.position = "none")
  
  g0a <- g0a + theme(axis.title.x = element_blank())
  g1a <- g1a + theme(axis.title.x = element_blank())
  g2a <- g2a + theme(axis.title.x = element_blank())
  g3a <- g3a + theme(axis.title.x = element_blank())
  
  # Setting plot length variables - 
  top_len <- length(unique(PDovrPC.BP$feature)) + 2
  bottom_len <- length(unique(PDovrHC.BP$feature)) + 2.5
  DAF_part1 <- cowplot::plot_grid(g1a, h1a, nrow = 2, align = "hv", labels = "AUTO",
                                  rel_heights = c(top_len, bottom_len))
  DAF_part2 <- cowplot::plot_grid(g2a, g3a, g0a, h2a, h3a, h0a, nrow = 2, ncol=3, align = "h", 
                                  rel_heights = c(top_len, bottom_len),
                                  rel_widths = c(1, 1, 1.5))
  DAF_final <- cowplot::plot_grid(DAF_part1, DAF_part2, ncol = 2, rel_widths = c(1.5, 1))
  print(DAF_final)
  
  ggsave(DAF_final, filename = paste0("data/DAF_Analysis/", cohort, "/DAF_", obj.name, tag, ".svg"),
         width = 20, height = (top_len+bottom_len)/3)
}


#-------------------------------------------------------------------------------
#                   DAF fold change summary plot
#-------------------------------------------------------------------------------

plot_daf_summary <- function(obj.name, obj, cohort = "TBC", tag = "", 
                             repelyup = 0, repelydn = 0){
  
  # # TROUBLESHOOTING
  # obj.name = "KOs.slim"
  # obj = dat.KOs.slim 
  # cohort = "Merged"
  # repelyup = 0.0025
  # repelydn = 0.005
  
  
  # Normalization and Transformation
  dat_obj <- obj %>% 
    subset_samples(donor_id %ni% low_qc[[1]]) %>% 
    microbiome::transform("compositional")
  otu_table(dat_obj) <- asin(sqrt(otu_table(dat_obj)))
  
  # PD v PC
  dat_pdpc = subset_samples(dat_obj, donor_group !="HC")
  abun.pdpc <- as.data.frame.matrix(abundances(dat_pdpc))
  # PD v HC
  dat_pdhc = subset_samples(dat_obj, paired !="No")
  abun.pdhc <- as.data.frame.matrix(abundances(dat_pdhc))
  
  ### Read-in MaAsLin2 output
  Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/", cohort, "/", 
                                obj.name, "_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>%
    filter(value == "Population Control")
  Maas.pd.pc.sig <- Maas.pd.pc %>%
    filter(qval < 0.25)
  
  Maas.pd.hc <- read_tsv(paste0("data/MaAsLin2_Analysis/",  cohort, "/", 
                                obj.name, "_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>%
    filter(value == "Household Control")
  Maas.pd.hc.sig <- Maas.pd.hc %>%
    filter(qval < 0.25)
  
  # Join MaAsLin stats to select most significant features
  q.stats <- left_join(Maas.pd.pc, Maas.pd.hc, by = "feature") %>% 
    dplyr::mutate(q.average = (qval.x + qval.y)/2 ) %>% 
    select(feature, q.average)
  
  #--------------------------------------------------------------------------
  #                              PD v PC gFC
  #--------------------------------------------------------------------------
  
  # # Subset of phyloseq obj subset to get samples of interest
  # dat_pdpc.PD = subset_samples(dat_pdpc, donor_group =="PD")
  # dat_pdpc.PDabun <- as.data.frame.matrix(abundances(dat_pdpc.PD))
  # dat_pdpc.PC = subset_samples(dat_pdpc, donor_group =="PC")
  # dat_pdpc.PCabun <- as.data.frame.matrix(abundances(dat_pdpc.PC))
  # 
  # ######  Generalized or pseudo-fold calculation  ######
  # gfc_data <- generalized_fold_change(dat_pdpc.PDabun, dat_pdpc.PCabun)
  # PDovrPC <- tibble("feature" = rownames(dat_pdpc.PCabun), "gFC_PDvsPC" = gfc_data)
  # 
  pdpc.stats1 <- 
    dat_pdpc %>% 
    abundances() %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "id")
  PDovrPC <- 
    pdpc.stats1 %>% 
    group_col_from_ids(ids = pdpc.stats1$id) %>% 
    gather(-c(id, group), key=feature, value=Abundance) %>% 
    dplyr::group_by(group, feature) %>% 
    dplyr::summarise(mean = mean(Abundance), median = median(Abundance)) %>% 
    pivot_wider(names_from = group, values_from = c(mean, median)) %>% 
    mutate(pc_median_l2fc = log2( (median_PD + 1) / (median_PC + 1) )) %>% 
    mutate(pc_mean_l2fc = log2( (mean_PD + 1) / (mean_PC + 1) ))
  
  #--------------------------------------------------------------------------
  #                              PD v HC gFC
  #--------------------------------------------------------------------------
  
  # # Subset of phyloseq obj subset to get samples of interest
  # dat_pdhc.PD = subset_samples(dat_pdhc, donor_group =="PD")
  # dat_pdhc.PDabun <- as.data.frame.matrix(abundances(dat_pdhc.PD))
  # dat_pdhc.HC = subset_samples(dat_pdhc, donor_group =="HC")
  # dat_pdhc.HCabun <- as.data.frame.matrix(abundances(dat_pdhc.HC)) 
  # 
  # ######  Generalized or pseudo-fold change
  # gfc_data <- generalized_fold_change(pd_abundance=dat_pdhc.PDabun,
  #                                     ctrl_abundances=dat_pdhc.HCabun)
  # PDovrHC <- tibble("feature" = rownames(dat_pdhc.PDabun), "gFC_PDvsHC" = gfc_data)
  
  pdhc.stats1 <- 
    dat_pdhc %>% 
    abundances() %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "id")
  PDovrHC <- 
    pdhc.stats1 %>% 
    group_col_from_ids(ids = pdhc.stats1$id) %>% 
    gather(-c(id, group), key=feature, value=Abundance) %>% 
    dplyr::group_by(group, feature) %>% 
    dplyr::summarise(mean = mean(Abundance), median = median(Abundance)) %>% 
    pivot_wider(names_from = group, values_from = c(mean, median)) %>% 
    mutate(hc_median_l2fc = log2( (median_PD + 1) / (median_HC + 1) )) %>% 
    mutate(hc_mean_l2fc = log2( (mean_PD + 1) / (mean_HC + 1) ))
  
  
  #--------------------------------------------
  #           Plots
  #--------------------------------------------
  
  # Table for shared depleted or enriched features labels
  features <- taxa_names(obj)
  features <- tibble("Node" = features)
  
  obj_cols <- daf_group_colors(features = features, 
                               Maas.pd.pc.sig = Maas.pd.pc.sig, 
                               Maas.pd.hc.sig = Maas.pd.hc.sig) %>% 
    dplyr::rename(feature = Node)
  
  obj.plot <- left_join(PDovrHC, PDovrPC,  by = "feature") %>% 
    left_join(obj_cols, by = "feature") %>% 
    left_join(q.stats, by = "feature")
  
  obj.plot.label <- 
    obj.plot %>% 
    dplyr::filter(group != "None") %>% 
    top_n(10, wt=-q.average) %>% 
    decode_rfriendly_rows(passed_column = "feature")
  
  pal.color <-
    c(
      "BOTH" = "#2ca02c",
      "PD" = "#bfbfbf",
      "PC" = "#ed7d31",
      "HC" = "#5b9bd5"
    )
  
  obj.plot.sig <- obj.plot %>% 
    dplyr::filter(group != "None")
  obj.plot.nonsig <- obj.plot %>% 
    dplyr::filter(group == "None")
    
  gFC_diff_plot <- 
    ggplot() +
    geom_point(data = obj.plot.nonsig,
               aes(x = hc_mean_l2fc, y = pc_mean_l2fc), 
               alpha = 0.2, color = "grey") +
    geom_point(data = obj.plot.sig,
               aes(x = hc_mean_l2fc, y = pc_mean_l2fc, color = group), 
               alpha = 0.7) +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, linetype = 3, color = "darkgrey") +
    geom_abline(intercept = 0, slope = -1, linetype = 3, color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = 1, color = "grey") +
    geom_hline(yintercept = 0, linetype = 1, color = "grey") +
    labs(x = expression(paste("log"[10] * '(PD/HC)')),
         y = expression(paste("log"[10] * '(PD/PC)')),
         title = obj.name,
         color = "Feature \nAssociation") +
    scale_color_manual(values = pal.color) +
    geom_text_repel(data = obj.plot.label,
                    aes(x = hc_mean_l2fc, y = pc_mean_l2fc, label = fullnames), 
                    segment.alpha = 0.5,
                    segment.size = 0.2, 
                    size = 1.75,
                    force = 1,
                    max.time = 1,
                    max.iter = Inf,
                    nudge_y = ifelse(obj.plot.label$pc_mean_l2fc > 0, repelyup, 
                                     ifelse(obj.plot.label$pc_mean_l2fc < 0, -repelydn, 0)),
                    max.overlaps = Inf) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  print(gFC_diff_plot)

  ggsave(gFC_diff_plot, filename = paste0("data/DAF_Analysis/", cohort, "/DAF_", obj.name, tag, ".png"),
         width = 6, height = 5)
}





