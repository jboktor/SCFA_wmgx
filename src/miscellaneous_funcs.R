### Miscellaneous Functions


#-----------------------------------------------------------------------------------------------------------
# Aesthetic variables

## Color Schemes
cols.pdpc <- c("PD"= "#bfbfbf", "PC" = "#ed7d31")
cols.pdhc <- c("PD"= "#bfbfbf", "HC" = "#5b9bd5")
cols.pdpchc <- c("PD"= "#bfbfbf", "PC" = "#ed7d31", "HC" = "#5b9bd5")
cols.pdpchc.dark <- c("PD"= "#494949", "PC" = "#ed7d31", "HC" = "#5b9bd5")

# Rims
cols.pdpc.rim <- c("PD"= "#494949", "PC" = "#c15811")
cols.pdhc.rim <- c("PD"= "#494949", "HC" = "#2e75b5")
cols.pdpchc.rim <- c("PD"= "#494949", "PC" = "#c15811", "HC" = "#2e75b5")

load_alphadiv_colors <- function(){ 
  cols.pdpchc <- c("PD"= "#494949", 
                   "PC" = "#ed7d31",
                   "HC" = "#5b9bd5")
  cols.pdpchc.rim <- c("PD"= "#494949", 
                       "PC" = "#c15811",
                       "HC" = "#2e75b5")
  assign("cols.pdpchc", cols.pdpchc, envir = .GlobalEnv)
  assign("cols.pdpchc.rim", cols.pdpchc.rim, envir = .GlobalEnv)
  
  }
  
load_betadiv_colors <- function(){
  cols.pdpchc <- c("PD Patient"= "#bfbfbf", 
                   "Population Control" = "#ed7d31",
                   "Household Control" = "#5b9bd5")
  cols.pdpchc.dark <- c("PD Patient"= "#494949", 
                        "Population Control" = "#ed7d31",
                        "Household Control" = "#5b9bd5")
  cols.pdpchc.rim <- c("PD Patient"= "#494949", 
                       "Population Control" = "#c15811",
                       "Household Control" = "#2e75b5") 
  
  assign("cols.pdpchc", cols.pdpchc, envir = .GlobalEnv)
  assign("cols.pdpchc.dark", cols.pdpchc.dark, envir = .GlobalEnv)
  assign("cols.pdpchc.rim", cols.pdpchc.rim, envir = .GlobalEnv)
}

#-----------------------------------------------------------------------------------------------------------
 # Remove Objects from environment

remove_dats <- function() {
  obj <- 
    c("dat.species",
      "dat.path",
      "dat.ec",
      "dat.KOs",
      "dat.EGGNOGs",
      "dat.PFAMs",
      "dat.path.slim",
      "dat.ec.slim",
      "dat.KOs.slim",
      "dat.EGGNOGs.slim",
      "dat.PFAMs.slim")
  withCallingHandlers(rm(list = obj), warning=function(w){invokeRestart("muffleWarning")})
  }


#--------------------------------------------------------------------------------------------------
# Load number of clean sample reads
#--------------------------------------------------------------------------------------------------

# load_reads <- function(){
#   func_reads <- read_tsv("files/humann2_read_and_species_count_table.tsv", col_names = T)
#   reads <- dplyr::select(func_reads, c("# samples","total reads")) %>% 
#     dplyr::rename( "id" = "# samples", "clean_total_reads" = "total reads")
#   return(reads)
# }

load_reads <- function(cohort){
  
  negative_controls <- c("S00A4-ATCC_MSA_1003_S96", 
                         "S00A4-neg2_S119", 
                         "S00A4-neg3_S125",
                         "S00A4-neg_S118", 
                         "S00A4NegExt_P00A4_S94", 
                         "S00A4NegH2O_P00A4_S95",
                         "S00A4_stagPos_S117", 
                         "BLANK")
  
  TBC_keys <- read.csv(file = "files/metadata_keys.csv", header= TRUE) %>% 
    dplyr::select(c(MBI_Sample_ID, id)) %>% 
    mutate(id = gsub("_", ".", id)) %>% 
    mutate(MBI_Sample_ID = as.character(MBI_Sample_ID)) %>% 
    mutate(id = as.character(id)) %>% 
    dplyr::rename(`# samples` = MBI_Sample_ID)
  
  RUSH_keys <- read.csv(file = "files/metadata_phyloseq_RUSH.csv", header= TRUE) %>% 
    dplyr::filter(study_group == "PD") %>% 
    dplyr::select(donor_id, host_subject_id) %>% 
    dplyr::mutate(donor_id = as.character(donor_id)) %>% 
    dplyr::mutate(host_subject_id = as.character(host_subject_id)) %>% 
    dplyr::rename(`# samples` = host_subject_id)
  
  if (cohort == "TBC"){
    func_reads_TBC <-
      read_tsv(
        "files/TBC_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
        col_names = T)
    reads <- 
      func_reads_TBC %>% 
      dplyr::filter(`# samples` %ni%  negative_controls) %>% 
      dplyr::mutate(`# samples` = substr(`# samples`, 1, 10)) %>% 
      left_join(TBC_keys, by = "# samples") %>% 
      dplyr::select(c("id","total reads")) %>% 
      dplyr::rename("donor_id" = "id", "clean_total_reads" = "total reads") %>% 
      return(reads)
    
  } else if (cohort == "RUSH") {
    
    func_reads_RUSH <-
      read_tsv(
        "files/RUSH_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
        col_names = T)
    reads <- 
      func_reads_RUSH %>% 
      dplyr::filter(str_detect(`# samples`, "BLANK", negate = TRUE)) %>% 
      dplyr::filter(str_detect(`# samples`, "MSA", negate = TRUE)) %>% 
      dplyr::mutate(`# samples` = substr(`# samples`, 1, 6)) %>% 
      left_join(RUSH_keys, by = "# samples") %>% 
      dplyr::select(c("donor_id", "total reads")) %>%
      dplyr::rename("clean_total_reads" = "total reads")
    return(reads)
    
  } else if (cohort == "Merged") {
    
    func_reads_TBC <-
      read_tsv(
        "files/TBC_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
        col_names = T
      )
    reads_TBC <-
      func_reads_TBC %>%
      dplyr::filter(`# samples` %ni%  negative_controls) %>%
      dplyr::mutate(`# samples` = substr(`# samples`, 1, 10)) %>%
      left_join(TBC_keys, by = "# samples") %>%
      dplyr::select(c("id", "total reads")) %>%
      dplyr::rename("donor_id" = "id", "clean_total_reads" = "total reads") %>% 
      dplyr::mutate(cohort = "TBC")
    func_reads_RUSH <-
      read_tsv(
        "files/RUSH_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
        col_names = T
      )
    reads_RUSH <-
      func_reads_RUSH %>%
      dplyr::filter(str_detect(`# samples`, "BLANK", negate = TRUE)) %>%
      dplyr::filter(str_detect(`# samples`, "MSA", negate = TRUE)) %>%
      dplyr::mutate(`# samples` = substr(`# samples`, 1, 6)) %>%
      left_join(RUSH_keys, by = "# samples") %>%
      dplyr::select(c("donor_id", "total reads")) %>%
      dplyr::rename("clean_total_reads" = "total reads") %>% 
      dplyr::mutate(cohort = "RUSH")
    reads <- rbind(reads_TBC, reads_RUSH)
    
    return(reads)
    
  }
}
#-----------------------------------------------------------------------------------------------------------
######## p-value significance (integer to symbol function)

sig_mapper <- function(pval, shh = F, porq = "p", symbols = T) {
  ###' Traditional mapping of p-value to symbol 
  ###' prints p-values if below significance
  
  if (symbols == T){
    if (is.na(pval)){
      sigvalue = ""
    } else if (pval <= .001) {
      sigvalue = "***"
    } else if (pval <= .01) {
      sigvalue = "**"
    } else if (pval <= .05) {
      sigvalue = "*"
    } else if (pval > .05 & shh == F) {
      sigvalue = paste0(porq, "=", format.pval(pval, digits=2)) 
    } else if (pval > .05 & shh == T) {
      sigvalue = ""
    }
  } else if (symbols == F){
    sigvalue = paste0(porq, "=", format.pval(pval, digits=2)) 
  }
  return(sigvalue)
}

#-----------------------------------------------------------------------------------------------------------

sig.symbol.generator <- function(Column, porq = "p", shh = F){
  sig.symbol <- c()
  for (i in Column){
    sig.symbol <- c(sig.symbol, sig_mapper(i, porq = porq, shh = shh))
  }
  return(sig.symbol)
}

#-----------------------------------------------------------------------------------------------------------

distribution_sanity <- function(df) {
  
  ###' input an abundance table and displays
  ###' a histogram and ECDF distribution of data 
  ###' colored by donor group
  
  abund.melt <- melt(df)
  abund.melt <- 
    mutate(abund.melt, group = if_else(grepl("HC", variable), "HC",
                                if_else(grepl("PC", variable), "PC","PD")))
  
  histo_plot <- ggplot(abund.melt, aes(x=value, fill = group), alpha = 0.4) + 
    theme_minimal() +
    geom_histogram(color = "black", position="dodge", boundary = 0) +
    theme(axis.title.x = element_blank(),
          legend.position = c(0.9, 0.5))
  
  ecdf_plot <- ggplot(abund.melt, aes(x=value, colour = group)) + stat_ecdf(geom = "step", pad = FALSE) +
    theme_minimal() +
    labs(y = "ECDF") +
    theme(legend.position = c(0.9, 0.5))
  
  cowplot::plot_grid(histo_plot, ecdf_plot, ncol = 1, align="v")
}

#--------------------------------- Version 2

distribution_sanity2 <- function(df, binN = 30) {
  
  ###' input an abundance table and displays
  ###' a histogram and ECDF distribution of data 
  ###' colored by donor group
  
  abund.melt <- melt(df)
  abund.melt <- 
    mutate(abund.melt, group = if_else(grepl("HC", Var2), "HC",
                                       if_else(grepl("PC", Var2), "PC","PD")))
  cols=c("PC"= "#bfbfbf", 
         "PD" = "#ed7d31", 
         "HC" = "#5b9bd5")
  
  histo_plot <- ggplot(abund.melt, aes(x=value, fill = group), alpha = 0.4) + 
    theme_minimal() +
    geom_histogram(color = "black", position="dodge", boundary = 0, bins = binN) +
    scale_fill_manual(values = cols) +
    theme(axis.title.x = element_blank(),
          legend.position = c(0.9, 0.5)) 
  
  ecdf_plot <- ggplot(abund.melt, aes(x=value, colour = group)) + stat_ecdf(geom = "step", pad = FALSE) +
    theme_minimal() +
    labs(y = "ECDF") +
    scale_colour_manual(values = cols) +
    theme(legend.position ="none")
  
  cowplot::plot_grid(histo_plot, ecdf_plot, ncol = 1, align="v")
}


#-----------------------------------------------------------------------------------------------------------
################################  Functions to adjust feature names  ########################################
#-----------------------------------------------------------------------------------------------------------


prep_species_names <- function(phylo_obj){
  taxa_names(phylo_obj) <- gsub("s__", "", taxa_names(phylo_obj))
  return(features)
}

prep_pathway_names <- function(phylo_obj){
  features <- paste0("PATHWAY_", taxa_names(phylo_obj))
  features <- gsub(":", ".", features)
  features <- gsub("\\|", ".", features)
  features <- gsub(" ", "_", features)
  features <- gsub("-", "_", features)
  return(features)
}
prep_enzyme_names <- function(phylo_obj){
  features <- paste0("ENZYME_", taxa_names(phylo_obj))
  features <- gsub(":", ".", features)
  features <- gsub("\\|", ".", features)
  features <- gsub(" ", "_", features)
  features <- gsub("-", "_", features)
  return(features)
}
prep_ko_names <- function(phylo_obj){
  features <- taxa_names(phylo_obj)
  features <- gsub(":", ".", features)
  features <- gsub("\\|", ".", features)
  features <- gsub(" ", "_", features)
  features <- gsub("-", "_", features)
  return(features)
}

#-----------------------------------------------------------------------------------------------------------

group_col_from_ids <- function(df.in, ids){
  df.out <- dplyr::mutate(df.in, group = if_else(grepl("HC", ids), "HC",
                                   if_else(grepl("PC", ids), "PC","PD")))
  # rownames(df.out) <- rownames(df.in)
  return(df.out)
}

group_col_from_ids2 <- function(df.in){
  df.out <- 
    df.in %>% 
    dplyr::mutate(group = if_else(grepl("HC", donor_id), "HC",
                                if_else(grepl("PC", donor_id), "PC", "PD")))
  
  rownames(df.out) <- rownames(df.in)
  return(df.out)
}


#-----------------------------------------------------------------------------------------------------------

boxplot_all <- function(df, x, y, cols = group.cols, title = blank.title, ylabel = blank.ylabel){
  
  blank.title = " "; blank.ylabel = " "
  group.cols = c("PC"= "#ed7d31", "PD" = "#bfbfbf", "HC" = "#5b9bd5")
  ###' Basic all group boxplot function
  
  set.seed(123)
  ggplot(data=df, aes(x=x, y=y)) +
    geom_boxplot(aes(color = x), outlier.alpha = 0, width = 0.9) +
    geom_point(aes(fill = x), position = position_jitterdodge(jitter.width = 1), 
               shape=21, size=1.5, alpha = 0.8) +
    theme_classic() +
    ggtitle(title) +
    labs(y = ylabel) +
    scale_color_manual(values = cols, name ="Group") +
    scale_fill_manual(values = cols, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank())
  
}

boxplot_all_facet <- function(df, x, y, cols = group.cols, title = blank.title, ylabel = blank.ylabel){
  
  blank.title = " "; blank.ylabel = " "
  group.cols = c("PC"= "#ed7d31", "PD" = "#bfbfbf", "HC" = "#5b9bd5")
  ###' Basic all group boxplot function with a cohort facet wrap
  
  set.seed(123)
  plot <- 
    ggplot(data=df, aes(x=x, y=y)) +
    geom_boxplot(aes(color = x), outlier.alpha = 0, width = 0.9) +
    geom_point(aes(fill = x), position = position_jitterdodge(jitter.width = 1), 
               shape=21, size=1.5, alpha = 0.8) +
    theme_classic() +
    ggtitle(title) +
    labs(y = ylabel) +
    facet_wrap(~cohort) +
    scale_color_manual(values = cols, name ="Group") +
    scale_fill_manual(values = cols, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank())
  
}

#-----------------------------------------------------------------------------------------------------------

boxplot_all_generic <- function(df, x, y, cols = group.cols, title = blank.title, ylabel = blank.ylabel){
  
  blank.title = " "; blank.ylabel = " "
  ###' Basic all group boxplot function
  
  set.seed(123)
  ggplot(data=df, aes(x=x, y=y)) +
    geom_boxplot(aes(color = x), outlier.alpha = 0, width = 0.9) +
    geom_point(aes(fill = x), position = position_jitterdodge(jitter.width = 1), 
               shape=21, size=1.5, alpha = 0.8) +
    theme_classic() +
    ggtitle(title) +
    labs(y = ylabel, color = "Group") +
    guides(fill = FALSE) +
    scale_color_d3()+
    scale_fill_d3() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank())
  
}

#-----------------------------------------------------------------------------------------------------------

alpha_div_boxplots <- function(df, x, y, 
                               df.pairs, df.pairs.x, df.pairs.y, pairs.column, 
                               cols, cols.rim, ylabel,PDvPC.stat, PDvHC.stat){
  
  ###' Function for alpha diveristy
  ###' boxplots with paired lines
  
  set.seed(123)
  p <- ggplot(data = df, aes(x = x, y = y)) + 
    theme_minimal() + 
    geom_point(aes(fill = x, color = x), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 1) +
    geom_boxplot(aes(fill = x, color = x), width=0.3, alpha = 0.1, outlier.alpha = 0) +
    geom_line(data = df.pairs, aes(x = df.pairs.x, y = df.pairs.y, group = pairs.column), 
              linetype = 'solid', color = "grey", alpha = 0.7) +
    theme(axis.title.x=element_blank(),
          legend.position = "none") +
    ylab(ylabel) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    labs(fill="Group") +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols.rim) +
    geom_signif(comparisons = list(c("PD", "HC")), annotations = sig_mapper(PDvHC.stat)) +
    geom_signif(comparisons = list(c("PC", "PD")), annotations = sig_mapper(PDvPC.stat))
  return(p)
}


#-----------------------------------------------------------------------------------------------------------
#                           Functions to explore select features in dataset
#-----------------------------------------------------------------------------------------------------------

explore_table <- function(obj){
  data.table <- obj %>%
    microbiome::abundances() %>% 
    as.data.frame() %>% 
    rownames_to_column() 
  return(data.table)
}

plot_feature <- function(obj, feature){
  
  d <- obj %>% abundances() %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    filter(rowname == feature) %>% 
    melt()
  
  dm <-  group_col_from_ids(d, id=d$variable)
  
  dm$value <- asin(sqrt(dm$value))
  dm$group <- factor(dm$group, levels = c("PC", "PD", "HC"))
  
  boxplot_all(dm, x=dm$group, y=dm$value, 
              cols=c("PC"= "#bfbfbf", 
                     "PD" = "#ed7d31", 
                     "HC" = "#5b9bd5"), 
              title=" ", 
              ylabel= paste(unique(dm$rowname), "Abundance"))
}

#-----------------------------------------------------------------------------------------------------------
#                                            ML FUNCTIONS
#-----------------------------------------------------------------------------------------------------------
# from: http://jaehyeon-kim.github.io/2015/05/Setup-Random-Seeds-on-Caret-Package.html

setSeeds <- function(method = "repeatedcv", numbers = 1, repeats = 1, tunes = NULL, seed = 42) {
  
  #' function to set up random seeds for repeated ML cross validation models
  
  #B is the number of resamples and integer vector of M (numbers + tune length if any)
  
  B <- if (method == "cv") numbers
  else if(method == "repeatedcv") numbers * repeats
  else NULL
  
  if(is.null(length)) {
    seeds <- NULL
  } else {
    set.seed(seed = seed)
    seeds <- vector(mode = "list", length = B)
    seeds <- lapply(seeds, function(x) sample.int(n = 1000000, size = numbers + ifelse(is.null(tunes), 0, tunes)))
    seeds[[length(seeds) + 1]] <- sample.int(n = 1000000, size = 1)
  }
  # return seeds
  seeds
}

#-----------------------------------------------------------------------------------------------------------

unregister <- function() {
  #' Unregister a foreach backend 
  #' Use when looping functions that run 
  #' analyses in parallel
  
  enviorn <- foreach:::.foreachGlobals
  rm(list=ls(name=enviorn), pos=enviorn)
}


#-----------------------------------------------------------------------------------------------------------

group_col_from_ids_ML <- function(df.in, ids){
  df.out <- mutate(df.in, group = if_else(grepl("control", ids), "control", "disease"))
  rownames(df.out) <- rownames(df.in)
  return(df.out)
}

#-----------------------------------------------------------------------------------------------------------


prep.CMD.Species.ML <- function(study, metafilter = NA){
  
  #' Funtion that preps input for ML analysis 
  #' Input: Study name from curatedMetagenomicData
  #' Output: dataframe of Species abundnace inluding an 
  #' identifying group column
  
  # alt.disease <- curatedMetagenomicData(paste0(study, ".metaphlan_bugs_list.stool"), dryrun=F)
  # 
  # df.spec <- alt.disease[[1]] %>%
  #   ExpressionSet2phyloseq() %>% 
  #   subset_taxa(!is.na(Species)) %>% 
  #   subset_taxa(is.na(Strain)) 
  
  # study <- study

  df.abund <- study %>% 
    # microbiome::transform("compositional") %>%
    microbiome::abundances() %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column()
  # df.abund[-1] <- asin(sqrt(df.abund[-1]))
  
  m <- microbiome::meta(study) %>% 
    select(study_condition) %>% 
    rownames_to_column()
  
  model.input <- left_join(df.abund, m) %>% 
    column_to_rownames() 
  
  if (!is.na(metafilter)){
    model.input <- filter(model.input, study_condition != metafilter)
  }
  
  count(model.input$study_condition)
  unique(model.input$study_condition)
  
  model.input <- model.input %>% 
    group_col_from_ids_ML(ids = model.input$study_condition) %>% 
    dplyr::select(-study_condition)
  
  model.input$group <- factor(model.input$group)
  
  return(model.input)
}

#-----------------------------------------------------------------------------------------------------------

# helper function for vizualizing Xgboost plots
# https://www.kaggle.com/pelkoja/visual-xgboost-tuning-with-caret/report
tuneplot <- function(x, probs = .90) {
  ggplot(x) +
    theme_bw()
}

#-----------------------------------------------------------------------------------------------------------


# Plot feature of interest (foi)
plot_foi <- function(datObj, foi){
  
  asin.input <- datObj %>%
    microbiome::transform("compositional") %>% 
    microbiome::abundances()
  df1 <- asin(sqrt(asin.input)) %>% 
    as.data.frame() %>% 
    rownames_to_column() 
  d <- df1 %>% 
    filter(rowname == foi) %>% 
    melt()
  dm <-  group_col_from_ids(d, id=d$variable)
  dm$group <- factor(dm$group, levels = c("PC", "PD", "HC"))
  boxplot_all(dm, x=dm$group, y=dm$value, 
              title=" ", 
              ylabel= paste(unique(dm$rowname), "Abundance"))
  
}

#-----------------------------------------------------------------------------------------------------------

# Plot feature of interest (foi) for any Dataset from CMG
plot_foi_general <- function(datObj, foi, title="" ){
  
  asin.input <- datObj %>%
    microbiome::transform("compositional") %>% 
    microbiome::abundances()
  df1 <- asin(sqrt(asin.input)) %>% 
    as.data.frame() %>% 
    rownames_to_column() 
  met <- meta(df.spec) %>% 
    dplyr::select(c(study_condition)) %>% 
    rownames_to_column(var = "variable")
  d <- df1 %>% 
    filter(rowname == foi) %>% 
    melt()
  d2 <- left_join(d, met)
  d2$study_condition <- factor(d2$study_condition)
  boxplot_all_generic(d2, x=d2$study_condition, y=d2$value, 
                      title="", 
                      ylabel= paste(unique(d2$rowname), "Abundance"))
}

#-------------------------------------------------------------------


# 2) CLR Log-odds ratio : facultative anaerobic / obligate anaerobic bacteria

#' Function returns normalized values for feature(s) of interest
#' 1) Adaptive to multiple types of normalization
#' 2) Will sum values if a list of foi's are input

fois <- function(datObj, foi) {
  
  if (length(foi) > 1) {
    
    # Init
    cnt <- 1
    
    for (i in foi) {
      cat("selecting feature: ", i, "\n")
      d <- datObj %>%
        microbiome::transform("compositional") %>%
        microbiome::abundances() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        filter(rowname == i) %>%
        melt() %>%
        dplyr::select(c(variable, value))
      if (cnt == 1){
        summed.vars <- d
      } else {
        summed.vars$value <- summed.vars$value + d$value
      }
      cnt <- cnt + 1
    }
    return(summed.vars)
    
  } else {
    for (i in foi){
      cat("selecting feature: ", i, "\n")
      d <- datObj %>%
        microbiome::transform("compositional") %>%
        microbiome::abundances() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        filter(rowname == i) %>%
        melt() %>%
        dplyr::select(c(variable, value))
    }
    return(d)
  }
}

#-------------------------------------------------------------------
# Negate %in% function

'%ni%' <- Negate('%in%')

#-------------------------------------------------------------------

# SGV Stats Function

nmean_summary <- function(df){
  output.stats <- tibble()
  vars = colnames(df)[1:ncol(df)-1]
  for (sv in vars) {
    print(sv)
    stat.col <- df %>%
      dplyr::select(group, sv) %>%
      filter(!is.na(df[[sv]])) %>%
      group_by(group) %>%
      dplyr::summarise(n = n(),across(where(is.numeric),
                              ~ mean(.x, na.rm = TRUE)))
    colnames(stat.col)[3] <- "ratio_detected"

    row2add <- cbind(stat.col[1:3], sv)
    output.stats <- rbind(output.stats, row2add)
  }
  return(output.stats)
}


#-------------------------------------------------------------------

# Functions to clean up column names 

clean.cols.tax <- function(x) {
  colnames(x) <- gsub("_taxonomic_profile", "", colnames(x)); x } 

clean.cols.abund <- function(x) {
  colnames(x) <- gsub("_Abundance", "", colnames(x)); x } 

clean.cols.abund_RPK <- function(x) {
  colnames(x) <- gsub("_Abundance-RPKs", "", colnames(x)); x } 

clean.cols.abund_CPM <- function(x) {
  colnames(x) <- gsub("_Abundance-CPM", "", colnames(x)); x } 


trim_cols <- function(x, cohort) {
  if (cohort == "TBC"){
    colnames(x) <- substr(colnames(x), 1, 10); x 
  } else if (cohort == "RUSH") {
    colnames(x) <- substr(colnames(x), 1, 6); x
    }
  } 


make_rfriendly_rows <- function(df, passed_column) {
  features <- rlang::sym(passed_column)
  df %>% 
    mutate(simplenames = gsub(":", ".gc.", !!features)) %>% 
    mutate(simplenames = gsub("\\|", ".gp.", simplenames)) %>% 
    mutate(simplenames = gsub(" ", ".gs.", simplenames)) %>% 
    mutate(simplenames = gsub("-", ".gh.", simplenames)) %>% 
    mutate(simplenames = gsub("/", ".gd.", simplenames)) %>% 
    mutate(simplenames = gsub("\\]", ".gsqrr.", simplenames)) %>% 
    mutate(simplenames = gsub("\\[", ".gsqrl.", simplenames)) %>% 
    mutate(simplenames = gsub("\\)", ".gpr.", simplenames)) %>% 
    mutate(simplenames = gsub("\\(", ".gpl.", simplenames)) %>% 
    mutate(simplenames = gsub(",", ".gm.", simplenames)) %>% 
    mutate(simplenames = gsub("\\+", ".gplus.", simplenames)) %>% 
    mutate(simplenames = gsub("\\'", ".gpar.", simplenames)) %>% 
    mutate(simplenames = paste0("feat_", simplenames)) %>% 
    tibble::column_to_rownames(var = "simplenames") %>% 
    dplyr::select(-!!features)
}


decode_rfriendly_rows <- function(df, passed_column) {
  features <- rlang::sym(passed_column)
  df %>% 
    mutate(fullnames = gsub("\\.gc.", ":", !!features)) %>% 
    mutate(fullnames = gsub("\\.gsqrr.", "\\]", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gsqrl.", "\\[", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gplus.", "\\+", fullnames)) %>%
    mutate(fullnames = gsub("\\.gpar.", "\\'", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gpr.", ")", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gpl.", "(", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gm.", ",", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gp.", "\\|", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gs.", " ", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gh.", "-", fullnames)) %>% 
    mutate(fullnames = gsub("\\.gd.", "/", fullnames)) %>% 
    mutate(fullnames = gsub("feat_", "", fullnames))
}


#----------------------------------------------------------------------------------------------------------

maaslin_prep <- function(dat){
  dat <- dat %>% 
    subset_samples(donor_id %ni% low_qc[[1]]) %>% 
    core(detection = 0, prevalence = 0.1)
  return(dat)
}

#-----------------------------------------------------------------------------------------------------------


binarize <- function(x) {
  ifelse(x == 0, 0, 1)
}
#-----------------------------------------------------------------------------------------------------------
