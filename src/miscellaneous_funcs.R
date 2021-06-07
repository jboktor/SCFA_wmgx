### Miscellaneous Functions


#_______________________________________________________________________________
# Aesthetic variables

# Color Schemes
genotype.cols <- c("ASO"= "#bfbfbf", "WT" = "#ed7d31")
treatment.cols <- c("PD"= "#Prebiotic", "HC" = "#Control")
treatment.group.cols <-
  c("ASO_Prebiotic" = "#91a6cc",
    "ASO_Control" = "#f1f1f1",
    "WT_Prebiotic" = "#1D2C4E",
    "WT_Control" = "#434343")

# Rims
treatment.group.cols.rims <-
  c("ASO_Prebiotic" = "#4f70ab",
    "ASO_Control" = "#404040",
    "WT_Prebiotic" = "#1d2c4e",
    "WT_Control" = "#080808")

treatment_group_labels <- c(
  "ASO_Control" = "ASO\nControl",
  "ASO_Prebiotic" = "ASO\nPrebiotic",
  "WT_Control" = "WT\nControl",
  "WT_Prebiotic" = "WT\nPrebiotic"
)


#_______________________________________________________________________________

print_line <- function(){
  cat("---------------------------------------------------------------------\n")
}
#_______________________________________________________________________________
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


#_______________________________________________________________________________
# Load number of clean sample reads
#_______________________________________________________________________________

load_reads <- function() {
  func_reads <-
    read_tsv(
      "files/SCFA_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
      col_names = T
    )
  reads <-
    func_reads %>%
    # dplyr::filter(`# samples` %ni%  negative_controls) %>%
    dplyr::mutate(`# samples` = substr(`# samples`, 1, 5)) %>%
    dplyr::mutate(`# samples` = str_remove(`# samples`, "_$")) %>%
    dplyr::rename(donor_id = `# samples`, clean_total_reads = "total reads")
  
  return(reads)
}

#_______________________________________________________________________________
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

#_______________________________________________________________________________

sig.symbol.generator <- function(Column, porq = "p", shh = F){
  sig.symbol <- c()
  for (i in Column){
    sig.symbol <- c(sig.symbol, sig_mapper(i, porq = porq, shh = shh))
  }
  return(sig.symbol)
}

#_______________________________________________________________________________

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

#_______________________________________________________________________________

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


#_______________________________________________________________________________

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

#_______________________________________________________________________________

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

#_______________________________________________________________________________

alpha_div_boxplots <- function(df, x, y, cols, cols.rim, ylabel,tukeystats){
  
  ###' Function for plotting alpha diversity data
  
  set.seed(123)
  comps <- get_comparisons(df, treatment_group, ref.group = NULL)
  comps$V1
  p <- ggplot(data = df, aes(x = x, y = y)) + 
    theme_minimal() + 
    geom_point(aes(fill = x, color = x), position = position_jitterdodge(dodge.width=1), 
               shape=21, size=1, alpha = 1) +
    geom_boxplot(aes(fill = x, color = x), width=0.5, alpha = 0.1, outlier.alpha = 0) +
    labs(y=ylabel, fill="Group") +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols.rim) +
    scale_x_discrete(labels = treatment_group_labels) +
    theme(axis.title.x=element_blank(),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle=45, hjust=1))
  return(p)
}

#_______________________________________________________________________________
#                           Functions to explore select features in dataset
#_______________________________________________________________________________

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

#_______________________________________________________________________________


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

#_______________________________________________________________________________

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

# Functions to clean up column names 

clean.cols.tax <- function(x) {
  colnames(x) <- gsub("_taxonomic_profile", "", colnames(x)); x } 

clean.cols.abund <- function(x) {
  colnames(x) <- gsub("_Abundance", "", colnames(x)); x } 

clean.cols.abund_RPK <- function(x) {
  colnames(x) <- gsub("_Abundance-RPKs", "", colnames(x)); x } 

clean.cols.abund_CPM <- function(x) {
  colnames(x) <- gsub("_Abundance-CPM", "", colnames(x)); x } 

trim_cols <- function(x) {
  colnames(x) <- substr(colnames(x), 1, 5); x }

trim_underscore <- function(x) {
  colnames(x) <- str_remove(colnames(x), "_$"); x }

clean.cols.qiime <- function(x) {
  colnames(x) <- gsub("13244.", "", colnames(x)); x}
clean.cols.qiime2 <- function(x) {
  colnames(x) <- gsub("\\.", "\\_", colnames(x)); x} 


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


#_______________________________________________________________________________--------

maaslin_prep <- function(dat){
  dat <- dat %>% 
    core(detection = 0, prevalence = 0.1)
  return(dat)
}

#_______________________________________________________________________________

binarize <- function(x) {
  ifelse(x == 0, 0, 1)
}

#_______________________________________________________________________________

my_clean_theme <- function() {
  th <- ggplot2::theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(1, 1, 1, 2), "cm")
    )
  return(th)
}

#_______________________________________________________________________________

rarefunc <- function(dat) {
  rarefy_even_depth(dat,
                    replace = F,
                    rngseed = 42,
                    verbose = T)
}


#_______________________________________________________________________________
#####                      Correlation Functions                           ##### 
#_______________________________________________________________________________
# 
# corr_abund_prep <- function(obj,  obj.name, cohort, sigfilter = T) {
#   
#   ### Read-in MaAsLin2 output
#   Maas.pd.pc.sig <- read_tssv(
#     paste0("data/MaAsLin2_Analysis/", cohort, "/", obj.name, 
#            "_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
#     filter(value == "Population Control") %>% 
#     filter(qval < 0.25)
#   
#   Maas.pd.hc.sig <- read_tsv(
#     paste0("data/MaAsLin2_Analysis/",  cohort, "/", obj.name, 
#            "_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
#     filter(value == "Household Control") %>% 
#     filter(qval < 0.25)
#   
#   features <- full_join(Maas.pd.pc.sig, Maas.pd.hc.sig, by = "feature") %>% 
#     dplyr::select("feature")
#   
#   if(sigfilter){
#     abundance_tbl <- 
#       obj %>% 
#       subset_samples(donor_id %ni% low_qc[[1]]) %>% 
#       core(detection = 0, prevalence = 0.1) %>% 
#       abundances() %>% 
#       as.data.frame() %>% 
#       rownames_to_column() %>% 
#       filter(rowname %in% features[[1]]) %>% 
#       column_to_rownames() %>% 
#       t() %>%
#       as.data.frame()
#   } else {
#     abundance_tbl <- 
#       obj %>% 
#       subset_samples(donor_id %ni% low_qc[[1]]) %>% 
#       core(detection = 0, prevalence = 0.1) %>% 
#       abundances() %>% 
#       as.data.frame() %>% 
#       t() %>%
#       as.data.frame()
#   }
#   return(abundance_tbl)
# }


corr_loop_parallel <- function(metadata, abundance, obj.name) {
  
  # # TROUBLE
  # metadata = corr.meta
  # abundance = corr.abund
  # obj.name = obj.name
  
  require(foreach)
  require(doParallel)
  
  #setup parallel back-end to use multiple threads
  start_time <- Sys.time()
  cores = detectCores()
  cl <- makeCluster(cores[1] - 1) #not to overload your computer
  registerDoParallel(cl)
  
  
  corr_output <- tibble()
  looped_df <-
    foreach(metavar = colnames(metadata), .combine = 'rbind') %:%
    foreach(feature = colnames(abundance), .combine = 'rbind') %dopar% {
      # Calculate Spearman's Correlation
      spearman <-
        cor.test(
          x = metadata[[metavar]],
          y = abundance[[feature]],
          method = "spearman",
          na.action = na.exclude,
          alternative = "two.sided"
        )
      
      data.frame(
        "metadata" = metavar,
        "feature" = feature,
        "object_name" = obj.name,
        "rho" = spearman$estimate[[1]],
        "S" = spearman$statistic[[1]],
        "n" = length(na.omit(metadata[[metavar]])),
        "p" = spearman$p.value[[1]]
      )
    }
  stopCluster(cl)
  # Remove NAs and add FDR (Benjamini Hochberg)
  statvars <- c("rho", "S", "n", "p")
  corr_output <-
    looped_df %>%
    na.omit() %>%
    mutate(across(all_of(statvars), as.character),
           across(all_of(statvars), as.numeric)) %>%
    group_by(object_name, metadata) %>%
    mutate(q = p.adjust(p, method = 'BH')) %>%
    ungroup()
  end_time <- Sys.time()
  cat("Correlations calculated in : ", 
      end_time - start_time, attr(end_time - start_time, "units"), "\n")
  return(corr_output)
}

corr_loop <- function(metadata, abundance, obj.name) {
  corr_output <- tibble()
  for (metavar in colnames(metadata)) {
    cat("Calculating correlations for: ", metavar[[1]], "\n")
    for (feature in colnames(abundance)) {
      # Calculate Spearman's Correlation
      spearman <-
        cor.test(
          x = metadata[[metavar]],
          y = abundance[[feature]],
          method = "spearman",
          na.action = na.exclude,
          alternative = "two.sided"
        )
      row2add <-
        cbind(
          "metadata" = metavar,
          "feature" = feature,
          "object_name" = obj.name,
          "rho" = spearman$estimate[[1]],
          "S" = spearman$statistic[[1]],
          "n" = length(na.omit(metadata[[metavar]])),
          "p" = spearman$p.value[[1]]
        )
      corr_output <- rbind(corr_output, row2add)
    }
  }
  # Remove NAs and add FDR (Benjamini Hochberg)
  statvars <- c("rho", "S", "n", "p")
  corr_output <-
    corr_output %>%
    na.omit() %>%
    mutate(across(all_of(statvars), as.character),
           across(all_of(statvars), as.numeric)) %>%
    decode_rfriendly_rows(passed_column = "feature") %>%
    dplyr::select(-feature) %>%
    dplyr::rename("feature"="fullnames") %>%
    dplyr::relocate(feature, .after = metadata) %>%
    group_by(object_name, metadata) %>%
    mutate(q = p.adjust(p, method = 'BH')) %>%
    ungroup()
  return(corr_output)
}


#_______________________________________________________________________________
#####                        Correlation XY Plot                          ##### 
#_______________________________________________________________________________


corr_xy <- function(obj, corr_obj, feature_var, metadata_var){
  
  #' Function creates a scatter plot of a given feature and a metadata column
  abund <- obj %>% 
    abundances() %>% 
    as.data.frame() %>%
    rownames_to_column() %>% 
    mutate(rowname = gsub("\\-", ".", rowname)) %>% 
    mutate(rowname = gsub("\\:", ".", rowname)) %>% 
    mutate(rowname = gsub("\\+", ".", rowname)) %>% 
    column_to_rownames() %>% 
    t() %>% as.data.frame() %>%
    rownames_to_column(var = "donor_id")
  df.plot <- obj %>% 
    meta() %>% 
    process_meta() %>% 
    left_join(abund, by = "donor_id")
  
  stat_col <-
    corr_obj %>%
    dplyr::filter(feature == sym(feature_var)) %>%
    dplyr::filter(metadata == sym(metadata_var))
  stat_title <-
    paste0(
      "Spearman's Rho: ", round(stat_col$rho, digits = 3), "\n",
      "P-value: ", format(stat_col$p, digits = 3, scientific =T),
      "  FDR: ", format(stat_col$q, digits = 3, scientific =T)
    )
  cat("Rho: ", stat_col$rho, ", ")
  cat("P-value: ", stat_col$p, ",  ")
  cat("Q-value: ", stat_col$q, "\n")
  
  df.plot %>%
    drop_na(metadata_var) %>%
    ggplot(aes(x = .data[[feature_var]], y = .data[[metadata_var]])) +
    geom_point(aes(fill = treatment_group, color = treatment_group), shape=21, alpha = 1) +
    geom_smooth(method = lm, color="darkgrey", linetype="dotted", se = F) +
    theme_bw() +
    labs(x = feature_var, y = metadata_var, title = stat_title) +
    scale_fill_manual(values = treatment.group.cols) +
    scale_color_manual(values = treatment.group.cols.rims) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 12))
}

#_______________________________________________________________________________
#####                      Correlation Heatmap                           ##### 
#_______________________________________________________________________________

corr_heatmap <- function(corr.df){
  
  # Calculate number of significant (P-VALUE < 0.05) associations
  #  and select top 30 features
  corr.df.top <- corr.df %>% 
    dplyr::filter(q < 0.25) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n)) %>%
    top_n(n = 30, wt = n) %>% 
    slice_head(n = 30)
  
  # Filter correlation df for top 30 features
  corr.df.trim <- corr.df %>% 
    mutate(metadata = as.character(metadata)) %>% 
    filter(feature %in% corr.df.top$feature)
  
  # create dataframe matrix of Rho correlation values for distance functions
  rho.df <- 
    corr.df.trim %>% 
    dplyr::select(feature, metadata, rho) %>% 
    pivot_wider(names_from = feature, values_from = rho, values_fill = NA) %>% 
    column_to_rownames(var = "metadata")
  
  # Hierarchical clustering of Rho values for features & metadata
  meta.dendro <- 
    as.dendrogram(hclust(d = dist(x = rho.df), method = "complete"))
  meta.dendro.plot <- ggdendrogram(data = meta.dendro, rotate = TRUE)
  feature.dendro <- 
    as.dendrogram(hclust(d = dist(x = as.data.frame(t(rho.df))), method = "complete"))
  feature.dendro.plot <- ggdendrogram(data = feature.dendro, rotate = TRUE)
  
  ### Reorder Heatmap axis using order of dendrograms
  feature.order <- order.dendrogram(feature.dendro)
  metadata.order <- order.dendrogram(meta.dendro)
  corr.df.trim.ordered <- corr.df.trim %>% 
    dplyr::mutate(feature = factor(feature, 
                                   ordered = TRUE,
                                   levels = unique(corr.df.trim$feature)[feature.order] )) %>% 
    dplyr::mutate(metadata = factor(metadata, 
                                    ordered = TRUE,
                                    levels = unique(corr.df.trim$metadata)[metadata.order] )) 
  
  ### Plot Heatmap
  h1 <- 
    corr.df.trim.ordered %>% 
    mutate(siglabel = if_else(q < 0.20, "*", "")) %>%
    ggplot(aes(x = metadata, y = feature, fill = rho)) +
    geom_tile() + 
    geom_text(aes(label=siglabel), size=5,vjust = 0.77, color = "white") +
    labs(fill = "Spearman correlation") +
    scale_y_discrete(position = "right") +
    scale_fill_distiller(palette = "RdBu") +
    theme(axis.text.x = element_text(angle=45, hjust =1),
          legend.position = "top",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(1, 1, 1, 5), "cm"))
  
  print(h1)
  return(h1)
}
#_______________________________________________________________________________
#####                Correlation Plot helper functions                       ##### 
#_______________________________________________________________________________

trim_sig_helper <- function(df_cors, n=10){
  output <- 
    df_cors %>%
    filter(q < 0.25) %>%
    slice_min(n = n, order_by = q)
  return(output)
}

top_n_scatterplots <- function(dat, obj.name, 
                               df.cors, metaclass,
                               n_plots = 10){
  
  df.cors.sig <- trim_sig_helper(df.cors, n=n_plots)
  
  for (corr in 1:nrow(df.cors.sig)) {
    cor_row <- df.cors.sig[corr, ]
    p <- corr_xy(obj = dat, df.cors,
                 feature_var = cor_row$feature[[1]],
                 metadata_var = as.character(cor_row$metadata[[1]]) )
    
    # print(p)
    obj.name.out <- gsub('/', '_', obj.name)
    obj.name.out <- gsub(',', '', obj.name.out)
    outputdir <- paste0("data/Correlations/", obj.name.out)
    dir.create(outputdir)
    plot.name = paste0(outputdir, "/", cor_row$feature[[1]], "_VS_",
                       as.character(cor_row$metadata[[1]]), ".svg")
    ggsave(p, filename = plot.name, height = 4, width = 5)
  }
}
#_______________________________________________________________________________






