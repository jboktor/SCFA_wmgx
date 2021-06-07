# Alpha Diversity adaptable

##### Alpha Diversity Boxplots Script
######## Load Data & functions
source("src/load_packages.R")
source("src/metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/community_composition_funcs.R")
source("src/stats_funcs.R")

alpha_diversity_summary <- function(x, z, tag){
  
  # TROUBLE
  # load("files/Phyloseq_objects_Woltka/Species_PhyloseqObj.RData")
  # x <- c(dat.species)
  # z <- c("Species")
  
  # --------------------------------------------------- 
  #         Alpha Diversity Plotting Loop
  # --------------------------------------------------  
  cnt <- 1
  for (i in x){
    
    cat("Processing input: ", z[cnt], "\n")
    dat_alpha <- i
    print(dat_alpha)
    cat("\n")
    
    # Run Metadata pre-processing function
    env <- meta(dat_alpha)
    alpha_metrics <- microbiome::alpha(abundances(dat_alpha)) %>% 
      rownames_to_column(var="donor_id")
    env <- left_join(env, alpha_metrics, by="donor_id")
    comps <- get_comparisons(env, treatment_group, ref.group = NULL)
    
    # --------------------------------------------------- 
    #                  Observed Species
    # --------------------------------------------------- 
    
    ### STATS
    tukeystats <-  anova_tukey(test_col = "observed", df = env)
    p1 <- alpha_div_boxplots(df=env, x=env$treatment_group, y=env$observed,
                             cols=treatment.group.cols, cols.rim=treatment.group.cols.rims,
                             ylabel = paste0("Observed Counts: ", z[cnt]))
    
    yval <- ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range[2]
    p1 <- p1 +
      geom_signif(comparisons = list(comps$V1), 
                  annotations = sig_mapper(tukeystats$treatment_group[1, 4]), 
                  y_position = yval, tip_length = 0) +
      geom_signif(comparisons = list(comps$V2), 
                  annotations = sig_mapper(tukeystats$treatment_group[2, 4]), 
                  y_position = yval*1.1, tip_length = 0) +
      geom_signif(comparisons = list(comps$V3), 
                  annotations = sig_mapper(tukeystats$treatment_group[3, 4]), 
                  y_position = yval*1.2, tip_length = 0) +
      geom_signif(comparisons = list(comps$V4), 
                  annotations = sig_mapper(tukeystats$treatment_group[4, 4]), 
                  y_position = yval*1.05, tip_length = 0) +
      geom_signif(comparisons = list(comps$V5), 
                  annotations = sig_mapper(tukeystats$treatment_group[5, 4]), 
                  y_position = yval*1.15, tip_length = 0) +
      geom_signif(comparisons = list(comps$V6), 
                  annotations = sig_mapper(tukeystats$treatment_group[6, 4]), 
                  y_position = yval, tip_length = 0)
    p1
    
    # --------------------------------------------------- 
    #                  Shannon Diversity
    # --------------------------------------------------- 
    
    tukeystats <-  anova_tukey(test_col = "diversity_shannon", df = env)
    p2 <- alpha_div_boxplots(df=env, x=env$treatment_group, y=env$diversity_shannon, 
                             cols=treatment.group.cols, cols.rim=treatment.group.cols.rims, 
                             ylabel = paste0("Shannon's Diversity: ", z[cnt]))
    yval <- ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range[2]
    p2 <- p2 +
      geom_signif(comparisons = list(comps$V1), 
                  annotations = sig_mapper(tukeystats$treatment_group[1, 4]), 
                  y_position = yval, tip_length = 0) +
      geom_signif(comparisons = list(comps$V2), 
                  annotations = sig_mapper(tukeystats$treatment_group[2, 4]), 
                  y_position = yval*1.1, tip_length = 0) +
      geom_signif(comparisons = list(comps$V3), 
                  annotations = sig_mapper(tukeystats$treatment_group[3, 4]), 
                  y_position = yval*1.2, tip_length = 0) +
      geom_signif(comparisons = list(comps$V4), 
                  annotations = sig_mapper(tukeystats$treatment_group[4, 4]), 
                  y_position = yval*1.05, tip_length = 0) +
      geom_signif(comparisons = list(comps$V5), 
                  annotations = sig_mapper(tukeystats$treatment_group[5, 4]), 
                  y_position = yval*1.15, tip_length = 0) +
      geom_signif(comparisons = list(comps$V6), 
                  annotations = sig_mapper(tukeystats$treatment_group[6, 4]), 
                  y_position = yval, tip_length = 0)
    p2
    
    # --------------------------------------------------- 
    #                  Simpsons Evenness
    # --------------------------------------------------- 
    
    tukeystats <-  anova_tukey(test_col = "evenness_simpson", df = env)
    p3 <- alpha_div_boxplots(df=env, x=env$treatment_group, y=env$evenness_simpson, 
                             cols=treatment.group.cols, cols.rim=treatment.group.cols.rims,
                             ylabel = paste0("Simpson's Evenness: ", z[cnt]))
    yval <- ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range[2]
    # yval <- 0.025
    p3 <- p3 +
      geom_signif(comparisons = list(comps$V1), 
                  annotations = sig_mapper(tukeystats$treatment_group[1, 4]), 
                  y_position = yval, tip_length = 0) +
      geom_signif(comparisons = list(comps$V2), 
                  annotations = sig_mapper(tukeystats$treatment_group[2, 4]), 
                  y_position = yval*1.1, tip_length = 0) +
      geom_signif(comparisons = list(comps$V3), 
                  annotations = sig_mapper(tukeystats$treatment_group[3, 4]), 
                  y_position = yval*1.2, tip_length = 0) +
      geom_signif(comparisons = list(comps$V4), 
                  annotations = sig_mapper(tukeystats$treatment_group[4, 4]), 
                  y_position = yval*1.05, tip_length = 0) +
      geom_signif(comparisons = list(comps$V5), 
                  annotations = sig_mapper(tukeystats$treatment_group[5, 4]), 
                  y_position = yval*1.15, tip_length = 0) +
      geom_signif(comparisons = list(comps$V6), 
                  annotations = sig_mapper(tukeystats$treatment_group[6, 4]), 
                  y_position = yval, tip_length = 0)
    p3
    
    # --------------------------------------------------- 
    #                  Gini's Dominance 
    # --------------------------------------------------- 
    
    tukeystats <-  anova_tukey(test_col = "evenness_simpson", df = env)
    p4 <- alpha_div_boxplots(df=env, x=env$treatment_group, y=env$dominance_gini, 
                             cols=treatment.group.cols, cols.rim=treatment.group.cols.rims,
                             ylabel = paste0("Gini's Dominance: ", z[cnt]))
    yval <- ggplot_build(p4)$layout$panel_scales_y[[1]]$range$range[2]
    p4 <- p4 +
      geom_signif(comparisons = list(comps$V1), 
                  annotations = sig_mapper(tukeystats$treatment_group[1, 4]), 
                  y_position = yval, tip_length = 0) +
      geom_signif(comparisons = list(comps$V2), 
                  annotations = sig_mapper(tukeystats$treatment_group[2, 4]), 
                  y_position = yval*1.002, tip_length = 0) +
      geom_signif(comparisons = list(comps$V3), 
                  annotations = sig_mapper(tukeystats$treatment_group[3, 4]), 
                  y_position = yval*1.004, tip_length = 0) +
      geom_signif(comparisons = list(comps$V4), 
                  annotations = sig_mapper(tukeystats$treatment_group[4, 4]), 
                  y_position = yval*1.001, tip_length = 0) +
      geom_signif(comparisons = list(comps$V5), 
                  annotations = sig_mapper(tukeystats$treatment_group[5, 4]), 
                  y_position = yval*1.003, tip_length = 0) +
      geom_signif(comparisons = list(comps$V6), 
                  annotations = sig_mapper(tukeystats$treatment_group[6, 4]), 
                  y_position = yval, tip_length = 0)
    p4
    
    cat("Alpha Div Plotted \n")
    
    
    ### MERGE PLOTS ### 
    alpha_cow <- cowplot::plot_grid(p1, p2, p3, p4, nrow = 1, align = "v")
    alpha_cow
    print(alpha_cow)
    ggsave(alpha_cow, filename =
             paste0("data/Community_Composition/Alpha_Diversity/AlphaDiversity_BoxPlots_",  
                    z[cnt], "_",
                    tag,
                    "_Summary.svg"),
           height = 5, width =12)
    
    cnt <- cnt + 1
  }
}
  


  