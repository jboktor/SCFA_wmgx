### QC Functions


'%ni%' <- Negate('%in%')

#--------------------------------------------------------------------------------------------------
# Get Pseduo-Counts from clean reads and abundance table
#--------------------------------------------------------------------------------------------------

PseudoCounts <- function(dat, reads){
  
  psudocnts <- dat %>% 
    microbiome::abundances() %>% 
    as.data.frame()
  cat("TSS \n")
  
  for (i in colnames(psudocnts)){
    donor_reads <- reads[[which(reads$donor_id == i), 2]]
    psudocnts[i] <- psudocnts[i] * donor_reads
  }
  cat("Pseudocount Transformation Complete\n")
  return(psudocnts)
}


#--------------------------------------------------------------------------------------------------
# Rarefaction Analysis 
#--------------------------------------------------------------------------------------------------

feature_accumulation_plot <- function(dat, featuretype="Species", reads, cohort){

  cat("\n\n\n"); cat(featuretype, "Rarefaction  \n")

  # Get Pseuo-counts
  cat("Pseudocount Estimation \n")
  psudocnts <- dat %>% 
    subset_samples(donor_id %ni% low_qc[[1]]) %>% 
    PseudoCounts(reads)
  
  # Filter 
  psudocnts.HC <- psudocnts %>% dplyr::select(contains("HC")) %>% 
    t() %>% as.data.frame()
  psudocnts.PC <- psudocnts %>% dplyr::select(contains("PC")) %>% 
    t() %>% as.data.frame()
  psudocnts.PD <- psudocnts %>% dplyr::select(!contains(c("HC", "PC"))) %>% 
    t() %>% as.data.frame()
  
  cat("Calculating Rarefaction Estimate for HC : This may take a second -   \n")
  acc.HC <- specaccum(psudocnts.HC, method="random", permutations=1000)
  cat("Calculating Rarefaction Estimate for PC - Almost there - ༼ つ ◕_◕ ༽つ  \n")
  acc.PC <- specaccum(psudocnts.PC, method="random", permutations=1000)
  cat("Calculating Rarefaction Estimate for PD - Homestretch -  ༼ つ ಥ_ಥ ༽つ  \n\n")
  acc.PD <- specaccum(psudocnts.PD, method="random", permutations=1000)
  cat(featuretype, " Rarefaction Calculations Complete:  ヽ༼ຈل͜ຈ༽ﾉ  \n\n")
  
  df.acc.HC <- data.frame(Sites=acc.HC$sites, Richness=acc.HC$richness, SD=acc.HC$sd)
  df.acc.PC <- data.frame(Sites=acc.PC$sites, Richness=acc.PC$richness, SD=acc.PC$sd)
  df.acc.PD <- data.frame(Sites=acc.PD$sites, Richness=acc.PD$richness, SD=acc.PD$sd)

  PD.col <- "#bfbfbf"; PD.col2 <- "#494949"; PD.col3 <- "#848484"
  PC.col <- "#ed7d31"; PC.col2 <- "#f3a977"
  HC.col <- "#5b9bd5"; HC.col2 <- "#99c1e5"
  
  p1 <- ggplot() +
    theme_bw() +
    # PD data
    # geom_point(data=df.acc.PD, aes(x=Sites, y=Richness), alpha=1.5, color = PD.col2) +
    geom_line(data=df.acc.PD, aes(x=Sites, y=Richness), size = 1, alpha=0.6, color = PD.col2) +
    geom_ribbon(data=df.acc.PD, aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2, fill = PD.col3) +
    # PC data
    # geom_point(data=df.acc.PC, aes(x=Sites, y=Richness), alpha=1.5, color = PC.col) +
    geom_line(data=df.acc.PC, aes(x=Sites, y=Richness), size = 1, alpha=0.6, color = PC.col) +
    geom_ribbon(data=df.acc.PC, aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2, fill = PC.col) +
    # HC data
    # geom_point(data=df.acc.HC, aes(x=Sites, y=Richness), alpha=1.5, color = HC.col) +
    geom_line(data=df.acc.HC, aes(x=Sites, y=Richness), size = 1, alpha=0.6, color = HC.col) +
    geom_ribbon(data=df.acc.HC, aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2, fill = HC.col) +
    labs(x = "Number of Samples", y = paste0("Number of ", featuretype)) +
    theme(strip.background = element_blank(),
          panel.grid = element_blank())
  
  ggsave(p1, filename = 
           paste0("data/Community_Composition/Rarefaction_Plots/Rarefaction_Plot_", cohort, "_", featuretype, ".svg"),
         width = 4, height = 3)
  
  return(p1)
  
}

#--------------------------------------------------------------------------------------------------
# Linear Regression Analysis - Alpha Diversity
#--------------------------------------------------------------------------------------------------


AlphaLinearRegressionQuantilePlot <- function(df, x, x2, y, color, fill, ylabel, title){
  
  PD.col <- "#bfbfbf"
  PC.col <- "#ed7d31"
  HC.col <- "#5b9bd5"
  
  p1 <- ggplot(df, aes(x=x, y=y, color=color, fill=fill)) +
    geom_smooth(method="lm", se=F) +
    geom_point(aes(fill=fill),shape=21, size=2, alpha = 0.9) +
    labs(x="Filtered Sample Reads", y=ylabel, fill="group") +
    theme_bw() +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), 
                 color = color), label.x = 15e6, label.y.npc=0.25, hjust = 0) +
    ggtitle(title) +
    scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
    scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.position = c(0.9, 0.5))
  p1
  
  p2 <- ggplot(df, aes(x=x2, y=y, color=color)) +
    geom_boxplot() +
    labs(x="Filtered Sample Read Quantiles", y=ylabel) +
    theme_bw() +
    scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())
  p2
  
  c1 <- cowplot::plot_grid(p1, p2, nrow = 2, align = "v")
  return(c1)
}

#--------------------------------------------------------------------------------------------------
# Linear Regression Analysis - Beta Diversity PCoA 1 & 2
#--------------------------------------------------------------------------------------------------

BetaLinearRegressionPlot <- function(df, x, y, y2, color, fill, feature, title){
  
  PD.col <- "#bfbfbf"
  PC.col <- "#ed7d31"
  HC.col <- "#5b9bd5"
  
  
  p1 <- ggplot(df, aes(x=x, y=y, color=color, fill=fill)) +
    geom_smooth(method="lm", se=F) +
    geom_point(aes(fill=fill),shape=21, size=2, alpha = 0.9) +
    labs(x="Filtered Sample Reads", y=paste0("PCoA Axis 1: ", feature)) +
    theme_bw() +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"),
                 color = color), label.x = 1.75e7, label.y.npc=0.25, hjust = 0) +
    ggtitle(title) +
    scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
    scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  
  p2 <- ggplot(df, aes(x=x, y=y2, color=color, fill=fill)) +
    geom_smooth(method="lm", se=F) +
    geom_point(aes(fill=fill),shape=21, size=2, alpha = 0.9) +
    labs(x="Filtered Sample Reads", y=paste0("PCoA Axis 2: ", feature), fill="group") +
    theme_bw() +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"),
                 color = color), label.x = 1.75e7, label.y.npc=0.25, hjust = 0) +
    scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
    scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE)  +
    theme(legend.position = c(0.85, 0.8),
          panel.grid = element_blank())
  
  c1 <- cowplot::plot_grid(p1, p2, nrow = 2, align = "v")
  return(c1)
}


#--------------------------------------------------------------------------------------------------
# DMM Analysis - Determine Optimal Clusters
#--------------------------------------------------------------------------------------------------


DMM_fit <- function(datObj, nmax){
  
  #' Function fits a DMM model to a given phlyoseq Obj
  #' Input Phyloseq Obj and the max number of clusters to fit
  #' Outputs' Laplace value and Model Fit plot
  
  # Add prevalence threshold
  dat.DMM <- datObj %>% 
    core(detection = 0, prevalence = 0.1)
  # Calculate Pseudocounts 
  count <- PseudoCounts(dat.DMM, reads) %>% 
    t() %>% as.matrix()
  # Fit the DMM model. Set the maximum allowed number of community types to 6 to speed up the analysis
  set.seed(42)
  fit <- lapply(1:nmax, dmn, count = count, verbose=TRUE)
  # Check model fit with different number of mixture components using standard information criteria
  lplc <- sapply(fit, laplace)
  aic  <- sapply(fit, AIC) 
  bic  <- sapply(fit, BIC) 
  
  dmm.fit <- 
    data.frame(cluster = 1:length(lplc), 
               Laplace = lplc, AIC = aic, BIC = bic) %>% 
    pivot_longer(-cluster, names_to = "model.metrics")
  
  dmm.fit.plot <- 
    ggplot(data = dmm.fit, aes(x = cluster, y = value, color = model.metrics)) +
    geom_point() +
    geom_line() +
    scale_color_d3() +
    labs(x="Number of Dirichlet Components", y="Model Fit", 
         color = "Model Metrics") +
    theme_classic() +
    theme(legend.position = c(0.2, 0.7))
  dmm.fit.plot
  
  output <- list(fit, lplc, dmm.fit.plot)
  names(output) <- c("fit", "laplace", "plot")
  return(output)
}

#--------------------------------------------------------------------------------------------------
# DMM Analysis - Samples per cluster DF 
#--------------------------------------------------------------------------------------------------

DMM_stats <- function(best.fit){
  
  # Sample-component assignments
  sas <- apply(mixture(best.fit), 1, which.max) %>% 
    as.data.frame() %>% dplyr::rename(cluster = ".")
  sas <- group_col_from_ids(sas, rownames(sas))
  # Summary Stats by group
  sas.stats <- sas %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::count(group) %>% 
    transmute(n, group, Percentage=n/sum(n)*100)
  
  output <- list(sas, sas.stats)
  names(output) <- c("sas", "sas.stats")
  return(output)
}

#--------------------------------------------------------------------------------------------------
# DMM Analysis - Sample Cluster Distribution
#--------------------------------------------------------------------------------------------------

DMM_cluster_plot <- function(sas.stats){
  
  #'Function plots sample count and percentage of group per cluster
  #'Input df of group stats 
  #'Output plots
  
  cluster_distribution <- 
    ggplot(data = sas.stats, aes(x = as.factor(cluster), y = n, fill = group)) + 
    geom_bar(stat = "identity", width = 0.75) +
    theme_bw() +
    labs(x = "Cluster", y = "Number of Samples") +
    scale_fill_manual(values = cols.pdpchc) +
    theme(panel.grid = element_blank())
  cluster_distribution2 <- 
    ggplot(data = sas.stats, aes(x = as.factor(cluster), y = Percentage, fill = group)) + 
    geom_bar(stat = "identity", width = 0.75) +
    theme_bw() +
    labs(x = "Cluster", y = "Percentage of Cluster") +
    scale_fill_manual(values = cols.pdpchc) +
    theme(panel.grid = element_blank())
  

  legend <- cowplot::plot_grid(get_legend(cluster_distribution))
  cluster_distribution <- cluster_distribution + theme(legend.position = "none")
  cluster_distribution2 <- cluster_distribution2 + theme(legend.position = "none")
  cluster_distribution_final <- 
    cowplot::plot_grid(cluster_distribution, cluster_distribution2, legend,
                       ncol = 3, rel_widths = c(4,4,1))
  
  return(cluster_distribution_final)
  
}

#--------------------------------------------------------------------------------------------------
# DMM Analysis - cluster filter
#--------------------------------------------------------------------------------------------------

DMM_select_cluster <- function(df, cluster_n){
  
  clust_label <- paste0("value.k", cluster_n)
  colnames(df) <- c("feature", "cluster", clust_label)
  df <- subset(df, cluster == cluster_n) %>%
    mutate(feature = factor(feature, levels = unique(feature))) %>% 
    dplyr::select(-cluster)
  return(df)
}

#--------------------------------------------------------------------------------------------------
# DMM Analysis - distinguishing features plotting function
#--------------------------------------------------------------------------------------------------

DMM_cluster_driver_plot <- function(df, yval, comparison){
  
  # # Troubleshoot
  # df = K2vK1.df
  # yval =  K2vK1.df$K2vK1
  # comparison = "2/1"
  
  LFC.KvK <- ggplot(df, aes(x = reorder(feature, yval), y = yval)) +
    geom_bar(stat = "identity",  width = 0.5) +
    labs(title = paste("Distinguishing features: Clusters ", comparison), y = expression(log[2]*" fold change")) +
    geom_rangeframe() + 
    theme_tufte() +
    coord_flip() +
    theme(axis.title.y = element_blank()) 
  
  return(LFC.KvK)
  
}

#--------------------------------------------------------------------------------------------------
#                                        Manual LM models
#--------------------------------------------------------------------------------------------------
# 
# lm.PdPc <- function(metadf, metric){
#   ###' Function conducts Linear Model for PD vs PC
#   env.PdPc <- filter(metadf, donor_group != "HC")
#   formula <- as.formula(
#     paste(metric, "~", paste(c("description", "host_age_factor", "host_body_mass_index", "sex"), collapse="+") ) )
#   
#   linear.model <- lm(formula, data=env.PdPc, na.action = na.omit)
#   plot_model(linear.model, show.values = TRUE, value.offset = .3)
#   qqnorm(resid(linear.model))
#   qqline(resid(linear.model))
#   dev.off()
#   return(linear.model)
# }
# 
# 
# lmm.PdHc <- function(metadf, metric){
#   ###' Function conducts Linear Mixed Model for PD vs HC
#   env.PdHc <- filter(metadf, paired != "No")
#   
#   formula <- as.formula(paste(metric, "~", paste(c("description"))))
#   lmm <- lme(formula, random= ~ 1 | paired, data=env.PdHc, na.action = na.omit)
#   qqnorm(resid(lmm))
#   qqline(resid(lmm))
#   return(lmm)
# }
#--------------------------------------------------------------------------------------------------
