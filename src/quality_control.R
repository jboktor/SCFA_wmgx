rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")



#----------------------------------------------
#  Plot Seq Depth Distributions by group
#----------------------------------------------

cols.pdpchc <- c("PD"= "#bfbfbf", 
                 "PC" = "#ed7d31",
                 "HC" = "#5b9bd5")
cols.pdpchc_dark <- c("PD"= "#494949", 
                      "PC" = "#ed7d31",
                      "HC" = "#5b9bd5")
cols.pdpchc.rim <- c("PD"= "#494949", 
                     "PC" = "#c15811",
                     "HC" = "#2e75b5")


# reads <- load_reads("Merged") %>% 
#   na.omit()

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

knead_reads_TBC <-
  read_tsv(
    "files/TBC_biobakery_output_slim/kneaddata/merged/kneaddata_read_count_table.tsv",
    col_names = T ) %>% 
  dplyr::rename("# samples" = "Sample")
func_reads_TBC <-
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
    col_names = T
  )
reads_TBC <-
  func_reads_TBC %>%
  left_join(knead_reads_TBC, by = "# samples") %>%
  dplyr::filter(`# samples` %ni%  negative_controls) %>%
  dplyr::mutate(`# samples` = substr(`# samples`, 1, 10)) %>%
  left_join(TBC_keys, by = "# samples") %>%
  dplyr::rename("donor_id" = "id", "clean_total_reads" = "total reads") %>% 
  dplyr::mutate(cohort = "TBC")

knead_reads_RUSH <-
  read_tsv(
    "files/RUSH_biobakery_output_slim/kneaddata/merged/kneaddata_read_count_table.tsv",
    col_names = T ) %>% 
  dplyr::rename("# samples" = "Sample")
func_reads_RUSH <-
  read_tsv(
    "files/RUSH_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
    col_names = T
  )
reads_RUSH <-
  func_reads_RUSH %>%
  left_join(knead_reads_RUSH, by = "# samples") %>%
  dplyr::filter(str_detect(`# samples`, "BLANK", negate = TRUE)) %>%
  dplyr::filter(str_detect(`# samples`, "MSA", negate = TRUE)) %>%
  dplyr::mutate(`# samples` = substr(`# samples`, 1, 6)) %>%
  left_join(RUSH_keys, by = "# samples") %>%
  dplyr::rename("clean_total_reads" = "total reads") %>% 
  dplyr::mutate(cohort = "RUSH")

reads <- rbind(reads_TBC, reads_RUSH)
df.reads <- group_col_from_ids(reads, reads$donor_id)
df.reads$group <- factor(df.reads$group, levels=c("PC", "PD", "HC"))

#--------------------------------------------------------------------
#                      Low Quality Reads
#--------------------------------------------------------------------

low_qc <-
  df.reads %>% 
  filter(clean_total_reads < 1000000 | is.na(clean_total_reads)) %>% 
  select(donor_id)
  
save(low_qc, file="files/low_quality_samples.RData")

# Samples to re-sequence
# df.reads %>% 
#   filter(clean_total_reads < 3000000 | is.na(clean_total_reads)) %>% 
#   select(`# samples`)
# select(donor_id)

#--------------------------------------------------------------------
#                          Plotting
#--------------------------------------------------------------------
histo_plot <- 
  ggplot(df.reads, aes(x=clean_total_reads, color = group), alpha = 0.3) + 
  theme_bw() +
  geom_density() +
  facet_wrap(~cohort) +
  scale_color_manual(values = cols.pdpchc)  +
  theme(axis.title.x = element_blank(),
        legend.position = c(0.9, 0.5),
        panel.grid = element_blank())

ecdf_plot <- 
  ggplot(df.reads, aes(x=clean_total_reads, colour = group)) + 
  stat_ecdf(geom = "step", pad = FALSE, alpha = 0.5) +
  stat_ecdf(geom = "point", pad = FALSE, alpha = 0.9, size = 0.75) +
  theme_bw() +
  facet_wrap(~cohort) +
  labs(x ="Clean Reads", y = "ECDF") +
  scale_color_manual(values = cols.pdpchc, guide = FALSE)  +
  theme(legend.position = c(0.9, 0.5),
        panel.grid = element_blank())
c1 <- cowplot::plot_grid(histo_plot, ecdf_plot, ncol = 1, align="v")
c1
ggsave(c1, filename = "data/Community_Composition/Quality_Control/Sequencing_Depth_Density_&_ECDF_Merged_Facet.png",
       width = 6, height = 5, dpi = 600)




################################################################################################
################################# Seq Depth Boxplots by group #################################  

# Test for normality - Distb is non-normal
shapiro.test(df.reads$clean_total_reads)
# KW AND Wilcoxon test
KW <- kruskal.test(clean_total_reads ~ group, data = df.reads)
DunnB <- dunnTest(clean_total_reads ~ group, data = df.reads, method="bh")
PCvsHC.stat <- DunnB$res$P.adj[1]
PCvsPD.stat <- DunnB$res$P.adj[3]
### Paired Wilcoxon Test between PD and HC 
m <- meta(dat.species) %>% dplyr::select(paired) %>% rownames_to_column(var = "donor_id")
paired.donors <- left_join(df.reads, m, var = "donor_id") %>% 
  filter(paired != "No")
stat.test <- wilcox.test(clean_total_reads ~ group, 
                         paired = F, data = paired.donors)
HCvsPD.stat <- stat.test$p.value


p2 <- ggplot(df.reads, aes(x=group, y=clean_total_reads)) +
  geom_point(aes(fill = group, color = group), position = position_jitterdodge(dodge.width=1),shape=21, size=1.25, alpha = 1) +
  geom_boxplot(aes(fill = group), outlier.alpha = 0, alpha = 0.3, width = 0.45) +
  theme_minimal() +
  labs(y="Total Clean Reads per Sample") +
  scale_fill_manual(values = cols.pdpchc) +
  scale_color_manual(values = cols.pdpchc.rim) +
  geom_signif(comparisons = list(c("HC", "PC")), tip_length = 0, 
              annotations = sig_mapper(PCvsHC.stat,  symbols = F, porq = "p")) +
  geom_signif(comparisons = list(c("PC", "PD")), tip_length = 0.01, 
              y_position = 1.95e7, annotations = sig_mapper(PCvsPD.stat,  symbols = F, porq = "p")) +
  geom_signif(comparisons = list(c("HC", "PD")), tip_length = 0.01, 
              y_position = 1.95e7, annotations = sig_mapper(HCvsPD.stat, symbols = F, porq = "p")) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank())

p3 <- ggplot(df.reads, aes(x=group, y=clean_total_reads)) +
  geom_point(aes(fill = group, color = group), position = position_jitterdodge(dodge.width=1),shape=21, size=1.25, alpha = 1) +
  geom_boxplot(aes(fill = group), outlier.alpha = 0, alpha = 0.3, width = 0.45) +
  theme_minimal() +
  labs(y="Total Clean Reads per Sample") +
  scale_fill_manual(values = cols.pdpchc) +
  scale_color_manual(values = cols.pdpchc.rim) +
  facet_wrap(~cohort) +
  # geom_signif(comparisons = list(c("HC", "PC")), tip_length = 0, 
  #             annotations = sig_mapper(PCvsHC.stat,  symbols = F, porq = "p")) +
  # geom_signif(comparisons = list(c("PC", "PD")), tip_length = 0.01, 
  #             y_position = 1.95e7, annotations = sig_mapper(PCvsPD.stat,  symbols = F, porq = "p")) +
  # geom_signif(comparisons = list(c("HC", "PD")), tip_length = 0.01, 
  #             y_position = 1.95e7, annotations = sig_mapper(HCvsPD.stat, symbols = F, porq = "p")) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank())

ggsave(p2, filename = "data/Community_Composition/Quality_Control/Sequencing_Depth_Boxplt.svg",
       width = 3, height = 6)




################################################################################################
#################################  Distrubution Plots by Quantile  ################################# 

# Goal is to see if there are any taxa detected at deeper sequencing at one group and not others
# In other words, is a deeper sequencing responsible for detection of rare taxa that are significant?

############## Species Abundnace by Read Depth : Distributions between donor group
## Pull abundance data - merge with clean_total_reads_factor & donor_group -> melt df

env <- mutate(env, clean_total_reads_factor_tight = 
                as.factor(cut(env$clean_total_reads, 
                              breaks = quantile(env$clean_total_reads, probs=seq(0, 1, 0.2)),
                              labels= c(1:5), include.lowest=TRUE) ))

seqdat <- dplyr::select(env, c("donor_id", "clean_total_reads_factor", "clean_total_reads_factor_tight", "donor_group"))

ab <- dat %>% 
  microbiome::transform("compositional") %>% abundances() %>% 
  t() %>% as.data.frame.matrix() %>% 
  rownames_to_column() %>% cbind(seqdat)
abm <- melt(ab)

# ArcSinSqrt Transform Abundances
abm$value <- asin(sqrt(abm$value))

# Bin Abundance by Quantiles
abm <-  mutate(abm, abundance_factor = 
                 if_else(abm$value > 0.4, 4,
                         if_else(abm$value > 0.1, 3,
                                 if_else(abm$value > 0, 2, 1 
                                 ))))

abm$abundance_factor <- as.factor(abm$abundance_factor)
PD.col <- "#bfbfbf"; PC.col <- "#ed7d31"; HC.col <- "#5b9bd5"

z1 <- abm %>% filter(abundance_factor != 1 ) %>% 
  ggplot(aes(x=value, color = donor_group, fill = donor_group), alpha = 0.4) +
  geom_histogram(position="dodge", boundary = 0) +
  facet_grid(cols=vars(clean_total_reads_factor_tight), 
             rows=vars(abundance_factor), scales = "free") +
  scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE)  +
  theme_classic() +
  labs(x = "Arcsin(Sqrt(Abundance))", y = "Histogram by Binned Feature Abundance") +
  ggtitle("Sequencing Depth Quantiles") +
  theme(plot.title = element_text(hjust = 0.5))
# z1
ggsave(z1, filename = "data/Community_Composition/Quality_Control/SeqDepth_&_Abundance_Quantile_Histogram_Matrix.png",
       width = 20, height = 8, dpi = 1200)

z2 <- abm %>% filter(abundance_factor != 1 ) %>% 
  ggplot(aes(x=value, color = donor_group)) +
  stat_ecdf(geom = "point", pad = FALSE, alpha =0.1) +
  stat_ecdf(geom = "step", pad = FALSE) +
  facet_grid(cols=vars(clean_total_reads_factor_tight), 
             rows=vars(abundance_factor), scales = "free") +
  labs(x = "Arcsin(Sqrt(Abundance))", y = "Binned Feature Abundance - ECDF") +
  ggtitle("Sequencing Depth Quantiles") +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col))  +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
# z2
ggsave(z2, filename = "data/Community_Composition/Quality_Control/SeqDepth_&_Abundance_Quantile_ECDF_matrix.png",
       width = 20, height = 8, dpi = 1200)




z1b <- abm %>% filter(abundance_factor == 2 & clean_total_reads_factor_tight == 1 ) %>% 
  ggplot(aes(x=value, color = donor_group, fill = donor_group), alpha = 0.4) +
  geom_histogram(position="dodge", boundary = 0) +
  facet_grid(cols=vars(clean_total_reads_factor_tight), 
             rows=vars(abundance_factor), scales = "free") +
  scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE)  +
  theme_classic() +
  ggtitle("Sequencing Depth: First Quantile vs AST transformed Abundance (0 < x < 0.1)") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
# z1b

z2b <- abm %>% filter(abundance_factor == 2 & clean_total_reads_factor_tight == 1 ) %>%
  ggplot(aes(x=value, color = donor_group)) +
  stat_ecdf(geom = "step", pad = FALSE) +
  facet_grid(cols=vars(clean_total_reads_factor_tight), 
             rows=vars(abundance_factor), scales = "free") +
  labs(x = "Arcsin(Sqrt(Abundance))", y = "Binned Feature Abundance - ECDF") +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col))  +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
# z2b

QC_speciesAbundanceVsreads <- cowplot::plot_grid(z1b, z2b, align = "h", nrow = 2)
# QC_speciesAbundanceVsreads

ggsave(QC_speciesAbundanceVsreads, filename = "data/Community_Composition/Quality_Control/SeqDepth_Quantile_Distribution.png",
       width = 8, height = 6, dpi = 1200)




################################################################################################
#####################  Seq-Depth x Alpha-Diversity : Linear Regression  ######################## 

# Prep Reads 
reads <- dplyr::select(func_reads, c("# samples","total reads")) %>% 
  dplyr::rename( "id" = "# samples", "clean_total_reads" = "total reads")

# Prep Metadata 
env <- meta(dat)
env <- left_join(env, reads, by = "id")
env <- mutate(env, clean_total_reads_factor = as.factor(cut(env$clean_total_reads, 
                                                            breaks = quantile(env$clean_total_reads), 
                                                            labels=c(1,2,3,4), include.lowest=T) ))

env$description <- factor(env$description, levels=c("PD Patient", "Population Control", "Household Control"))
env$donor_group <- factor(env$donor_group, levels=c("PC", "PD", "HC"))


# Objects to loop over
obj <- c(dat, dat.path.slim, dat.ec.slim, dat.KOs.slim, dat.Pfams.slim, dat.Eggnogs.slim)
obj.label <- c("Species", "Pathways", "Enzymes", "KOs", "Pfams", "Eggnogs")

# Initiate count
cnt <- 1

for (i in obj) {
  
  cat("\n\n"); cat("Processing", obj.label[cnt], "for Alpha Diversity x Seq-Depth Linear Regression Analysis \n\n")
  
  # Prep abundnace table
  dat_alpha <-  transform(i, "compositional")
  ## Calculate Alpha Diversity Metrics and add cols to df
  env$observed <- alpha(abundances(dat_alpha), 'observed')$observed
  env$shannon <- alpha(abundances(dat_alpha), 'shannon')$diversity_shannon
  env$evenness <- evenness(abundances(dat_alpha), 'simpson')$simpson
  
  
  observed.plot <- AlphaLinearRegressionQuantilePlot(df=env, x=env$clean_total_reads, x2=env$clean_total_reads_factor, y=env$observed,
                                                     color=env$donor_group, fill=env$donor_group, 
                                                     ylabel=paste("Observed", obj.label[cnt]), title= paste("Observed", obj.label[cnt], "by Read Depth"))
  
  shannon.plot <- AlphaLinearRegressionQuantilePlot(df=env, x=env$clean_total_reads, x2=env$clean_total_reads_factor, y=env$shannon,
                                                    color=env$donor_group, fill=env$donor_group,
                                                    ylabel=paste("Shannon Diversity:", obj.label[cnt]), title="Shannon Index by Read Depth")
  
  evenness.plot <- AlphaLinearRegressionQuantilePlot(df=env, x=env$clean_total_reads, x2=env$clean_total_reads_factor, y=env$evenness,  
                                                     color=env$donor_group, fill=env$donor_group,
                                                     ylabel= paste("Simpsons Evenness:", obj.label[cnt]), title="Simpsons Evenness by Read Depth")
  
  QC_alphaVsreads <- cowplot::plot_grid(observed.plot, shannon.plot, evenness.plot, 
                                        ncol = 3, align = "h")
  
  ggsave(QC_alphaVsreads, filename = paste0("data/Community_Composition/Quality_Control/Read_Depth_vs_Alpha_Diversity/Read_Depth_vs_",  obj.label[cnt], "_AlphaDiversity.png"),
         width = 16, height = 8)
  
  cnt <- cnt + 1
}


################################################################################################
#####################  Seq-Depth x Beta-Diversity : Linear Regression  ######################## 


# Initiate count
cnt <- 1

for (i in obj) {
  
  ## AITCHISONS DISTANCE
  cat("\n\n"); cat("Processing", obj.label[cnt], "for Beta Diversity x Seq-Depth Linear Regression Analysis \n\n")
  
  obj_clr <- microbiome::transform(i, "clr")
  iDist <- distance(obj_clr, method="euclidean")
  iMDS  <- ordinate(obj_clr, "MDS", distance=iDist)
  #  PCoA for Axis 1 and 2
  
  p <- plot_ordination(obj_clr, iMDS, color="description", axes = c(1, 2))
  pcoa.data = p$data
  df.pcoa <- left_join(pcoa.data, reads, by = "id")
  df.pcoa$donor_group <- factor(df.pcoa$donor_group, levels=c("PC", "PD", "HC"))
  
  QC_betaVsreads <- BetaLinearRegressionPlot(df=df.pcoa, x=df.pcoa$clean_total_reads, y=df.pcoa$Axis.1, y2=df.pcoa$Axis.2,
                                             color=df.pcoa$donor_group, fill=df.pcoa$donor_group,
                                             feature=obj.label[cnt], title=paste0(obj.label[cnt], ": Beta Diversity by Read Depth"))
  
  ggsave(QC_betaVsreads, filename = paste0("data/Community_Composition/Quality_Control/Read_Depth_vs_Beta_Diversity/Read_Depth_vs_",  obj.label[cnt], "_BetaDiversity.png"),
         width = 8, height = 8)
  
  cnt <- cnt + 1
}



