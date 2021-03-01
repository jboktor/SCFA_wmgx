### Figure 2) Differentially Abundant Taxa

rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")
source("src/daf_functions.R")
load("files/low_quality_samples.RData")
wkd <- getwd()

######### INPUT SPECIES 
load_all_cohorts()
LEV <- dat.species %>% 
  subset_samples(donor_id %ni% low_qc[[1]])
lev <- "Species"

############# Visualization Transformations ############# 
taxa_names(LEV) <- gsub("s__", "", taxa_names(LEV))
dat_obj <- microbiome::transform(LEV, "compositional")
## ArcSinSqrt() Transformation
otu_table(dat_obj) <- asin(sqrt(otu_table(dat_obj)))
## Color Schemes
cols.pdpc <- c("PD"= "#bfbfbf", "PC" = "#ed7d31")
cols.pdhc <- c("PD"= "#bfbfbf", "HC" = "#5b9bd5")
# Rims
cols.pdpc.rim <- c("PD"= "#494949", "PC" = "#c15811")
cols.pdhc.rim <- c("PD"= "#494949", "HC" = "#2e75b5")

############# data prep ############# 
# PD v PC
dat_pdpc = subset_samples(dat_obj, donor_group !="HC")
abun.pdpc <- as.data.frame.matrix(abundances(dat_pdpc))
# PD v HC PAIRED
dat_pdhc = subset_samples(dat_obj, paired !="No")
abun.pdhc <- as.data.frame.matrix(abundances(dat_pdhc))


############# Read-in Maaslin Files - all features used in significance testing ############# 
 
Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/Merged/", lev, "_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Population Control")
Maas.pd.pc$feature <- gsub("s__", "", Maas.pd.pc$feature)
Maas.pd.pc.sig <- Maas.pd.pc %>% filter(qval < 0.25)

Maas.pd.hc <- read_tsv(paste0("data/MaAsLin2_Analysis/Merged/", lev, "_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Household Control")
Maas.pd.hc$feature <- gsub("s__", "", Maas.pd.hc$feature)
Maas.pd.hc.sig <- Maas.pd.hc %>% filter(qval < 0.25) 


#'#########  select significant features from abundance tables  ##############
#'#########  Pull Genus/Phylum annotations for each significant feature #########  
####### PD v PC
abun.pdpc <- rownames_to_column(abun.pdpc)
abun.pdpc.filtered <- filter(abun.pdpc, rowname %in% Maas.pd.pc.sig$feature)
abun.pdpc.inpt <- abun.pdpc.filtered %>% column_to_rownames(var="rowname") %>% 
  t() %>% melt() %>% mutate(group = if_else(grepl(".PC", Var1), "PC", "PD"))
## Geom Tile data - Genus and Phylum levels
abun.pdpc.filtered$speciesname <- paste0("s__", abun.pdpc.filtered$rowname)
abun.pdpc.inpt.phylo <- taxa_genus_phlyum_annotation(dat.species, abun.pdpc.filtered$speciesname) 


####### PD v HC PAIRED
abun.pdhc <- rownames_to_column(abun.pdhc)
abun.pdhc.filtered <- filter(abun.pdhc, rowname %in% Maas.pd.hc.sig$feature) 
abun.pdhc.inpt <- abun.pdhc.filtered %>% column_to_rownames(var="rowname")  %>% 
  t() %>% melt() %>% mutate(group = if_else(grepl(".HC", Var1), "HC", "PD"))
## Geom Tile data - Genus and Phylum levels
abun.pdhc.filtered$speciesname <- paste0("s__", abun.pdhc.filtered$rowname)
abun.pdhc.inpt.phylo <- taxa_genus_phlyum_annotation(dat.species, abun.pdhc.filtered$speciesname) 


## Create BarTile Color Palette - Manual process 
all.phylo <- rbind(abun.pdpc.inpt.phylo, abun.pdhc.inpt.phylo)
unique(all.phylo$Phylum)
phylum.cols <-
  c(
    "Actinobacteria" = "#d62728",
    "Firmicutes" = "#a6cee3",
    "Bacteroidetes" = "#2CA02C",
    "Proteobacteria" = "#6a3d9a",
    "Verrucomicrobia" = "#6a3d1a"
  )
tile.cols <- phylum.cols

######################################################################## 
##########################    PD v PC Plots   ##########################
######################################################################## 

# Subset of phlyoseq obj subset to get samples of interest
dat_pdpc.PD = subset_samples(dat_pdpc, donor_group =="PD")
dat_pdpc.PDabun <- as.data.frame.matrix(abundances(dat_pdpc.PD)) %>% 
  rownames_to_column() %>%  filter(rowname %in% Maas.pd.pc.sig$feature) %>% column_to_rownames()

dat_pdpc.PC = subset_samples(dat_pdpc, donor_group =="PC")
dat_pdpc.PCabun <- as.data.frame.matrix(abundances(dat_pdpc.PC)) %>% 
  rownames_to_column() %>%  filter(rowname %in% Maas.pd.pc.sig$feature) %>% column_to_rownames()

###### PHYLUM TILES ###### 
qsig <- Maas.pd.pc.sig %>% dplyr::select(c(feature, qval)) %>%  dplyr::rename("Species" = "feature")
phylo.pc <- boxplot_phylobars(inpt.phylo = abun.pdpc.inpt.phylo, sigvals = qsig, tile.cols = tile.cols)

###### Generalized or pseudo-fold calculation ###### 
gfc_data <- generalized_fold_change(dat_pdpc.PDabun, dat_pdpc.PCabun)
PDovrPC <- tibble("feature" = rownames(dat_pdpc.PCabun), "gFC" = gfc_data)
###### Generalized Fold Change (gFC) BarPlot ###### 
PDovrPC.BP <- PDovrPC
PDovrPC.BP <- mutate(PDovrPC.BP, direction = if_else(PDovrPC.BP$gFC > 0, "PD",
                                                     if_else(PDovrPC.BP$gFC < 0, "PC",  "error")))
PDovrPC.BP$feature <- factor(PDovrPC.BP$feature, levels = rev(phylo.pc$Axis.order)) 
g0 <- gfc_plot(PDovrPC.BP, cols.pdpc, alfa = 0.8)

###### Significance Plot ###### 
sigplot.df.pdpc <- dplyr::select(Maas.pd.pc.sig, c("feature", "pval", "qval")) %>%  melt()
sigplot.df.pdpc$feature <- factor(sigplot.df.pdpc$feature, levels = rev(phylo.pc$Axis.order)) 
g2 <- significance_barplot(sigplot.df.pdpc)

###### BoxPlots ###### 
## Prepping Significance labels
abun.pdpc.inpt <- daf_boxplot_sigvalues(sigplot.df.pdpc, abun.pdpc.inpt)
abun.pdpc.inpt$Var2 <- factor(abun.pdpc.inpt$Var2, levels = rev(phylo.pc$Axis.order)) 
g1 <-
  daf_boxplots(
    abun.pdpc.inpt,
    fill_cols = cols.pdpc,
    rim_cols = cols.pdpc.rim,
    alfa = 0.2,
    obj.name = lev
  )

###### Prevalence Plot ######
# Subset of phlyoseq obj subset to get samples of interest
dat_pdpc.PDprev <- tibble::enframe(prevalence(dat_pdpc.PD)) %>% filter(name %in% Maas.pd.pc.sig$feature) 
colnames(dat_pdpc.PDprev) <- c("feature", "PD")
dat_pdpc.PCprev <- tibble::enframe(prevalence(dat_pdpc.PC)) %>% filter(name %in% Maas.pd.pc.sig$feature)
colnames(dat_pdpc.PCprev) <- c("feature", "PC")
dat_pdpc.PREV <- left_join(dat_pdpc.PDprev, dat_pdpc.PCprev, by = "feature") %>% melt()
dat_pdpc.PREV$feature <- factor(dat_pdpc.PREV$feature,  rev(phylo.pc$Axis.order))
dat_pdpc.PREV$variable <- factor(dat_pdpc.PREV$variable, levels = c("PC", "PD"))
g3 <- prevalence_barplot(dat_pdpc.PREV, cols.pdpc, alfa = 0.8)



######################################################################## 
##########################    PD v HC Plots   ########################## 
######################################################################## 

# Subset of phlyoseq obj subset to get samples of interest
dat_pdhc.PD = subset_samples(dat_pdhc, donor_group =="PD")
dat_pdhc.PDabun <- as.data.frame.matrix(abundances(dat_pdhc.PD))  %>% 
  rownames_to_column() %>%  filter(rowname %in% Maas.pd.hc.sig$feature) %>% column_to_rownames()

dat_pdhc.HC = subset_samples(dat_pdhc, donor_group =="HC")
dat_pdhc.HCabun <- as.data.frame.matrix(abundances(dat_pdhc.HC)) %>% 
  rownames_to_column() %>%  filter(rowname %in% Maas.pd.hc.sig$feature) %>% column_to_rownames()


###### PHYLUM TILES ###### 
qsig <- Maas.pd.hc.sig %>% dplyr::select(c(feature, qval)) %>%  dplyr::rename("Species" = "feature")
phylo.hc <- boxplot_phylobars(inpt.phylo = abun.pdhc.inpt.phylo, sigvals = qsig, tile.cols = tile.cols)


######  Generalized or pseudo-fold change 
gfc_data <- generalized_fold_change(dat_pdhc.PDabun, dat_pdhc.HCabun)
PDovrHC <- tibble("feature" = rownames(dat_pdhc.PDabun), "gFC" = gfc_data)
###### Generalized Fold Change (gFC) BarPlot ###### 
PDovrHC.BP <- PDovrHC
PDovrHC.BP <- mutate(PDovrHC.BP, direction = if_else(PDovrHC.BP$gFC > 0, "PD",
                     if_else(PDovrHC.BP$gFC < 0, "HC",  "error")))
PDovrHC.BP$feature <- factor(PDovrHC.BP$feature, levels = rev(phylo.hc$Axis.order)) 
h0 <- gfc_plot(PDovrHC.BP, cols.pdhc, alfa = 1)

###### Significance Plot ###### 
sigplot.df.pdhc <- dplyr::select(Maas.pd.hc.sig, c("feature", "pval", "qval")) %>%  melt()
sigplot.df.pdhc$feature <- factor(sigplot.df.pdhc$feature, levels = rev(phylo.hc$Axis.order)) 
h2 <- significance_barplot(sigplot.df.pdhc)

###### BoxPlots ###### 
## Prepping Significance labels
abun.pdhc.inpt <- daf_boxplot_sigvalues(sigplot.df.pdhc, abun.pdhc.inpt)
abun.pdhc.inpt$Var2 <- factor(abun.pdhc.inpt$Var2, levels = rev(phylo.hc$Axis.order)) 
h1 <-
  daf_boxplots(
    abun.pdhc.inpt,
    fill_cols = cols.pdhc,
    rim_cols = cols.pdhc.rim,
    alfa = 0.2,
    obj.name = lev
  )

###### Prevalence Plot ######
# Subset of phlyoseq obj subset to get samples of interest
dat_pdhc.PDprev <- tibble::enframe(prevalence(dat_pdhc.PD)) %>% filter(name %in% Maas.pd.hc.sig$feature) 
colnames(dat_pdhc.PDprev) <- c("feature", "PD")
dat_pdhc.HCprev <- tibble::enframe(prevalence(dat_pdhc.HC)) %>% filter(name %in% Maas.pd.hc.sig$feature)
colnames(dat_pdhc.HCprev) <- c("feature", "HC")
dat_pdhc.PREV <- left_join(dat_pdhc.PDprev, dat_pdhc.HCprev, by = "feature") %>% melt()
dat_pdhc.PREV$feature <- factor(dat_pdhc.PREV$feature, levels = rev(phylo.hc$Axis.order))
dat_pdhc.PREV$variable <- factor(dat_pdhc.PREV$variable, levels = c("HC", "PD"))
h3 <- prevalence_barplot(dat_pdhc.PREV, cols.pdhc, alfa = 0.8)


######################## Merge Panels ######################## 

g.bars <- phylo.pc$Bars + theme(axis.text.x = element_blank()) + ggtitle(" ")
g.legend <- phylo.pc$Legends
g0a <- g0 + theme(axis.title.x = element_blank(), 
                  axis.text.y = element_blank())
g1a <- g1 + theme(axis.title.x = element_blank(), 
                  legend.position = "none", axis.text.y = element_blank())
g3a <- g3 + theme(axis.title.x = element_blank(), 
                  axis.text.y = element_blank(), legend.position = "none")

h.bars <- phylo.hc$Bars + ggtitle(" ")
h.legend <- phylo.hc$Legends
h0a <- h0 + theme(axis.text.y = element_blank())
h1a <- h1 + theme(legend.position = "none", axis.text.y = element_blank())
h3a <- h3 + theme(axis.text.y = element_blank(), legend.position = "none")


# Setting plot length variables 
top_len <- length(unique(PDovrPC.BP$feature)) + 1
bottom_len <- length(unique(PDovrHC.BP$feature)) + 2

DAF_part1 <- cowplot::plot_grid(g.bars, h.bars, ncol = 1, align = "hv", 
                       labels = "AUTO",
                       rel_heights = c(top_len, bottom_len))

DAF_part2 <- cowplot::plot_grid(g1a, g3a, g0a, 
                                h1a, h3a, h0a, 
                                nrow = 2, ncol=3, align = "h", 
                                rel_heights = c(top_len, bottom_len),
                                rel_widths = c(3, 1, 1.5))


DAF_final <- cowplot::plot_grid(DAF_part1, DAF_part2, ncol = 2, align = "hv", 
                                rel_widths = c(1, 5.25))
DAF_final

# ggsave(DAF_final, filename = paste0("figures/Figure_2/DAF_Figure_2_new.svg"),
#        width = 14, height = 7)
# 
# # Legends
# ggsave(phylo.pc$Legends, filename = paste0("figures/Figure_2/DAF_Figure_2_PC.legend_new.svg"),
#        width = 7, height = 7)
# ggsave(phylo.hc$Legends, filename = paste0("figures/Figure_2/DAF_Figure_2_HC.legend_new.svg"),
#        width = 7, height = 7)
# 

