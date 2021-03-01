### QC_Seq_Depth Analysis

rm(list = ls())
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")
source("src/alpha_diversity.R")
source("src/beta_diversity.R")
load("files/low_quality_samples.RData")

#------------------------------------------
#             Alpha/Beta Diversity
#------------------------------------------

remove_dats()
load_tbc()
# alpha_diversity_summary("TBC")
withCallingHandlers(beta_diversity_summary("TBC"), warning=function(w){invokeRestart("muffleWarning")})
withCallingHandlers(beta_diversity_summary("TBC", dist = "bray"), warning=function(w){invokeRestart("muffleWarning")})
withCallingHandlers(beta_diversity_summary("TBC", dist = "jaccard"), warning=function(w){invokeRestart("muffleWarning")})

remove_dats()
load_rush()
# alpha_diversity_summary("RUSH")
withCallingHandlers(beta_diversity_summary("RUSH"), warning=function(w){invokeRestart("muffleWarning")})
withCallingHandlers(beta_diversity_summary("RUSH", dist = "bray"), warning=function(w){invokeRestart("muffleWarning")})
withCallingHandlers(beta_diversity_summary("RUSH", dist = "jaccard"), warning=function(w){invokeRestart("muffleWarning")})

remove_dats()
load_all_cohorts()
# alpha_diversity_summary("Merged")
withCallingHandlers(beta_diversity_summary("Merged"), warning=function(w){invokeRestart("muffleWarning")})
withCallingHandlers(beta_diversity_summary("Merged", dist = "bray"), warning=function(w){invokeRestart("muffleWarning")})
withCallingHandlers(beta_diversity_summary("Merged", dist = "jaccard"), warning=function(w){invokeRestart("muffleWarning")})

#----------------------------------------------
#  Plot Rarefaction Curves per group
#----------------------------------------------

##' Manual Input here: to create the species, genus, phylum, pathways, enzyme,
##' KO, Eggnog, or Pfam Rarefaction plot input the appropriate __ dat __ frame
##' below This analysis does 1000 random permutations per number of samples
##' (1-(n of group)) and plots richness along with it's 95% confidence interval

remove_dats()
load_tbc()
reads <- load_reads("TBC")
acc.1.TBC <- feature_accumulation_plot(dat.species, featuretype = "Species", reads, cohort = "TBC")
acc.2.TBC <- feature_accumulation_plot(dat.genus, featuretype = "Genera", reads, cohort = "TBC")
acc.3.TBC <- feature_accumulation_plot(dat.path.slim, featuretype = "Pathways", reads, cohort = "TBC")
acc.4.TBC <- feature_accumulation_plot(dat.KOs.slim, featuretype = "KOs", reads, cohort = "TBC")
acc.5.TBC <- feature_accumulation_plot(dat.PFAMs.slim, featuretype = "Pfams", reads, cohort = "TBC")
acc.6.TBC <- feature_accumulation_plot(dat.EGGNOGs.slim, featuretype = "Eggnogs", reads, cohort = "TBC")

remove_dats()
load_rush()
reads <- load_reads("RUSH")
acc.1.RUSH <- feature_accumulation_plot(dat.species, featuretype = "Species", reads, cohort = "RUSH")
acc.2.RUSH <- feature_accumulation_plot(dat.genus, featuretype = "Genera", reads, cohort = "RUSH")
acc.3.RUSH <- feature_accumulation_plot(dat.path.slim, featuretype = "Pathways", reads, cohort = "RUSH")
acc.4.RUSH <- feature_accumulation_plot(dat.KOs.slim, featuretype = "KOs", reads, cohort = "RUSH")
acc.5.RUSH <- feature_accumulation_plot(dat.PFAMs.slim, featuretype = "Pfams", reads, cohort = "RUSH")
acc.6.RUSH <- feature_accumulation_plot(dat.Eggnogs.slim, featuretype = "Eggnogs", reads, cohort = "RUSH")

remove_dats()
load_all_cohorts()
reads <- load_reads("Merged")
acc.1.Merged <- feature_accumulation_plot(dat.species, featuretype = "Species", reads, cohort = "Merged")
acc.2.Merged <- feature_accumulation_plot(dat.genus, featuretype = "Genera", reads, cohort = "Merged")
acc.3.Merged <- feature_accumulation_plot(dat.path.slim, featuretype = "Pathways", reads, cohort = "Merged")
acc.4.Merged <- feature_accumulation_plot(dat.KOs.slim, featuretype = "KOs", reads, cohort = "Merged")
acc.5.Merged <- feature_accumulation_plot(dat.PFAMs.slim, featuretype = "Pfams", reads, cohort = "Merged")
acc.6.Merged <- feature_accumulation_plot(dat.EGGNOGs.slim, featuretype = "Eggnogs", reads, cohort = "Merged")



################################################################################################
#################################  Plot Abundance Bars by group #################################  

remove_dats()
load_all_cohorts()

dat.obj <- 
  dat.genus %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) %>%
  core(detection = 0, prevalence = 0.1)

# Create Metadata Column for Cohort x Donor Group
sample_data(dat.obj)$cohort_donor_group <- 
  paste(sample_data(dat.obj)$cohort, sample_data(dat.obj)$donor_group)
# Abundance filter for top 30 Genera
dat.top.30 <- dat.obj %>% 
  get_top_taxa(n=30, relative = TRUE, discard_other = F, other_label = "Other")
dat.top.20 <- dat.obj %>% 
  get_top_taxa(n=20, relative = TRUE, discard_other = F, other_label = "Other")
dat.top.15 <- dat.obj %>% 
  get_top_taxa(n=15, relative = TRUE, discard_other = F, other_label = "Other")

barcols <- c(
  "#386cb0",
  "#7fc97f",
  "#beaed4",
  "#fdc086",
  "#ffff99",
  "#f0027f",
  "#bf5b17",
  "#666666",
  "#7fc97f",
  "#beaed4"
)

# Plot all Samples
barplt1 <- 
  fantaxtic_bar(
    dat.top.30,
    color_by = "Order",
    label_by = "Genus",
    other_label = "Other",
    facet_by = "donor_group",
    grid_by = "cohort",
    facet_cols = 3,
    order_alg = "hclust",
    # base_color = "#5b9bd5", 
    palette = barcols
    # color_levels = barcol_ID
    ) +
  labs(y = "Relative Abundance") +
  theme(axis.text.x = element_blank())

ggsave(barplt1, filename = "data/Community_Composition/Stacked_Barplots/Top30_Genera_Cohort_Facet.png",
       width = 12, height = 6)

# Merge cohort specific donor groups
dat.top2 <- merge_samples(dat.top.30, "cohort_donor_group")

# Summary Barplot
barplt2 <- 
  fantaxtic_bar(
    dat.top2,
    color_by = "Order",
    label_by = "Genus",
    other_label = "Other",
    facet_by = "cohort",
    facet_cols = 2,
    order_alg = "hclust",
    # base_color = "#5b9bd5"
    palette = barcols
  ) +
  labs(y = "Relative Abundance", x = NULL)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_blank())


ggsave(barplt2, filename = "data/Community_Composition/Stacked_Barplots/Top30_Genera_Cohort_Facet_Summary.png",
       width = 4.5, height = 6)
  


# Abundance filter for top 15 Genera
dat.top3 <- dat.obj %>% 
  get_top_taxa(n=15, relative = TRUE, discard_other = F, other_label = "Other")

# Merge cohort specific donor groups
dat.top3 <- merge_samples(dat.top3, "cohort_donor_group")

# Summary Barplot
barplt2 <- 
  fantaxtic_bar(
    dat.top2,
    color_by = "Order",
    label_by = "Genus",
    other_label = "Other",
    facet_by = "cohort",
    facet_cols = 2,
    order_alg = "hclust",
    base_color = "#5b9bd5"
  ) +
  labs(y = "Relative Abundance")
ggsave(barplt2, filename = "data/Community_Composition/Stacked_Barplots/Top30_Genera_Cohort_Facet_Summary.svg",
       width = 4.5, height = 6)





