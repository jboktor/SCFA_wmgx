### QC_Seq_Depth Analysis

rm(list = ls())
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")
source("src/alpha_diversity.R")
source("src/beta_diversity.R")
base::load("files/Phyloseq_objects_Woltka.RData")
base::load("files/Phyloseq_objects_Woltka_Rarefied.RData")


#------------------------------------------
#             Alpha Diversity
#------------------------------------------

alpha_diversity_summary(x = phylo_dats_rare, z = names(phylo_dats_rare), 
                       tag = "Rarefied")
alpha_diversity_summary(x = phylo_dats, z = names(phylo_dats), 
                       tag = "Non-Rarefied")


#------------------------------------------
#             Beta Diversity
#------------------------------------------

beta_diversity_summary(x = phylo_dats_rare, z = names(phylo_dats_rare), 
                       tag = "Rarefied",  dist = "Aitchisons")
beta_diversity_summary(x = phylo_dats, z = names(phylo_dats), 
                       tag = "Non-Rarefied",  dist = "Aitchisons")
beta_diversity_summary(x = phylo_dats_rare, z = names(phylo_dats_rare), 
                       tag = "Rarefied",  dist = "bray")
beta_diversity_summary(x = phylo_dats, z = names(phylo_dats), 
                       tag = "Non-Rarefied",  dist = "bray")


#______________________________________________________________________________
#                      Taxonomy Rel Abundance Bars   -----  
#______________________________________________________________________________

dat.obj <- 
  dat.genus

dat.top.30 <- dat.obj %>% 
  get_top_taxa(n=30, relative = TRUE, discard_other = F, other_label = "Other")
dat.top.20 <- dat.obj %>% 
  get_top_taxa(n=20, relative = TRUE, discard_other = F, other_label = "Other")
dat.top.15 <- dat.obj %>% 
  get_top_taxa(n=15, relative = TRUE, discard_other = F, other_label = "Other")


dat.obj %>% 
  abundance_heatmap(treatment = "treatment_group")

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
    facet_by = "genotype",
    grid_by = "diet",
    facet_cols = 2,
    order_alg = "hclust",
    base_color = "#5b9bd5",
    palette = barcols,
    # color_levels = barcol_ID
    ) +
  labs(y = "Relative Abundance") +
  theme(axis.text.x = element_blank())

ggsave(barplt1, filename = "data/Community_Composition/Stacked_Barplots/Top30_Genera_Cohort_Facet.png",
       width = 8, height = 5.75)

# Merge cohort specific donor groups
dat.top2 <- merge_samples(dat.top.30, "treatment_group")

# Summary Barplot
barplt2 <- 
  fantaxtic_bar(
    dat.top2,
    color_by = "Order",
    label_by = "Genus",
    other_label = "Other",
    order_alg = "hclust",
    palette = barcols
  ) +
  labs(y = "Relative Abundance", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_blank(), 
        plot.margin = unit(c(1, 1, 1, 1), "cm")
  )
å
ggsave(barplt2, filename = "data/Community_Composition/Stacked_Barplots/Top30_Genera_Group_Facet_Summary.png",
       width = 4, height = 6)
  


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





