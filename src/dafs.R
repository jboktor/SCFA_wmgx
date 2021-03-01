# Differentially Abundant Features (DAFs)

######## Load Data & functions
rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")
source("src/daf_functions.R")
load("files/low_quality_samples.RData")
wkd <- getwd()

#--------------------------------------------------------------------------------
#####                              All Cohorts                                #### 
#--------------------------------------------------------------------------------
remove_dats()
load_all_cohorts()

# Metadata Selection
corr.meta.clinical <-
  dat.species %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  process_meta(cohort = "Merged") %>%
  select(donor_id, contains(c(motor_severity_scores_summary,
      clinical_variables))) %>% 
  column_to_rownames(var = "donor_id")

#-------------------------
#####  Species #### 
#-------------------------

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(dat.species)
# PD v PC abundance data
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# Plot Variance Estimate
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.1)
# PD v HC PAIRED abundance data
dat_pdhc = subset_samples(dat.object, paired !="No")
# Plot Variance Estimate
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.1)
fit_models(dat = dat.object, 
           obj.name = "Species", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 1,
           plot_scatter = T, 
           cohort = "TBC")

#-------------------------
#####  Genera #### 
#-------------------------

dat.object <- maaslin_prep(dat.genus)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
fit_models(dat = dat.object, 
           obj.name = "Genera", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 1,
           plot_scatter = T)

# # Motor Severity Score Calculation
# df_input_data <- dat.object %>% 
#   microbiome::transform("compositional") %>% 
#   variance_filter(filter.percent = 0)
# obj.name = "Genera"
# fit_data = Maaslin2(
#   input_data = df_input_data, 
#   input_metadata = corr.meta.clinical, 
#   output = paste0(wkd, "/data/MaAsLin2_Analysis/Motor_Severity_Associations/", 
#                   obj.name, "_maaslin2_output"), 
#   min_prevalence = 0,
#   random_effects = c("paired", "cohort"),
#   fixed_effects = c(motor_severity_scores_summary, clinical_variables),
#   analysis_method = "LM",
#   normalization = "NONE",
#   transform = "AST",
#   cores = 6,
#   plot_scatter = T
# )

#-------------------------
#####  Phylum  #### 
#-------------------------

dat.object <- maaslin_prep(dat.phylum)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
fit_models(dat = dat.object, 
           obj.name = "Phylum", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 1,
           plot_scatter = T)

#-------------------------
##### Pathways slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.path.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
fit_models(dat = dat.object, 
           obj.name = "Pathways.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F)

#-------------------------
##### Enzymes slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.ec.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
fit_models(dat = dat.object, 
           obj.name = "Enzymes.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F)

#-------------------------
#### Kegg Orthology slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.KOs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(0.2)
fit_models(dat = dat.object, 
           obj.name = "KOs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F)

#-------------------------
#### Gene Ontology ####
#-------------------------

dat.object <- maaslin_prep(dat.GOs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
fit_models(dat = dat.object, 
           obj.name = "GOs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F)


# GO Molecular Function
abund_rename <-
  dat.object %>% 
  abundances() %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  filter(grepl(".MF.", rowname) ) %>% 
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows=T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.MF <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.MF, cohort = "Merged", tag = "_MF")

# GO Biological Processes
abund_rename <-
  dat.object %>% 
  abundances() %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  filter(grepl(".BP.", rowname) ) %>% 
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows=T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.BP <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.BP, cohort = "Merged", tag = "_BP")

# GO Cellular Compartments
abund_rename <-
  dat.object %>% 
  abundances() %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  filter(grepl(".CC.", rowname) ) %>% 
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows=T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.CC <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.CC, cohort = "Merged", tag = "_CC")


#-------------------------
#   EGGNOGS slim
#-------------------------

dat.object <- maaslin_prep(dat.EGGNOGs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.7)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.7)
fit_models(dat = dat.object, 
           obj.name = "EGGNOGs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F)

#-------------------------
#####  PFAMS slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.PFAMs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
fit_models(dat = dat.object, 
           obj.name = "PFAMs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F)


#--------------------------------------------------------------------------------
#####                              TBC ONLY                                #### 
#--------------------------------------------------------------------------------

remove_dats()
load_tbc()

#-------------------------
#####  Species #### 
#-------------------------

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(dat.species)
# PD v PC abundance data
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# Plot Variance Estimate
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.1)
# PD v HC PAIRED abundance data
dat_pdhc = subset_samples(dat.object, paired !="No")
# Plot Variance Estimate
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.1)
fit_models(dat = dat.object, 
           obj.name = "Species", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 1,
           plot_scatter = T,
           cohort = "TBC")

#-------------------------
#####  Genera #### 
#-------------------------

dat.object <- maaslin_prep(dat.genus)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
fit_models(dat = dat.object, 
           obj.name = "Genera", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 1,
           plot_scatter = T,
           cohort = "TBC")

#-------------------------
#####  Phylum  #### 
#-------------------------

dat.object <- maaslin_prep(dat.phylum)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
fit_models(dat = dat.object, 
           obj.name = "Phylum", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 1,
           plot_scatter = T,
           cohort = "TBC")

#-------------------------
##### Pathways slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.path.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
fit_models(dat = dat.object, 
           obj.name = "Pathways.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "TBC")

#-------------------------
##### Enzymes slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.ec.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
fit_models(dat = dat.object, 
           obj.name = "Enzymes.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "TBC")

#-------------------------
#### Kegg Orthology slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.KOs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(0.2)
fit_models(dat = dat.object, 
           obj.name = "KOs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "TBC")

#-------------------------
#### Gene Ontology ####
#-------------------------

dat.object <- maaslin_prep(dat.GOs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
fit_models(dat = dat.object, 
           obj.name = "GOs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "TBC")


# GO Molecular Function
abund_rename <-
  dat.object %>% 
  abundances() %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  filter(grepl(".MF.", rowname) ) %>% 
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows=T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.MF <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.MF, cohort = "TBC", tag = "_MF")

# GO Biological Processes
abund_rename <-
  dat.object %>% 
  abundances() %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  filter(grepl(".BP.", rowname) ) %>% 
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows=T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.BP <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.BP, cohort = "TBC", tag = "_BP")

# GO Cellular Compartments
abund_rename <-
  dat.object %>% 
  abundances() %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  filter(grepl(".CC.", rowname) ) %>% 
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows=T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.CC <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.CC, cohort = "TBC", tag = "_CC")


#-------------------------
#   EGGNOGS slim
#-------------------------

dat.object <- maaslin_prep(dat.EGGNOGs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.7)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.7)
fit_models(dat = dat.object, 
           obj.name = "EGGNOGs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "TBC")

#-------------------------
#####  PFAMS slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.PFAMs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
fit_models(dat = dat.object, 
           obj.name = "PFAMs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "TBC")


#--------------------------------------------------------------------------------
#####                              RUSH ONLY                                #### 
#--------------------------------------------------------------------------------

remove_dats()
load_rush()

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(dat.species)
# PD v PC abundance data
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# Plot Variance Estimate
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.1)
# PD v HC PAIRED abundance data
dat_pdhc = subset_samples(dat.object, paired !="No")
# Plot Variance Estimate
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.1)
fit_models(dat = dat.object, 
           obj.name = "Species", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 1,
           plot_scatter = T,
           cohort = "RUSH")

#-------------------------
#####  Genera #### 
#-------------------------

dat.object <- maaslin_prep(dat.genus)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
fit_models(dat = dat.object, 
           obj.name = "Genera", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 1,
           plot_scatter = T,
           cohort = "RUSH")

#-------------------------
#####  Phylum  #### 
#-------------------------

dat.object <- maaslin_prep(dat.phylum)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0)
fit_models(dat = dat.object, 
           obj.name = "Phylum", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 1,
           plot_scatter = T,
           cohort = "RUSH")

#-------------------------
##### Pathways slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.path.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
fit_models(dat = dat.object, 
           obj.name = "Pathways.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "RUSH")

#-------------------------
##### Enzymes slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.ec.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
fit_models(dat = dat.object, 
           obj.name = "Enzymes.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "RUSH")

#-------------------------
#### Kegg Orthology slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.KOs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(0.2)
fit_models(dat = dat.object, 
           obj.name = "KOs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "RUSH")

#-------------------------
#### Gene Ontology ####
#-------------------------

dat.object <- maaslin_prep(dat.GOs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.3)
fit_models(dat = dat.object, 
           obj.name = "GOs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "RUSH")


# GO Molecular Function
abund_rename <-
  dat.object %>% 
  abundances() %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  filter(grepl(".MF.", rowname) ) %>% 
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows=T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.MF <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.MF, cohort = "RUSH", tag = "_MF")

# GO Biological Processes
abund_rename <-
  dat.object %>% 
  abundances() %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  filter(grepl(".BP.", rowname) ) %>% 
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows=T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.BP <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.BP, cohort = "RUSH", tag = "_BP")

# GO Cellular Compartments
abund_rename <-
  dat.object %>% 
  abundances() %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  filter(grepl(".CC.", rowname) ) %>% 
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows=T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.CC <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.CC, cohort = "RUSH", tag = "_CC")


#-------------------------
#   EGGNOGS slim
#-------------------------

dat.object <- maaslin_prep(dat.EGGNOGs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.7)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.7)
fit_models(dat = dat.object, 
           obj.name = "EGGNOGs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "RUSH")

#-------------------------
#####  PFAMS slim #### 
#-------------------------

dat.object <- maaslin_prep(dat.PFAMs.slim)
dat_pdpc = subset_samples(dat.object, donor_group !="HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
dat_pdhc = subset_samples(dat.object, paired !="No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% 
  microbiome::transform("compositional") %>% 
  variance_filter(filter.percent = 0.2)
fit_models(dat = dat.object, 
           obj.name = "PFAMs.slim", 
           dat_pdpc = dat_pdpc, 
           dat_pdhc = dat_pdhc,
           df_input_data_pdhc = df_input_data_pdhc, 
           df_input_data_pdpc = df_input_data_pdpc,
           cores = 6,
           plot_scatter = F,
           cohort = "RUSH")






#--------------------------------------------
########   Plot Summaries MERGED    ########   
#--------------------------------------------

remove_dats()
load_all_cohorts()

plot_dafs(obj.name = "Phylum", obj = dat.phylum, cohort = "Merged")
plot_dafs(obj.name = "Genera", obj = dat.genus, cohort = "Merged")
plot_dafs(obj.name = "Species", obj = dat.species, cohort = "Merged")
plot_daf_summary(obj.name = "Species", obj = dat.species, cohort = "Merged")

plot_dafs(obj.name = "Pathways.slim", obj = dat.path.slim, cohort = "Merged")
plot_daf_summary(obj.name = "Pathways.slim", obj = dat.path.slim, cohort = "Merged", 
                 repelyup = 0, repelydn = 0)

plot_dafs(obj.name = "Enzymes.slim", obj = dat.ec.slim, cohort = "Merged")
plot_daf_summary(obj.name = "Enzymes.slim", obj = dat.ec.slim, cohort = "Merged", 
                 repelyup = 0, repelydn = 0)

plot_dafs(obj.name = "KOs.slim", obj = dat.KOs.slim, cohort = "Merged")
plot_daf_summary(obj.name = "KOs.slim", obj = dat.KOs.slim, cohort = "Merged", 
                 repelyup = 0.0025, repelydn = 0)

plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim, cohort = "Merged", tag = "")
plot_daf_summary(obj.name = "GOs.slim", obj = dat.GOs.slim, cohort = "Merged", 
                 repelyup = 0.001, repelydn = 0.001)

plot_dafs(obj.name = "EGGNOGs.slim", obj = dat.EGGNOGs.slim, cohort = "Merged")
plot_daf_summary(obj.name = "EGGNOGs.slim", obj = dat.EGGNOGs.slim, cohort = "Merged", 
                 repelyup = 0.002, repelydn = 0.002)

plot_dafs(obj.name = "PFAMs.slim", obj = dat.PFAMs.slim, cohort = "Merged")
plot_daf_summary(obj.name = "PFAMs.slim", obj = dat.PFAMs.slim, cohort = "Merged", 
                 repelyup = 0.001, repelydn = 0.001)


#--------------------------------------------
##########   Plot Summaries TBC  ##########   
#--------------------------------------------

remove_dats()
load_rush()

plot_dafs(obj.name = "Phylum", obj = dat.phylum, cohort = "TBC")
plot_dafs(obj.name = "Genera", obj = dat.genus, cohort = "TBC")
plot_dafs(obj.name = "Species", obj = dat.species, cohort = "TBC")
plot_daf_summary(obj.name = "Species", obj = dat.species, cohort = "TBC")

plot_dafs(obj.name = "Pathways.slim", obj = dat.path.slim, cohort = "TBC")
plot_daf_summary(obj.name = "Pathways.slim", obj = dat.path.slim, cohort = "TBC", 
                 repelyup = 0, repelydn = 0)

plot_dafs(obj.name = "Enzymes.slim", obj = dat.ec.slim, cohort = "TBC")
plot_daf_summary(obj.name = "Enzymes.slim", obj = dat.ec.slim, cohort = "TBC", 
                 repelyup = 0, repelydn = 0)

plot_dafs(obj.name = "KOs.slim", obj = dat.KOs.slim, cohort = "TBC")
plot_daf_summary(obj.name = "KOs.slim", obj = dat.KOs.slim, cohort = "TBC", 
                 repelyup = 0.0025, repelydn = 0.005)

plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim, cohort = "TBC", tag = "")
plot_daf_summary(obj.name = "GOs.slim", obj = dat.GOs.slim, cohort = "TBC", 
                 repelyup = 0.001, repelydn = 0.001)

plot_dafs(obj.name = "EGGNOGs.slim", obj = dat.EGGNOGs.slim, cohort = "TBC")
plot_daf_summary(obj.name = "EGGNOGs.slim", obj = dat.EGGNOGs.slim, cohort = "TBC", 
                 repelyup = 0.002, repelydn = 0.002)

plot_dafs(obj.name = "PFAMs.slim", obj = dat.PFAMs.slim, cohort = "TBC")
plot_daf_summary(obj.name = "PFAMs.slim", obj = dat.PFAMs.slim, cohort = "TBC", 
                 repelyup = 0.001, repelydn = 0.001)


#--------------------------------------------
##########     Plot Summaries RUSH    ##########   
#--------------------------------------------

remove_dats()
load_rush()

plot_dafs(obj.name = "Phylum", obj = dat.phylum, cohort = "RUSH")
plot_dafs(obj.name = "Genera", obj = dat.genus, cohort = "RUSH")
plot_dafs(obj.name = "Species", obj = dat.species, cohort = "RUSH")
plot_daf_summary(obj.name = "Species", obj = dat.species, cohort = "RUSH")

plot_dafs(obj.name = "Pathways.slim", obj = dat.path.slim, cohort = "RUSH")
plot_daf_summary(obj.name = "Pathways.slim", obj = dat.path.slim, cohort = "RUSH", 
                 repelyup = 0, repelydn = 0)

plot_dafs(obj.name = "Enzymes.slim", obj = dat.ec.slim, cohort = "RUSH")
plot_daf_summary(obj.name = "Enzymes.slim", obj = dat.ec.slim, cohort = "RUSH", 
                 repelyup = 0, repelydn = 0)

plot_dafs(obj.name = "KOs.slim", obj = dat.KOs.slim, cohort = "RUSH")
plot_daf_summary(obj.name = "KOs.slim", obj = dat.KOs.slim, cohort = "RUSH", 
                 repelyup = 0.0025, repelydn = 0.005)

plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim, cohort = "RUSH", tag = "")
plot_daf_summary(obj.name = "GOs.slim", obj = dat.GOs.slim, cohort = "RUSH", 
                 repelyup = 0.001, repelydn = 0.001)

plot_dafs(obj.name = "EGGNOGs.slim", obj = dat.EGGNOGs.slim, cohort = "RUSH")
plot_daf_summary(obj.name = "EGGNOGs.slim", obj = dat.EGGNOGs.slim, cohort = "RUSH", 
                 repelyup = 0.002, repelydn = 0.002)

plot_dafs(obj.name = "PFAMs.slim", obj = dat.PFAMs.slim, cohort = "RUSH")
plot_daf_summary(obj.name = "PFAMs.slim", obj = dat.PFAMs.slim, cohort = "RUSH", 
                 repelyup = 0.001, repelydn = 0.001)






