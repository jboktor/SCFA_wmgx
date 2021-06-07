# Differentially Abundant Features (DAFs)

# Load data & functions ----
rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")
source("src/daf_functions.R")
datfiles <- list.files("files/Phyloseq_objects_Woltka/")
datfilepaths <- paste0("files/Phyloseq_objects_Woltka/", datfiles)
lapply(datfilepaths, base::load, .GlobalEnv)
wkd <- getwd()

df_metadata <- process_meta(dat.species, verbose = T) %>% 
  mutate(diet = factor(diet,  levels = c("Control", "Prebiotic"))) %>% 
  mutate(genotype = factor(genotype,  levels = c("WT", "ASO")))

#_________________________
#          Species  ---- 
#_________________________
## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(dat.species) %>% 
  microbiome::transform("compositional")
obj.name <- "Species"

# Plot Variance Estimate
variance_plot(dat.object)
df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.1) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = T)

#_________________________
#####  Genera #### 
#_________________________

dat.object <- maaslin_prep(dat.genus) %>% 
  microbiome::transform("compositional")
obj.name <- "Genus"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.1) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = T)

#_________________________
#         Phylum   ----
#_________________________

dat.object <- maaslin_prep(dat.phylum) %>% 
  microbiome::transform("compositional")
obj.name <- "Phylum"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.1) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = T)

#_________________________
#    Super pathways   ----
#_________________________

dat.object <- maaslin_prep(dat.Superpathways) %>% 
  microbiome::transform("compositional")
obj.name <- "Superpathways"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.1) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = T)

#_________________________
#       Pathways   ----
#_________________________

dat.object <- maaslin_prep(dat.Pathways) %>% 
  microbiome::transform("compositional")
obj.name <- "Pathways"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.1) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = T)

#_________________________
#       Reactions   ----
#_________________________

dat.object <- maaslin_prep(dat.Rxns) %>% 
  microbiome::transform("compositional")
obj.name <- "Rxns"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.2) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = F)

#_________________________
#        Enzrxns   ----
#_________________________

dat.object <- maaslin_prep(dat.Enzrxn) %>% 
  microbiome::transform("compositional")
obj.name <- "Enzrxn"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.2) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = F)

#_________________________
#       MC Proteins   ----
#_________________________

dat.object <- maaslin_prep(dat.MCProteins) %>% 
  microbiome::transform("compositional")
obj.name <- "MCProteins"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.2) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = F)


#_________________________
#          GO   ----
#_________________________

dat.object <- maaslin_prep(dat.GOs) %>% 
  microbiome::transform("compositional")
obj.name <- "GOs"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.3) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = F)

#_________________________
#          KOs  ----
#_________________________

dat.object <- maaslin_prep(dat.KOs) %>% 
  microbiome::transform("compositional")
obj.name <- "KOs"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.3) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = F)

#_________________________
#       eggNOGs ----
#_________________________

dat.object <- maaslin_prep(dat.eggNOGs) %>% 
  microbiome::transform("compositional")
obj.name <- "eggNOGs"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0.4) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = F)

#_____________________________
#     Plot Summaries 
#_____________________________


dafs_genotype_by_diet("Species")
dafs_genotype_by_diet("Genera")
dafs_genotype_by_diet("Phylum")
dafs_genotype_by_diet("Superpathways")
dafs_genotype_by_diet("Pathways")
dafs_genotype_by_diet("Rxns")
dafs_genotype_by_diet("Enzrxn")
dafs_genotype_by_diet("KOs")
dafs_genotype_by_diet("GOs")
dafs_genotype_by_diet("eggNOGs")




