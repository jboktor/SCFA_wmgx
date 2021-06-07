# Gut-Brain Module (GBM) and Gut-Metabolic Module (GMM) Analysis

rm(list = ls())
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")
source("src/daf_functions.R")
wkd <- getwd()

#------------------------------------------------------------------------------
#                              Setup
#------------------------------------------------------------------------------
# Download Binary from: https://github.com/omixer/omixer-rpmR
# Install using R CMD INSTALL omixeRpm_x.y.z.tar.gz (with appropriate version)
# VERSION USED: 0.3.2
# A Java Development Kit (JDK) is also required
library(omixerRpm)

#------------------------------------------------------------------------------
#             Generate Gut-Brain / Gut-Metabolic Modules
#------------------------------------------------------------------------------
datfiles <- list.files("files/Phyloseq_objects_Woltka/")
datfilepaths <- paste0("files/Phyloseq_objects_Woltka/", datfiles)
lapply(datfilepaths, base::load, .GlobalEnv)

# Filter only KO number from new annotated IDs
input_kos <- 
  dat.KOs %>% 
  # subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  microbiome::abundances() %>% 
  data.frame() %>% 
  rownames_to_column(var = "entry") #%>% 
  # dplyr::mutate(entry = substr(entry, 6, 11))

# Run the module mapping on the loaded table.
gmm <- rpm(input_kos, minimum.coverage=0.3, annotation = 1, 
           module.db	= loadDB("GMMs.v1.07"))
gmm.abundance <- gmm@abundance
gmm.annotation <- gmm@annotation
gmm.coverage <- gmm@coverage

gbm <- rpm(input_kos, minimum.coverage=0.3, annotation = 1, 
            module.db	= loadDB("GBMs.v1.0"))
gbm.abundance <- gbm@abundance
gbm.annotation <- gbm@annotation
gbm.coverage <- gbm@coverage

# Load the default mapping database
db_gmm <- loadDB("GMMs.v1.07")
db_gbm <- loadDB("GBMs.v1.0")


# Create data frames with abundances and modules along with translated name
####### GBM ####### 
translated_annotations_gbm <- c()
modname_list_gbm <- c()
for (mod in gbm@annotation){
  for (modname in mod){
    modname_list_gbm <- c(modname_list_gbm, modname)
    print(modname)
    print(as.character(omixerRpm::getNames(db_gbm, modname) ))
    translated_annotations_gbm <- c(translated_annotations_gbm, omixerRpm::getNames(db_gbm, modname) )
  }
}
df.gbm <- data.frame(row.names = translated_annotations_gbm,
                     "module"=modname_list_gbm,
                     gbm.abundance)

####### GMM ####### 
translated_annotations_gmm <- c()
modname_list_gmm <- c()
for (mod in gmm@annotation){
  for (modname in mod){
    modname_list_gmm <- c(modname_list_gmm, modname)
    print(modname)
    print(as.character(omixerRpm::getNames(db_gmm, modname) ))
    translated_annotations_gmm <- c(translated_annotations_gmm, omixerRpm::getNames(db_gmm, modname) )
  }
}
df.gmm <- data.frame(row.names=translated_annotations_gmm,
                     "module"=modname_list_gmm,
                     gmm.abundance)

#---------------------------------------------------------------------------------
#                            Create GM Phyloseq Objects
#---------------------------------------------------------------------------------

### Metadata 
my_sample_data <- meta(dat.species) %>% sample_data()

#------------------ Gut Brain Modules ------------------ 

my_GBM.ab_table <- df.gbm[-1] %>% 
  t() %>% 
  otu_table(taxa_are_rows=F)
dat.GBMs <- phyloseq(my_GBM.ab_table, my_sample_data)
print(dat.GBMs)
save(dat.GBMs, file = "files/Phyloseq_objects_Woltka/GBMs_PhyloseqObj.RData")

#------------------ Gut Metabolic Modules ------------------ 

my_GMM.ab_table <- df.gmm[-1] %>% 
  t() %>% 
  otu_table(taxa_are_rows=F)
dat.GMMs <- phyloseq(my_GMM.ab_table, my_sample_data)
print(dat.GMMs)
save(dat.GMMs, file = "files/Phyloseq_objects_Woltka/GMMs_PhyloseqObj.RData")

#_________________________
#  GBMs  MaAsLin2   ----
#_________________________

df_metadata <- process_meta(dat.species, verbose = T) %>% 
  mutate(diet = factor(diet,  levels = c("Control", "Prebiotic"))) %>% 
  mutate(genotype = factor(genotype,  levels = c("WT", "ASO")))

dat.object <- maaslin_prep(dat.GBMs) %>% 
  microbiome::transform("compositional")
obj.name <- "GBMs"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = T)


#_________________________
#  GMMs  MaAsLin2   ----
#_________________________

df_metadata <- process_meta(dat.species, verbose = T) %>% 
  mutate(diet = factor(diet,  levels = c("Control", "Prebiotic"))) %>% 
  mutate(genotype = factor(genotype,  levels = c("WT", "ASO")))

dat.object <- maaslin_prep(dat.GMMs) %>% 
  microbiome::transform("compositional")
obj.name <- "GMMs"
variance_plot(dat.object)

df_input <- dat.object %>% 
  variance_filter(filter.percent = 0) %>% 
  sqrt() %>% asin()

fit_models(obj.name,
           df_input = df_input,
           df_metadata = df_metadata,
           plot_scatter = T)


#-------------------------------------------------------------------------------
####                             Plotting                                #### 
#-------------------------------------------------------------------------------

load("files/Phyloseq_Merged/GBMs_PhyloseqObj.RData")
load("files/Phyloseq_Merged/GMMs_PhyloseqObj.RData")

plot_dafs(obj.name = "GBMs", obj = dat.GBMs, cohort = "Merged")
plot_daf_summary(obj.name = "GBMs", obj = dat.GBMs, cohort = "Merged", 
                 repelyup = 0.002, repelydn = 0.002)

plot_dafs(obj.name = "GMMs", obj = dat.GMMs, cohort = "Merged")
plot_daf_summary(obj.name = "GMMs", obj = dat.GMMs, cohort = "Merged", 
                 repelyup = 0.001, repelydn = 0.001)



