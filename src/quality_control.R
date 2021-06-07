rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")

#----------------------------------------------
#  Plot Seq Depth Distributions by group
#----------------------------------------------

datfiles <- list.files("files/Phyloseq_objects_Woltka/")
datfilepaths <- paste0("files/Phyloseq_objects_Woltka/", datfiles)
lapply(datfilepaths, base::load, .GlobalEnv)

all.dats <- c(dat.species, dat.genus, dat.phylum,
              dat.MCProteins, dat.Rxns, dat.Enzrxn,
              dat.Pathways, dat.Pathways, dat.Superpathways,
              dat.GOs, dat.KOs, dat.eggNOGs)

dat.names <- c("Species", "Genera", "Phylum",
              "MetaCyc Proteins", "Reactions", "Enzymes",
              "Pathways", "Super Pathways",
              "Gene Ontology", "Kegg Orthology", "eggNOGs")
names(all.dats) <- dat.names



df_abund <- tibble()
for (dat in dat.names){
  
  df_stats <- 
    all.dats[[dat]] %>% 
    abundances() %>% 
    colSums() %>% 
    as.data.frame() %>% 
    dplyr::rename("sums" = ".") %>% 
    rownames_to_column(var = "donor_id") %>% 
    mutate(data_type = dat)
  df_abund <- bind_rows(df_abund, df_stats)
}
metadat <- 
  dat.species %>% 
  meta()
df_abund <- df_abund %>% 
  left_join(metadat, by = "donor_id")


df_abund %>% 
  ggplot(aes(x=sums, colour = treatment_group)) + 
  stat_ecdf(geom = "step", pad = FALSE, alpha = 0.5) +
  stat_ecdf(geom = "point", pad = FALSE, alpha = 0.9, size = 0.75) +
  theme_bw() +
  facet_wrap(~data_type, scales = "free_x") +
  labs(x ="Clean Reads", y = "ECDF") +
  # scale_color_manual(values = cols.pdpchc, guide = FALSE)  +
  theme(panel.grid = element_blank())




