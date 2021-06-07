## PERMANOVA script

rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
load("files/Phyloseq_objects_Woltka/Species_PhyloseqObj.RData")


#-------------------------------------------------------------------------------
#####                       PERMANOVA Analysis                             ##### 
#-------------------------------------------------------------------------------

distlist <- c("euclidean", "bray")
sigmetadata <- vector("list", length(distlist))
names(sigmetadata) = distlist

for (j in distlist) {
  sigmetadata[[j]]$siglst <- c()
  sigmetadata[[j]]$almostSiglst <- c()
}
siglst <- c()
env.pa_perm <- tibble()

objs <-
  c(dat.species)
obj_name <-
  c("Species")

#-------------------------------------------------------------------------------
#####                       Aitchisons Distance                            ##### 
#-------------------------------------------------------------------------------

cnt <- 1
for (obj in objs) {
  obj_processed <- obj %>% 
    # subset_samples(donor_id %ni% low_qc[[1]]) %>% 
    microbiome::transform("compositional") %>% 
    microbiome::transform("clr")
  sp <- t(abundances(obj_processed))
  # Run Metadata pre-processing function
  env <- process_meta(obj_processed)
  # remove ID variable and duplicate grouping variable
  env.pa <- env 
  
  for (i in 1:length(env.pa)) { 
    a <- env.pa[,i]
    a.narm <- na.omit(a)
    if (any(is.na(a))) {
      sp.narm <- sp[-attr(a.narm, "na.action"), ]
    } else {
      sp.narm <- sp
    }
    cat("Aitchisons", obj_name[cnt], ": ", colnames(env.pa[i]), "\n")
    meta_ano = adonis(vegdist(sp.narm, method = "euclidean") ~ a.narm, permutations = 9999)
    row2add <-
      cbind(colnames(env.pa[i]),
            meta_ano$aov.tab[1, ],
            length(a.narm),
            obj_name[cnt], 
            "Aitchisons")
    env.pa_perm <- rbind(env.pa_perm, row2add)
    sig <- meta_ano$aov.tab$`Pr(>F)`[1]
  }
  cnt = cnt + 1
}
permanova_aitch <- env.pa_perm %>%
  remove_rownames() %>%
  dplyr::rename(
    "p_value" = "Pr(>F)",
    "vars" = "colnames(env.pa[i])",
    "n_meta" = "length(a.narm)",
    "data_type" = "obj_name[cnt]",
    "distance" = '"Aitchisons"'
  ) 



#-------------------------------------------------------------------------------
#####                       Bray-Curtis Distance                            ##### 
#-------------------------------------------------------------------------------

env.pa_perm <- tibble()
cnt <- 1
for (obj in objs) {
  obj_processed <- obj %>% 
    microbiome::transform("compositional")
  sp <- t(abundances(obj_processed))
  # Run Metadata pre-processing function
  env <- process_meta(obj_processed)
  # remove ID variable and duplicate grouping variable
  env.pa <- env
  ## Adding Alpha Diversity to PERMANOVA
  alpha_diversity <- obj_processed %>% abundances() %>% 
    microbiome::alpha('shannon')
  env.pa$diversity_shannon <- alpha_diversity$diversity_shannon
  
  for (i in 1:length(env.pa)) { 
    a <- env.pa[,i]
    a.narm <- na.omit(a)
    if (any(is.na(a))) {
      sp.narm <- sp[-attr(a.narm, "na.action"), ]
    } else {
      sp.narm <- sp
    }
    cat("Bray-Curtis", obj_name[cnt], ": ", colnames(env.pa[i]), "\n")
    meta_ano = adonis(vegdist(sp.narm, method = "bray") ~ a.narm, permutations = 9999)
    row2add <-
      cbind(colnames(env.pa[i]),
            meta_ano$aov.tab[1, ],
            length(a.narm),
            obj_name[cnt],
            "BrayCurtis")
    env.pa_perm <- rbind(env.pa_perm, row2add)
    sig <- meta_ano$aov.tab$`Pr(>F)`[1]
  }
  cnt = cnt + 1
}
permanova_bray <- env.pa_perm %>%
  remove_rownames() %>%
  dplyr::rename(
    "p_value" = "Pr(>F)",
    "vars" = "colnames(env.pa[i])",
    "n_meta" = "length(a.narm)",
    "data_type" = "obj_name[cnt]",
    "distance" = '"BrayCurtis"'
  )  


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

permdf <- 
  bind_rows(permanova_aitch, permanova_bray) %>% 
  dplyr::mutate(vars = as.character(vars)) %>%
  dplyr::mutate(R2 = R2 * 100) %>%
  group_by(data_type, distance) %>%
  mutate(FDR = p.adjust(p_value, method = 'BH')) %>%
  ungroup()


## Add Metadata category column
permdf <-
  mutate(
    permdf,
    metacat = if_else(
      permdf$vars %in% motortest_vars,
      "Motor Symptom Assessment",
      if_else(
        permdf$vars %in% grouping_vars,
        "Grouping",
        if_else(
          permdf$vars %in% mouse_vars,
          "Sample Characteristics",
          "AA_notmatched"
        )
      )
    )
  )
                          

write.csv(permdf, file = 
            paste0('data/Community_Composition/PERMANOVA/permanova_analysis_', 
                   Sys.Date(), '.csv'))



