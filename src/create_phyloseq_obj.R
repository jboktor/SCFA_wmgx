# Create phyloseq objects

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/metaphlanToPhyloseq_Waldron.R")

negative_controls <- c(
  "S00A4-ATCC_MSA_1003_S96",
  "S00A4-neg2_S119",
  "S00A4-neg3_S125",
  "S00A4-neg_S118",
  "S00A4NegExt_P00A4_S94",
  "S00A4NegH2O_P00A4_S95",
  "S00A4_stagPos_S117",
  "BLANK"
)

#---------------------------------------------------------- 
#-                 Metadata Prep - TBC   
#---------------------------------------------------------- 
metadata_TBC <- 
  read.csv(file = "files/metadata_phyloseq_TBC.csv", header= TRUE) %>% 
  mutate(TBC_Subject_ID = format(TBC_Subject_ID, scientific=F))

TBC_keys <- read.csv(file = "files/metadata_keys.csv", header= TRUE) %>% 
  dplyr::select(c(MBI_Sample_ID, id)) %>% 
  mutate(id = gsub("_", ".", id)) %>% 
  mutate(MBI_Sample_ID = as.character(MBI_Sample_ID)) %>% 
  mutate(id = as.character(id))

TBC_keymap <- TBC_keys$MBI_Sample_ID
names(TBC_keymap) <- TBC_keys$id

#load MDS-Survey data and merge with core metadata
select_val <- function(x, na.rm = FALSE) substr(x, start = 0, stop = 1)

MDSUPDRS <- 
  # read.csv(file = "files/MDSUPDRS_20210217.csv", header= T) %>%
  read.csv(file = "files/MDSUPDRS_20200831.csv", header= T) %>%
  dplyr::select(-who_is_filling_out_survey) %>% 
  mutate_at(vars(-TBC_Subject_ID), select_val) %>% 
  mutate_at(vars(-TBC_Subject_ID), as.numeric) %>% 
  mutate(TBC_Subject_ID = format(TBC_Subject_ID, scientific=F)) %>%
  # mutate(TBC_Subject_ID = gsub(" ", "", TBC_Subject_ID)) %>% 
  rowwise() %>% 
  mutate(mds_updrs_survey_total = sum(c_across(where(is.numeric)),na.rm = TRUE))

str(MDSUPDRS)
str(metadata_TBC)
reads.TBC <- load_reads("TBC")

# Merged MDS-UDPRS Part I & II Survey scores with TBC data
metadata_TBC <-
  tst <- 
  left_join(metadata_TBC, MDSUPDRS, by = "TBC_Subject_ID") %>% 
  left_join(reads.TBC, by = "donor_id")
rownames(metadata_TBC) <- metadata_TBC$donor_id
metadata_TBC <- as.data.frame(metadata_TBC)
metadata_TBC[is.na(metadata_TBC)] <- "not provided"

#---------------------------------------------------------- 
#-                TBC -  Taxonomy 
#---------------------------------------------------------- 
met.table <-read_tsv(
  file = "files/TBC_biobakery_output_slim/metaphlan/merged/metaphlan_taxonomic_profiles.tsv",
  col_names = T)

# Select only species rows from 
bugs.species <- met.table %>% 
  dplyr::rename("taxonomy" = `# taxonomy`) %>% 
  filter(grepl("s__", taxonomy)) %>% 
  filter(!grepl("t__", taxonomy)) %>% 
  column_to_rownames(var = "taxonomy") %>% 
  clean.cols.tax() %>%
  dplyr::select(-contains(negative_controls)) %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))

dat.species <- metaphlanToPhyloseq_Waldron(
  tax = bugs.species, metadat = metadata_TBC)
dat.species

#-------- Species Level Object --------
save(dat.species, file = "files/Phyloseq_TBC/Species_PhyloseqObj.RData")
#-------- Genus Level Object --------
dat.genus = tax_glom(dat.species, taxrank = "Genus", NArm = F)
taxa_names(dat.genus) <- tax_table(dat.genus)[,6]
save(dat.genus, file = "files/Phyloseq_TBC/Genus_PhyloseqObj.RData")
#-------- Family Level Object --------
dat.family = tax_glom(dat.species, taxrank = "Family", NArm = F)
taxa_names(dat.family) <- tax_table(dat.family)[,5]
save(dat.family, file = "files/Phyloseq_TBC/Family_PhyloseqObj.RData")
#-------- Order Level Object --------
dat.order = tax_glom(dat.species, taxrank = "Order", NArm = F)
taxa_names(dat.order) <- tax_table(dat.order)[,4]
save(dat.order, file = "files/Phyloseq_TBC/Order_PhyloseqObj.RData")
#-------- Class Level Object --------
dat.class = tax_glom(dat.species, taxrank = "Class", NArm = F)
taxa_names(dat.class) <- tax_table(dat.class)[,3]
save(dat.class, file = "files/Phyloseq_TBC/Class_PhyloseqObj.RData")
#-------- Phylum Level Object --------
dat.phylum = tax_glom(dat.species, taxrank = "Phylum", NArm = F)
taxa_names(dat.phylum) <- tax_table(dat.phylum)[,2]
save(dat.phylum, file = "files/Phyloseq_TBC/Phylum_PhyloseqObj.RData")
#-------- Kingdom Level Object --------
dat.kingdom = tax_glom(dat.species, taxrank = "Kingdom", NArm = F)
taxa_names(dat.kingdom) <- tax_table(dat.kingdom)[,1]
save(dat.kingdom, file = "files/Phyloseq_TBC/Kingdom_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   TBC -  Pathways 
#---------------------------------------------------------- 

path.abund <- 
  read_tsv(file = "files/TBC_biobakery_output_slim/humann/merged/pathabundance_relab.tsv", col_names = T) 
path.abund <- 
  path.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund() %>% 
  filter(!grepl("UNMAPPED", `# Pathway`)) %>% 
  filter(!grepl("UNINTEGRATED", `# Pathway`))
path.abund.slim <- path.abund %>% 
  filter(!grepl("g__", `# Pathway`)) %>% 
  filter(!grepl("unclassified", `# Pathway`)) %>%
  make_rfriendly_rows(passed_column = "# Pathway") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
path.abund <- path.abund %>% 
  make_rfriendly_rows(passed_column = "# Pathway") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
# All Pathway Data
my_pathab_table <- otu_table(path.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.path <- phyloseq(my_pathab_table, my_sample_data)
dat.path
save(dat.path, file = "files/Phyloseq_TBC/Pathways_PhyloseqObj.RData")
# Slim Pathway Data
my_pathab_table <- otu_table(path.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.path.slim <- phyloseq(my_pathab_table, my_sample_data)
dat.path.slim
save(dat.path.slim, file = "files/Phyloseq_TBC/Pathways.slim_PhyloseqObj.RData")



#---------------------------------------------------------- 
#                   TBC -  Enzymes 
#---------------------------------------------------------- 

ec.abund <- 
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/ecs_relab.tsv", 
    col_names = T)
ec.abund <- 
  ec.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_RPK()
ec.abund.slim <- 
  ec.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
ec.abund <- 
  ec.abund %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
# All Enzyme Data
my_EC.ab_table <- otu_table(ec.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.ec <- phyloseq(my_EC.ab_table, my_sample_data)
dat.ec
save(dat.ec, file = "files/Phyloseq_TBC/Enzymes_PhyloseqObj.RData")
# Slim Enzyme Data - no stratification
my_EC.ab_table <- otu_table(ec.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.ec.slim <- phyloseq(my_EC.ab_table, my_sample_data)
dat.ec.slim
# Save Phyloseq obj as .RData file 
save(dat.ec.slim, file = "files/Phyloseq_TBC/Enzymes.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   TBC -  Kegg Orthology
#---------------------------------------------------------- 

KOs.abund <- 
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/ko-cpm-named.tsv", 
    col_names = T)
KOs.abund <- 
  KOs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
KOs.abund.slim <- 
  KOs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
KOs.abund <- 
  KOs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
# All KOs 
my_KOs.ab_table <- otu_table(KOs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.KOs <- phyloseq(my_KOs.ab_table, my_sample_data)
dat.KOs
save(dat.KOs, file = "files/Phyloseq_TBC/KOs_PhyloseqObj.RData")
# Slim KOs - no stratification
my_KOs.ab_table.slim <- otu_table(KOs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.KOs.slim <- phyloseq(my_KOs.ab_table.slim, my_sample_data)
dat.KOs.slim
# Save Phyloseq obj as .RData file 
save(dat.KOs.slim, file = "files/Phyloseq_TBC/KOs.slim_PhyloseqObj.RData")


#---------------------------------------------------------- 
#                   TBC -  Gene Ontology
#---------------------------------------------------------- 

GOs.abund <- 
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/go-cpm-named.tsv", 
    col_names = T)
GOs.abund <- 
  GOs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>% 
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
GOs.abund.slim <- GOs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap)) %>% 
  as.data.frame.matrix()
GOs.abund <- 
  GOs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap)) %>% 
  as.data.frame.matrix()
# All GOs 
my_GOs.ab_table <- otu_table(GOs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.GOs <- phyloseq(my_GOs.ab_table, my_sample_data)
dat.GOs
save(dat.GOs, file = "files/Phyloseq_TBC/GOs_PhyloseqObj.RData")
# Slim GOs - no stratification
my_GOs.ab_table.slim <- otu_table(GOs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.GOs.slim <- phyloseq(my_GOs.ab_table.slim, my_sample_data)
dat.GOs.slim
# Save Phyloseq obj as .RData file 
save(dat.GOs.slim, file = "files/Phyloseq_TBC/GOs.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   TBC -  Pfam
#---------------------------------------------------------- 

PFAMs.abund <- 
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/pfam-cpm-named.tsv", 
    col_names = T)
PFAMs.abund <- 
  PFAMs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
PFAMs.abund.slim <- PFAMs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap)) %>% 
  as.data.frame.matrix()
PFAMs.abund <- 
  PFAMs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap)) %>% 
  as.data.frame.matrix()
# All PFAMs 
my_PFAMs.ab_table <- otu_table(PFAMs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.PFAMs <- phyloseq(my_PFAMs.ab_table, my_sample_data)
dat.PFAMs
save(dat.PFAMs, file = "files/Phyloseq_TBC/PFAMs_PhyloseqObj.RData")
# Slim PFAMs - no stratification
my_PFAMs.ab_table.slim <- otu_table(PFAMs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.PFAMs.slim <- phyloseq(my_PFAMs.ab_table.slim, my_sample_data)
dat.PFAMs.slim
# Save Phyloseq obj as .RData file 
save(dat.PFAMs.slim, file = "files/Phyloseq_TBC/PFAMs.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   TBC -  Eggnog
#---------------------------------------------------------- 

EGGNOGs.abund <- 
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/eggnog-cpm.tsv", 
    col_names = T)
EGGNOGs.abund <- 
  EGGNOGs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
EGGNOGs.abund.slim <- EGGNOGs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
EGGNOGs.abund <- 
  EGGNOGs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
# All EGGNOGs 
my_EGGNOGs.ab_table <- otu_table(EGGNOGs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.EGGNOGs <- phyloseq(my_EGGNOGs.ab_table, my_sample_data)
dat.EGGNOGs
save(dat.EGGNOGs, file = "files/Phyloseq_TBC/EGGNOGs_PhyloseqObj.RData")
# Slim EGGNOGs - no stratification
my_EGGNOGs.ab_table.slim <- otu_table(EGGNOGs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.EGGNOGs.slim <- phyloseq(my_EGGNOGs.ab_table.slim, my_sample_data)
dat.EGGNOGs.slim
# Save Phyloseq obj as .RData file 
save(dat.EGGNOGs.slim, file = "files/Phyloseq_TBC/EGGNOGs.slim_PhyloseqObj.RData")


# temporary rename
dat.kingdom.TBC <- dat.kingdom
dat.phylum.TBC <- dat.phylum
dat.class.TBC <- dat.class
dat.order.TBC <- dat.order
dat.family.TBC <- dat.family
dat.genus.TBC <- dat.genus
dat.species.TBC <- dat.species
dat.path.TBC <- dat.path
dat.path.slim.TBC <- dat.path.slim
dat.ec.TBC <- dat.ec
dat.ec.slim.TBC <- dat.ec.slim
dat.KOs.TBC <- dat.KOs
dat.KOs.slim.TBC <- dat.KOs.slim
dat.GOs.TBC <- dat.GOs
dat.GOs.slim.TBC <- dat.GOs.slim
dat.PFAMs.TBC <- dat.PFAMs
dat.PFAMs.slim.TBC <- dat.PFAMs.slim
dat.EGGNOGs.TBC <- dat.EGGNOGs
dat.EGGNOGs.slim.TBC <- dat.EGGNOGs.slim


#---------------------------------------------------------------------------------------
##---------------------------------   RUSH - Cohort  -----------------------------------
#---------------------------------------------------------------------------------------

#---------------------------------------------------------- 
#-                 Metadata Prep - RUSH   
#---------------------------------------------------------- 

metadata_RUSH <- read.csv(file = "files/metadata_phyloseq_RUSH.csv", header= TRUE) 

RUSH_keys <- 
  metadata_RUSH %>% 
  select(donor_id, host_subject_id) %>% 
  mutate(donor_id = as.character(donor_id)) %>% 
  mutate(host_subject_id = as.character(host_subject_id))

RUSH_keymap <- RUSH_keys$host_subject_id
names(RUSH_keymap) <- RUSH_keys$donor_id

reads.RUSH <- load_reads("RUSH")
metadata_RUSH <- left_join(metadata_RUSH, reads.RUSH, by = "donor_id")

rownames(metadata_RUSH) <- metadata_RUSH$donor_id
metadata_RUSH <- as.data.frame(metadata_RUSH)
metadata_RUSH[is.na(metadata_RUSH)] <- "not provided"

#---------------------------------------------------------- 
#-                RUSH -  Taxonomy 
#---------------------------------------------------------- 
met.table <-read_tsv(
  file = "files/RUSH_biobakery_output_slim/metaphlan/merged/metaphlan_taxonomic_profiles.tsv",
  col_names = T)

# Select only species rows from 
bugs.species <- met.table %>% 
  dplyr::rename("taxonomy" = `# taxonomy`) %>% 
  filter(grepl("s__", taxonomy)) %>% 
  filter(!grepl("t__", taxonomy)) %>% 
  column_to_rownames(var = "taxonomy") %>% 
  clean.cols.tax() %>%
  dplyr::select(-contains(negative_controls)) %>%
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap)

dat.species <- metaphlanToPhyloseq_Waldron(
  tax = bugs.species, metadat = metadata_RUSH) %>% 
  subset_samples(study_group != "MSA")

#-------- Species Level Object --------
save(dat.species, file = "files/Phyloseq_RUSH/Species_PhyloseqObj.RData")
#-------- Genus Level Object --------
dat.genus = tax_glom(dat.species, taxrank = "Genus", NArm = F)
taxa_names(dat.genus) <- tax_table(dat.genus)[,6]
save(dat.genus, file = "files/Phyloseq_RUSH/Genus_PhyloseqObj.RData")
#-------- Family Level Object --------
dat.family = tax_glom(dat.species, taxrank = "Family", NArm = F)
taxa_names(dat.family) <- tax_table(dat.family)[,5]
save(dat.family, file = "files/Phyloseq_RUSH/Family_PhyloseqObj.RData")
#-------- Order Level Object --------
dat.order = tax_glom(dat.species, taxrank = "Order", NArm = F)
taxa_names(dat.order) <- tax_table(dat.order)[,4]
save(dat.order, file = "files/Phyloseq_RUSH/Order_PhyloseqObj.RData")
#-------- Class Level Object --------
dat.class = tax_glom(dat.species, taxrank = "Class", NArm = F)
taxa_names(dat.class) <- tax_table(dat.class)[,3]
save(dat.class, file = "files/Phyloseq_RUSH/Class_PhyloseqObj.RData")
#-------- Phylum Level Object --------
dat.phylum = tax_glom(dat.species, taxrank = "Phylum", NArm = F)
taxa_names(dat.phylum) <- tax_table(dat.phylum)[,2]
save(dat.phylum, file = "files/Phyloseq_RUSH/Phylum_PhyloseqObj.RData")
#-------- Kingdom Level Object --------
dat.kingdom = tax_glom(dat.species, taxrank = "Kingdom", NArm = F)
taxa_names(dat.kingdom) <- tax_table(dat.kingdom)[,1]
save(dat.kingdom, file = "files/Phyloseq_RUSH/Kingdom_PhyloseqObj.RData")

#-------------------- -------------------------------------- 
#                   RUSH -  Pathways 
#---------------------------------------------------------- 

path.abund <- 
  read_tsv(file = "files/RUSH_biobakery_output_slim/humann/merged/pathabundance_relab.tsv", col_names = T) 
path.abund <- 
  path.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund() %>% 
  filter(!grepl("UNMAPPED", `# Pathway`)) %>% 
  filter(!grepl("UNINTEGRATED", `# Pathway`))
path.abund.slim <- path.abund %>% 
  filter(!grepl("g__", `# Pathway`)) %>% 
  filter(!grepl("unclassified", `# Pathway`)) %>%
  make_rfriendly_rows(passed_column = "# Pathway") %>% 
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap)
path.abund <- path.abund %>% 
  make_rfriendly_rows(passed_column = "# Pathway") %>% 
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap)
# All Pathway Data
my_pathab_table <- otu_table(path.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.path <- phyloseq(my_pathab_table, my_sample_data)
dat.path
save(dat.path, file = "files/Phyloseq_RUSH/Pathways_PhyloseqObj.RData")
# Slim Pathway Data
my_pathab_table <- otu_table(path.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.path.slim <- phyloseq(my_pathab_table, my_sample_data) %>% 
  subset_samples(study_group != "MSA")
dat.path.slim
save(dat.path.slim, file = "files/Phyloseq_RUSH/Pathways.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   RUSH -  Enzymes 
#---------------------------------------------------------- 

ec.abund <- 
  read_tsv(
    "files/RUSH_biobakery_output_slim/humann/merged/ecs_relab.tsv", 
    col_names = T)
ec.abund <- 
  ec.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_RPK()
ec.abund.slim <- 
  ec.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>% 
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap)
ec.abund <- 
  ec.abund %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>% 
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap)
# All Enzyme Data
my_EC.ab_table <- otu_table(ec.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.ec <- phyloseq(my_EC.ab_table, my_sample_data)
dat.ec
save(dat.ec, file = "files/Phyloseq_RUSH/Enzymes_PhyloseqObj.RData")
# Slim Enzyme Data - no stratification
my_EC.ab_table <- otu_table(ec.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.ec.slim <- phyloseq(my_EC.ab_table, my_sample_data) %>% 
  subset_samples(study_group != "MSA")
dat.ec.slim
# Save Phyloseq obj as .RData file 
save(dat.ec.slim, file = "files/Phyloseq_RUSH/Enzymes.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   RUSH -  Kegg Orthology
#---------------------------------------------------------- 

KOs.abund <- 
  read_tsv(
    "files/RUSH_biobakery_output_slim/humann/merged/ko-cpm-named.tsv", 
    col_names = T)
KOs.abund <- 
  KOs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
KOs.abund.slim <- 
  KOs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap)
KOs.abund <- 
  KOs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap)
# All KOs 
my_KOs.ab_table <- otu_table(KOs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.KOs <- phyloseq(my_KOs.ab_table, my_sample_data)
dat.KOs
save(dat.KOs, file = "files/Phyloseq_RUSH/KOs_PhyloseqObj.RData")
# Slim KOs - no stratification
my_KOs.ab_table.slim <- otu_table(KOs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.KOs.slim <- phyloseq(my_KOs.ab_table.slim, my_sample_data) %>% 
  subset_samples(study_group != "MSA")
dat.KOs.slim
# Save Phyloseq obj as .RData file 
save(dat.KOs.slim, file = "files/Phyloseq_RUSH/KOs.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   RUSH -  Gene Ontology
#---------------------------------------------------------- 

GOs.abund <- 
  read_tsv(
    "files/RUSH_biobakery_output_slim/humann/merged/go-cpm-named.tsv", 
    col_names = T)
GOs.abund <- 
  GOs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
GOs.abund.slim <- GOs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap) %>% 
  as.data.frame.matrix()
GOs.abund <- 
  GOs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap) %>% 
  as.data.frame.matrix()
# All GOs 
my_GOs.ab_table <- otu_table(GOs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.GOs <- phyloseq(my_GOs.ab_table, my_sample_data)
dat.GOs
save(dat.GOs, file = "files/Phyloseq_RUSH/GOs_PhyloseqObj.RData")
# Slim GOs - no stratification
my_GOs.ab_table.slim <- otu_table(GOs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.GOs.slim <- phyloseq(my_GOs.ab_table.slim, my_sample_data) %>% 
  subset_samples(study_group != "MSA")
dat.GOs.slim
# Save Phyloseq obj as .RData file 
save(dat.GOs.slim, file = "files/Phyloseq_RUSH/GOs.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   RUSH -  Pfam
#---------------------------------------------------------- 

PFAMs.abund <- 
  read_tsv(
    "files/RUSH_biobakery_output_slim/humann/merged/pfam-cpm-named.tsv", 
    col_names = T)
PFAMs.abund <- 
  PFAMs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
PFAMs.abund.slim <- PFAMs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap) %>% 
  as.data.frame.matrix()
PFAMs.abund <- 
  PFAMs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap) %>% 
  as.data.frame.matrix()
# All PFAMs 
my_PFAMs.ab_table <- otu_table(PFAMs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.PFAMs <- phyloseq(my_PFAMs.ab_table, my_sample_data)
dat.PFAMs
save(dat.PFAMs, file = "files/Phyloseq_RUSH/PFAMs_PhyloseqObj.RData")
# Slim PFAMs - no stratification
my_PFAMs.ab_table.slim <- otu_table(PFAMs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.PFAMs.slim <- phyloseq(my_PFAMs.ab_table.slim, my_sample_data) %>% 
  subset_samples(study_group != "MSA")
dat.PFAMs.slim
# Save Phyloseq obj as .RData file 
save(dat.PFAMs.slim, file = "files/Phyloseq_RUSH/PFAMs.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   RUSH -  Eggnog
#---------------------------------------------------------- 

EGGNOGs.abund <- 
  read_tsv(
    "files/RUSH_biobakery_output_slim/humann/merged/eggnog-cpm.tsv", 
    col_names = T)
EGGNOGs.abund <- 
  EGGNOGs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
EGGNOGs.abund.slim <- EGGNOGs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap)
EGGNOGs.abund <- 
  EGGNOGs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap)
# All EGGNOGs 
my_EGGNOGs.ab_table <- otu_table(EGGNOGs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.EGGNOGs <- phyloseq(my_EGGNOGs.ab_table, my_sample_data)
dat.EGGNOGs
save(dat.EGGNOGs, file = "files/Phyloseq_RUSH/EGGNOGs_PhyloseqObj.RData")
# Slim EGGNOGs - no stratification
my_EGGNOGs.ab_table.slim <- otu_table(EGGNOGs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.EGGNOGs.slim <- phyloseq(my_EGGNOGs.ab_table.slim, my_sample_data) %>% 
  subset_samples(study_group != "MSA")
dat.EGGNOGs.slim
# Save Phyloseq obj as .RData file 
save(dat.EGGNOGs.slim, file = "files/Phyloseq_RUSH/EGGNOGs.slim_PhyloseqObj.RData")


# temporary rename
dat.kingdom.RUSH <- dat.kingdom
dat.phylum.RUSH <- dat.phylum
dat.class.RUSH <- dat.class
dat.order.RUSH <- dat.order
dat.family.RUSH <- dat.family
dat.genus.RUSH <- dat.genus
dat.species.RUSH <- dat.species

dat.path.RUSH <- dat.path
dat.path.slim.RUSH <- dat.path.slim
dat.ec.RUSH <- dat.ec
dat.ec.slim.RUSH <- dat.ec.slim
dat.KOs.RUSH <- dat.KOs
dat.KOs.slim.RUSH <- dat.KOs.slim
dat.GOs.RUSH <- dat.GOs
dat.GOs.slim.RUSH <- dat.GOs.slim
dat.PFAMs.RUSH <- dat.PFAMs
dat.PFAMs.slim.RUSH <- dat.PFAMs.slim
dat.EGGNOGs.RUSH <- dat.EGGNOGs
dat.EGGNOGs.slim.RUSH <- dat.EGGNOGs.slim



#---------------------------------------------------------------------------------------
##--------------------------------   Merge - Cohorts  ----------------------------------
#---------------------------------------------------------------------------------------



dat.kingdom <- merge_phyloseq(dat.kingdom.TBC, dat.kingdom.RUSH)
dat.phylum <- merge_phyloseq(dat.phylum.TBC, dat.phylum.RUSH)
dat.class <- merge_phyloseq(dat.class.TBC, dat.class.RUSH)
dat.order <- merge_phyloseq(dat.order.TBC, dat.order.RUSH)
dat.family <- merge_phyloseq(dat.family.TBC, dat.family.RUSH)
dat.genus <- merge_phyloseq(dat.genus.TBC, dat.genus.RUSH)
dat.species <- merge_phyloseq(dat.species.TBC, dat.species.RUSH)

dat.path <- merge_phyloseq(dat.path.TBC, dat.path.RUSH)
dat.path.slim <- merge_phyloseq(dat.path.slim.TBC, dat.path.slim.RUSH)
dat.ec <- merge_phyloseq(dat.ec.TBC, dat.ec.RUSH)
dat.ec.slim <- merge_phyloseq(dat.ec.slim.TBC, dat.ec.slim.RUSH)
dat.KOs <- merge_phyloseq(dat.KOs.TBC, dat.KOs.RUSH)
dat.KOs.slim <- merge_phyloseq(dat.KOs.slim.TBC, dat.KOs.slim.RUSH)
dat.GOs <- merge_phyloseq(dat.GOs.TBC, dat.GOs.RUSH)
dat.GOs.slim <- merge_phyloseq(dat.GOs.slim.TBC, dat.GOs.slim.RUSH)
dat.PFAMs <- merge_phyloseq(dat.PFAMs.TBC, dat.PFAMs.RUSH)
dat.PFAMs.slim <- merge_phyloseq(dat.PFAMs.slim.TBC, dat.PFAMs.slim.RUSH)
dat.EGGNOGs <- merge_phyloseq(dat.EGGNOGs.TBC, dat.EGGNOGs.RUSH)
dat.EGGNOGs.slim <- merge_phyloseq(dat.EGGNOGs.slim.TBC, dat.EGGNOGs.slim.RUSH)




save(dat.species, file = "files/Phyloseq_Merged/Species_PhyloseqObj.RData")
save(dat.genus, file = "files/Phyloseq_Merged/Genus_PhyloseqObj.RData")
save(dat.family, file = "files/Phyloseq_Merged/Family_PhyloseqObj.RData")
save(dat.order, file = "files/Phyloseq_Merged/Order_PhyloseqObj.RData")
save(dat.class, file = "files/Phyloseq_Merged/Class_PhyloseqObj.RData")
save(dat.phylum, file = "files/Phyloseq_Merged/Phylum_PhyloseqObj.RData")
save(dat.kingdom, file = "files/Phyloseq_Merged/Kingdom_PhyloseqObj.RData")
save(dat.path, file = "files/Phyloseq_Merged/Pathways_PhyloseqObj.RData")
save(dat.path.slim, file = "files/Phyloseq_Merged/Pathways.slim_PhyloseqObj.RData")
save(dat.ec, file = "files/Phyloseq_Merged/Enzymes_PhyloseqObj.RData")
save(dat.ec.slim, file = "files/Phyloseq_Merged/Enzymes.slim_PhyloseqObj.RData")
save(dat.KOs, file = "files/Phyloseq_Merged/KOs_PhyloseqObj.RData")
save(dat.KOs.slim, file = "files/Phyloseq_Merged/KOs.slim_PhyloseqObj.RData")
save(dat.GOs, file = "files/Phyloseq_Merged/GOs_PhyloseqObj.RData")
save(dat.GOs.slim, file = "files/Phyloseq_Merged/GOs.slim_PhyloseqObj.RData")
save(dat.PFAMs, file = "files/Phyloseq_Merged/PFAMs_PhyloseqObj.RData")
save(dat.PFAMs.slim, file = "files/Phyloseq_Merged/PFAMs.slim_PhyloseqObj.RData")
save(dat.EGGNOGs, file = "files/Phyloseq_Merged/EGGNOGs_PhyloseqObj.RData")
save(dat.EGGNOGs.slim, file = "files/Phyloseq_Merged/EGGNOGs.slim_PhyloseqObj.RData")






#-------------------------------------------------------------------------------------------
### Create list for objects
Phylo_Objects <- vector(mode="list", length=19)
names(Phylo_Objects) <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom",
                          "Pathways", "Pathways.slim",
                          "Enzymes", "Enzymes.slim",
                          "KOs", "KOs.slim", 
                          "dat.GOs", "dat.GOs.slim",
                          "dat.PFAMs", "dat.PFAMs.slim",
                          "dat.EGGNOGs", "dat.EGGNOGs.slim")

Phylo_Objects$Species <- dat.species; Phylo_Objects$Genus <- dat.genus;
Phylo_Objects$Family <- dat.family; Phylo_Objects$Order <- dat.order; Phylo_Objects$Class <- dat.class;
Phylo_Objects$Phylum <- print(dat.phylum); Phylo_Objects$Kingdom <- dat.kingdom
Phylo_Objects$Pathways <- dat.path; Phylo_Objects$Pathways.slim <- dat.path.slim;
Phylo_Objects$Enzymes <- dat.ec; Phylo_Objects$Enzymes.slim <- dat.ec.slim;
Phylo_Objects$KOs <- dat.KOs; Phylo_Objects$KOs.slim <- dat.KOs.slim
Phylo_Objects$GOs <- dat.GOs; Phylo_Objects$GOs.slim <- dat.GOs.slim
Phylo_Objects$PFAMs <- dat.PFAMs; Phylo_Objects$PFAMs.slim <- dat.PFAMs.slim
Phylo_Objects$EGGNOGs <- dat.EGGNOGs; Phylo_Objects$EGGNOGs.slim <- dat.EGGNOGs.slim
save(Phylo_Objects, file = "files/Phyloseq_Merged/PhyloseqObj.RData")
#-------------------------------------------------------------------------------------------

