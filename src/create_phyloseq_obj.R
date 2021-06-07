# Create phyloseq objects

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/metaphlanToPhyloseq_Waldron.R")

negative_controls <- c("BLANK")

#---------------------------------------------------------- 
#-                 metadata prep     
#---------------------------------------------------------- 

motor_testing_vars <- c(
  "host_weight",
  "pole_test_s",
  "sticker_test_s",
  "x_of_steps",
  "errors_step",
  "beam_time_to_cross",
  "wire_hang_s",
  "hindlimb_score"
)


metadata_df <- 
  read.csv(file = "files/metadata.csv", header= TRUE) %>% 
  janitor::clean_names() %>% 
  dplyr::rename(donor_id = tube_id) %>% 
  dplyr::filter(str_detect(donor_id, negative_controls, negate = T)) %>%
  dplyr::select(host_weight, collection_timestamp, donor_id:description) %>% 
  dplyr::mutate(
    diet = str_replace(diet, "Acetate", "Prebiotic"),
    treatment_group = str_replace(treatment_group, "Acetate", "Prebiotic")) %>% 
  mutate(across(all_of(c(motor_testing_vars)), as.character),
         across(all_of(motor_testing_vars), as.numeric)) %>% 
  mutate_if(is.character, as.factor) 

reads <- load_reads()
metadata_df <- left_join(metadata_df, reads, by = "donor_id")
rownames(metadata_df) <- metadata_df$donor_id
metadata_df <- as.data.frame(metadata_df)



#---------------------------------------------------------- 
#-                     Taxonomy 
#---------------------------------------------------------- 
met.table <-
  read_tsv(file = "files/SCFA_biobakery_output_slim/metaphlan/merged/metaphlan_taxonomic_profiles.tsv",
           col_names = T)


# Select only species rows from 
bugs.species <- 
  met.table %>% 
  dplyr::rename("taxonomy" = `# taxonomy`) %>% 
  filter(grepl("s__", taxonomy)) %>% 
  filter(!grepl("t__", taxonomy)) %>% 
  column_to_rownames(var = "taxonomy") %>% 
  clean.cols.tax() %>%
  dplyr::select(-contains(negative_controls)) %>%
  trim_cols() %>%
  trim_underscore()

#-------- Species Level Object --------
dat.species <- metaphlanToPhyloseq_Waldron(
  tax = bugs.species, metadat = metadata_df)
#-------- Genus Level Object --------
dat.genus = tax_glom(dat.species, taxrank = "Genus", NArm = F)
taxa_names(dat.genus) <- tax_table(dat.genus)[,6]
#-------- Family Level Object --------
dat.family = tax_glom(dat.species, taxrank = "Family", NArm = F)
taxa_names(dat.family) <- tax_table(dat.family)[,5]
#-------- Order Level Object --------
dat.order = tax_glom(dat.species, taxrank = "Order", NArm = F)
taxa_names(dat.order) <- tax_table(dat.order)[,4]
#-------- Class Level Object --------
dat.class = tax_glom(dat.species, taxrank = "Class", NArm = F)
taxa_names(dat.class) <- tax_table(dat.class)[,3]
#-------- Phylum Level Object --------
dat.phylum = tax_glom(dat.species, taxrank = "Phylum", NArm = F)
taxa_names(dat.phylum) <- tax_table(dat.phylum)[,2]
#-------- Kingdom Level Object --------
dat.kingdom = tax_glom(dat.species, taxrank = "Kingdom", NArm = F)
taxa_names(dat.kingdom) <- tax_table(dat.kingdom)[,1]

#-------------------- -------------------------------------- 
#                    -  Pathways 
#---------------------------------------------------------- 

path.abund <- 
  read_tsv(file = "files/SCFA_biobakery_output_slim/humann/merged/pathabundance_relab.tsv", col_names = T)
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
  trim_cols() %>%
  trim_underscore()
path.abund <- path.abund %>% 
  make_rfriendly_rows(passed_column = "# Pathway") %>% 
  trim_cols() %>%
  trim_underscore()
# All Pathway Data
my_pathab_table <- otu_table(path.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.path <- phyloseq(my_pathab_table, my_sample_data)
dat.path
# Slim Pathway Data
my_pathab_table <- otu_table(path.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.path.slim <- phyloseq(my_pathab_table, my_sample_data)
dat.path.slim

#---------------------------------------------------------- 
#                        Enzymes 
#---------------------------------------------------------- 

ec.abund <- 
  read_tsv(
    "files/SCFA_biobakery_output_slim/humann/merged/ecs_relab.tsv", 
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
  trim_cols() %>%
  trim_underscore()
ec.abund <- 
  ec.abund %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>% 
  trim_cols() %>%
  trim_underscore()
# All Enzyme Data
my_EC.ab_table <- otu_table(ec.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.ec <- phyloseq(my_EC.ab_table, my_sample_data)
dat.ec
# Slim Enzyme Data - no stratification
my_EC.ab_table <- otu_table(ec.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.ec.slim <- phyloseq(my_EC.ab_table, my_sample_data)
dat.ec.slim

#---------------------------------------------------------- 
#                     Kegg Orthology
#---------------------------------------------------------- 

KOs.abund <- 
  read_tsv(
    "files/SCFA_biobakery_output_slim/humann/merged/ko-cpm-named.tsv", 
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
  trim_cols() %>%
  trim_underscore()
KOs.abund <- 
  KOs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols() %>%
  trim_underscore()
# All KOs 
my_KOs.ab_table <- otu_table(KOs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.KOs <- phyloseq(my_KOs.ab_table, my_sample_data)
dat.KOs
# Slim KOs - no stratification
my_KOs.ab_table.slim <- otu_table(KOs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.KOs.slim <- phyloseq(my_KOs.ab_table.slim, my_sample_data)
dat.KOs.slim

#---------------------------------------------------------- 
#                     Gene Ontology
#---------------------------------------------------------- 

GOs.abund <- 
  read_tsv(
    "files/SCFA_biobakery_output_slim/humann/merged/go-cpm-named.tsv", 
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
  trim_cols() %>%
  trim_underscore()
  # as.data.frame.matrix()
GOs.abund <- 
  GOs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols() %>%
  trim_underscore()
  # as.data.frame.matrix()
# All GOs 
my_GOs.ab_table <- otu_table(GOs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.GOs <- phyloseq(my_GOs.ab_table, my_sample_data)
dat.GOs

# Slim GOs - no stratification
my_GOs.ab_table.slim <- otu_table(GOs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.GOs.slim <- phyloseq(my_GOs.ab_table.slim, my_sample_data)
dat.GOs.slim

#---------------------------------------------------------- 
#                        Pfams
#---------------------------------------------------------- 

PFAMs.abund <- 
  read_tsv(
    "files/SCFA_biobakery_output_slim/humann/merged/pfam-cpm-named.tsv", 
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
  trim_cols() %>%
  trim_underscore()
  # as.data.frame.matrix()
PFAMs.abund <- 
  PFAMs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols() %>%
  trim_underscore()
  # as.data.frame.matrix()
# All PFAMs 
my_PFAMs.ab_table <- otu_table(PFAMs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.PFAMs <- phyloseq(my_PFAMs.ab_table, my_sample_data)
dat.PFAMs

# Slim PFAMs - no stratification
my_PFAMs.ab_table.slim <- otu_table(PFAMs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.PFAMs.slim <- phyloseq(my_PFAMs.ab_table.slim, my_sample_data)
dat.PFAMs.slim

#---------------------------------------------------------- 
#                      eggNOGs
#---------------------------------------------------------- 

eggNOGs.abund <- 
  read_tsv(
    "files/SCFA_biobakery_output_slim/humann/merged/eggnog-cpm.tsv", 
    col_names = T)
eggNOGs.abund <- 
  eggNOGs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
eggNOGs.abund.slim <- eggNOGs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols() %>%
  trim_underscore()
eggNOGs.abund <- 
  eggNOGs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols() %>%
  trim_underscore()
# All eggNOGs 
my_eggNOGs.ab_table <- otu_table(eggNOGs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.eggNOGs <- phyloseq(my_eggNOGs.ab_table, my_sample_data)
dat.eggNOGs

# Slim eggNOGs - no stratification
my_eggNOGs.ab_table.slim <- otu_table(eggNOGs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.eggNOGs.slim <- phyloseq(my_eggNOGs.ab_table.slim, my_sample_data)
dat.eggNOGs.slim


#-------------------------------------------------------------------------------
########                         Save Objects                          ########          
#-------------------------------------------------------------------------------


save(dat.species, file = "files/Phyloseq_objects/Species_PhyloseqObj.RData")
save(dat.genus, file = "files/Phyloseq_objects/Genus_PhyloseqObj.RData")
save(dat.family, file = "files/Phyloseq_objects/Family_PhyloseqObj.RData")
save(dat.order, file = "files/Phyloseq_objects/Order_PhyloseqObj.RData")
save(dat.class, file = "files/Phyloseq_objects/Class_PhyloseqObj.RData")
save(dat.phylum, file = "files/Phyloseq_objects/Phylum_PhyloseqObj.RData")
save(dat.kingdom, file = "files/Phyloseq_objects/Kingdom_PhyloseqObj.RData")
save(dat.path, file = "files/Phyloseq_objects/Pathways_PhyloseqObj.RData")
save(dat.path.slim, file = "files/Phyloseq_objects/Pathways.slim_PhyloseqObj.RData")
save(dat.ec, file = "files/Phyloseq_objects/Enzymes_PhyloseqObj.RData")
save(dat.ec.slim, file = "files/Phyloseq_objects/Enzymes.slim_PhyloseqObj.RData")
save(dat.KOs, file = "files/Phyloseq_objects/KOs_PhyloseqObj.RData")
save(dat.KOs.slim, file = "files/Phyloseq_objects/KOs.slim_PhyloseqObj.RData")
save(dat.GOs, file = "files/Phyloseq_objects/GOs_PhyloseqObj.RData")
save(dat.GOs.slim, file = "files/Phyloseq_objects/GOs.slim_PhyloseqObj.RData")
save(dat.PFAMs, file = "files/Phyloseq_objects/PFAMs_PhyloseqObj.RData")
save(dat.PFAMs.slim, file = "files/Phyloseq_objects/PFAMs.slim_PhyloseqObj.RData")
save(dat.eggNOGs, file = "files/Phyloseq_objects/eggNOGs_PhyloseqObj.RData")
save(dat.eggNOGs.slim, file = "files/Phyloseq_objects/eggNOGs.slim_PhyloseqObj.RData")


#-------------------------------------------------------------------------------------------
# ### Create list for objects
# Phylo_Objects <- vector(mode="list", length=19)
# names(Phylo_Objects) <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom",
#                           "Pathways", "Pathways.slim",
#                           "Enzymes", "Enzymes.slim",
#                           "KOs", "KOs.slim", 
#                           "dat.GOs", "dat.GOs.slim",
#                           "dat.PFAMs", "dat.PFAMs.slim",
#                           "dat.eggNOGs", "dat.eggNOGs.slim")
# 
# Phylo_Objects$Species <- dat.species; Phylo_Objects$Genus <- dat.genus;
# Phylo_Objects$Family <- dat.family; Phylo_Objects$Order <- dat.order; Phylo_Objects$Class <- dat.class;
# Phylo_Objects$Phylum <- print(dat.phylum); Phylo_Objects$Kingdom <- dat.kingdom
# Phylo_Objects$Pathways <- dat.path; Phylo_Objects$Pathways.slim <- dat.path.slim;
# Phylo_Objects$Enzymes <- dat.ec; Phylo_Objects$Enzymes.slim <- dat.ec.slim;
# Phylo_Objects$KOs <- dat.KOs; Phylo_Objects$KOs.slim <- dat.KOs.slim
# Phylo_Objects$GOs <- dat.GOs; Phylo_Objects$GOs.slim <- dat.GOs.slim
# Phylo_Objects$PFAMs <- dat.PFAMs; Phylo_Objects$PFAMs.slim <- dat.PFAMs.slim
# Phylo_Objects$eggNOGs <- dat.eggNOGs; Phylo_Objects$eggNOGs.slim <- dat.eggNOGs.slim
# save(Phylo_Objects, file = "files/Phyloseq_objects/PhyloseqObj.RData")
#-------------------------------------------------------------------------------------------

