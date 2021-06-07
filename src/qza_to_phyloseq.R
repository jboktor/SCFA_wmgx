# Create phyloseq objects
rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
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

rownames(metadata_df) <- metadata_df$donor_id
metadata_df <- as.data.frame(metadata_df)

#---------------------------------------------------------- 
#-                     Taxonomy 
#---------------------------------------------------------- 

#-------- Species Level Object --------

species.table.df <- otu_table(import_biom("files/03_species/109256_species.biom")) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame() 

species.table <- species.table.df %>% 
  otu_table(taxa_are_rows = T)

species_taxa_strat <- 
  species.table.df %>%  
  rownames_to_column() %>% select(rowname) %>% 
  separate(rowname, c("Kingdom", "Phylum", "Class", "Order",
                      "Family", "Genus", "Species"), 
           sep = ";", remove = F) %>%  column_to_rownames() %>% 
  as.matrix()

dat.species <- 
  phyloseq(
    features = species.table,
    metadata =sample_data(metadata_df), 
    tax_table(species_taxa_strat))
newnames <- gsub("^.*s__", "", taxa_names(dat.species))
newnames <- gsub(" ", "_", newnames)
taxa_names(dat.species) <- newnames
dat.species


# #-------- Genus Level Object --------
genus.table.df <- otu_table(import_biom("files/04_genus/109257_genus.biom")) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame() 

genus.table <- genus.table.df %>% 
  otu_table(taxa_are_rows = T)

genus_taxa_strat <- 
  genus.table.df %>%  
  rownames_to_column() %>% select(rowname) %>% 
  separate(rowname, c("Kingdom", "Phylum", "Class", "Order",
                      "Family", "Genus"), 
           sep = ";", remove = F) %>%  column_to_rownames() %>% 
  as.matrix()

dat.genus <- 
  phyloseq(
    features = genus.table,
    metadata =sample_data(metadata_df),
    tax_table(genus_taxa_strat))
newnames <- gsub("^.*g__", "", taxa_names(dat.genus))
newnames <- gsub(" ", "_", newnames)
taxa_names(dat.genus) <- newnames
dat.genus

# #-------- Phylum Level Object --------
phylum.table.df <- otu_table(import_biom("files/05_phylum/109254_phylum.biom")) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame()

phylum.table <- phylum.table.df %>% 
  otu_table(taxa_are_rows = T)

phylum_taxa_strat <- 
  phylum.table.df %>%  
  rownames_to_column() %>% select(rowname) %>% 
  separate(rowname, c("Kingdom", "Phylum"), 
           sep = ";", remove = F) %>%  column_to_rownames() %>% 
  as.matrix()

dat.phylum<- 
  phyloseq(
    features = phylum.table,
    metadata =sample_data(metadata_df),
    tax_table(phylum_taxa_strat))
newnames <- gsub("^.*p__", "", taxa_names(dat.phylum))
newnames <- gsub(" ", "_", newnames)
taxa_names(dat.phylum) <- newnames
dat.phylum

# # Rename Species after regrouping
# #-------- Family Level Object --------
# dat.family = tax_glom(dat.species, taxrank = "Rank5", NArm = F)
# taxa_names(dat.family) <- tax_table(dat.family)[,5]
# #-------- Order Level Object --------
# dat.order = tax_glom(dat.species, taxrank = "Rank4", NArm = F)
# taxa_names(dat.order) <- tax_table(dat.order)[,4]
# #-------- Class Level Object --------
# dat.class = tax_glom(dat.genus, taxrank = "Class", NArm = F)
# tst <- taxa_names(dat.class) <- tax_table(dat.class)[,3]
# #-------- Kingdom Level Object --------
# dat.kingdom = tax_glom(dat.phylum, taxrank = "Kingdom", NArm = F)
# taxa_names(dat.kingdom) <- tax_table(dat.kingdom)[,1]


#----------------------------------------------------------
#                        Functional
#----------------------------------------------------------

#-------- MetaCyc proteins  --------
MetaCycProteins.table.df <- read_tsv(
  file = "files/02_per_genome_predictions/woltka_subgrouping/protein.tsv", 
  skip = 1) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  column_to_rownames(var = "#OTU ID") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame()

dat.MCProteins <- 
  phyloseq(
    features = otu_table(MetaCycProteins.table.df, taxa_are_rows = T),
    metadata =sample_data(metadata_df))
dat.MCProteins

#-------- MetaCyc Rxns  --------
MetaCycRxns.table.df <- read_tsv(
  file = "files/02_per_genome_predictions/woltka_subgrouping/reaction.tsv", 
  skip = 1) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  column_to_rownames(var = "#OTU ID") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame()

dat.Rxns <- 
  phyloseq(
    features = otu_table(MetaCycRxns.table.df, taxa_are_rows = T),
    metadata =sample_data(metadata_df))
dat.Rxns


#-------- enzrxns  --------
ecs.table.df <- read_tsv(
  file = "files/02_per_genome_predictions/woltka_subgrouping/enzrxn.tsv", 
  skip = 1) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  column_to_rownames(var = "#OTU ID") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame()

dat.Enzrxn<- 
  phyloseq(
    features = otu_table(ecs.table.df, taxa_are_rows = T),
    metadata =sample_data(metadata_df))
dat.Enzrxn

#-------- MetaCyc Pathways  --------
MetaCycPathways.table.df <- read_tsv(
  file = "files/02_per_genome_predictions/woltka_subgrouping/pathway.tsv", 
  skip = 1) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  column_to_rownames(var = "#OTU ID") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame()

dat.Pathways <- 
  phyloseq(
    features = otu_table(MetaCycPathways.table.df, taxa_are_rows = T),
    metadata =sample_data(metadata_df))
dat.Pathways

#-------- MetaCyc Super Pathways  --------
MetaCycSuperPathways.table.df <- read_tsv(
  file = "files/02_per_genome_predictions/woltka_subgrouping/superpathway.tsv", 
  skip = 1) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  column_to_rownames(var = "#OTU ID") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame()

dat.Superpathways <- 
  phyloseq(
    features = otu_table(MetaCycSuperPathways.table.df, taxa_are_rows = T),
    metadata =sample_data(metadata_df))
dat.Superpathways

#-------- GOs  --------
gos.table.df <- read_tsv(
  file = "files/02_per_genome_predictions/woltka_subgrouping/go.tsv", 
  skip = 1) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  column_to_rownames(var = "#OTU ID") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame()

dat.GOs <- 
  phyloseq(
    features = otu_table(gos.table.df, taxa_are_rows = T),
    metadata =sample_data(metadata_df))
dat.GOs

#-------- KOs  --------
kos.table.df <- read_tsv(
  file = "files/02_per_genome_predictions/woltka_subgrouping/ko.tsv", 
  skip = 1) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  column_to_rownames(var = "#OTU ID") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame()

dat.KOs <- 
  phyloseq(
    features = otu_table(kos.table.df, taxa_are_rows = T),
    metadata =sample_data(metadata_df))
dat.KOs

#-------- eggnogs  --------
eggNOGs.table.df <- read_tsv(
  file = "files/02_per_genome_predictions/woltka_subgrouping/eggnog.tsv", 
  skip = 1) %>% 
  as.data.frame() %>% 
  clean.cols.qiime() %>% 
  clean.cols.qiime2() %>% 
  dplyr::select(!contains("BLANK")) %>% 
  column_to_rownames(var = "#OTU ID") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("tempID") %>%
  arrange(tempID) %>%
  column_to_rownames(var = "tempID") %>%
  t() %>% as.data.frame()

dat.eggNOGs <- 
  phyloseq(
    features = otu_table(eggNOGs.table.df, taxa_are_rows = T),
    metadata =sample_data(metadata_df))
dat.eggNOGs


# Save data
save(dat.species, file = "files/Phyloseq_objects_Woltka/Species_PhyloseqObj.RData")
save(dat.genus, file = "files/Phyloseq_objects_Woltka/Genus_PhyloseqObj.RData")
save(dat.phylum, file = "files/Phyloseq_objects_Woltka/Phylum_PhyloseqObj.RData")
save(dat.MCProteins, file = "files/Phyloseq_objects_Woltka/MCProteins_PhyloseqObj.RData")
save(dat.Rxns, file = "files/Phyloseq_objects_Woltka/Rxns_PhyloseqObj.RData")
save(dat.Enzrxn, file = "files/Phyloseq_objects_Woltka/Enzrxn_PhyloseqObj.RData")
save(dat.Pathways, file = "files/Phyloseq_objects_Woltka/Pathways_PhyloseqObj.RData")
save(dat.Superpathways, file = "files/Phyloseq_objects_Woltka/Superpathways_PhyloseqObj.RData")
save(dat.GOs, file = "files/Phyloseq_objects_Woltka/GOs_PhyloseqObj.RData")
save(dat.KOs, file = "files/Phyloseq_objects_Woltka/KOs_PhyloseqObj.RData")
save(dat.eggNOGs, file = "files/Phyloseq_objects_Woltka/eggNOGs_PhyloseqObj.RData")

phylo_names <- c(
  "Species",
  "Genus",
  "Phylum",
  "MCProteins",
  "Enzrxn",
  "Reactions",
  "Superpathways",
  "Pathways",
  "KOs",
  "GOs",
  "eggNOGs"
)

phylo_dats <- 
  list(dat.species,
    dat.genus,
    dat.phylum,
    dat.MCProteins,
    dat.Enzrxn,
    dat.Rxns,
    dat.Superpathways,
    dat.Pathways,
    dat.KOs,
    dat.GOs,
    dat.eggNOGs)

names(phylo_dats) <- phylo_names
save(phylo_dats, file = "files/Phyloseq_objects_Woltka.RData")

# Rarefy data
phylo_dats_rare <- list()
for (obj in names(phylo_dats)){
  print(obj)
  phylo_dats_rare[[obj]] <- rarefunc(phylo_dats[[obj]])
}
save(phylo_dats_rare, file = "files/Phyloseq_objects_Woltka_Rarefied.RData")


