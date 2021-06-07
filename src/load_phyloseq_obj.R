# Load Phyloseq Objects

load_all_cohorts <- function(){
  #### Taxa
  load("files/Phyloseq_objects/Species_PhyloseqObj.RData")
  cat("Loading: Species \n"); assign("dat.species", dat.species, envir = .GlobalEnv)
  load("files/Phyloseq_objects/Genus_PhyloseqObj.RData")
  cat("Loading: Genera \n"); assign("dat.genus", dat.genus, envir = .GlobalEnv)
  load("files/Phyloseq_objects/Family_PhyloseqObj.RData")
  cat("Loading: Families \n"); assign("dat.family", dat.family, envir = .GlobalEnv)
  load("files/Phyloseq_objects/Order_PhyloseqObj.RData")
  cat("Loading: Orders \n"); assign("dat.order", dat.order, envir = .GlobalEnv)
  load("files/Phyloseq_objects/Class_PhyloseqObj.RData")
  cat("Loading: Classes \n"); assign("dat.class", dat.class, envir = .GlobalEnv)
  load("files/Phyloseq_objects/Phylum_PhyloseqObj.RData")
  cat("Loading: Phyla \n"); assign("dat.phylum", dat.phylum, envir = .GlobalEnv)
  load("files/Phyloseq_objects/Kingdom_PhyloseqObj.RData")
  cat("Loading: Kindoms \n"); assign("dat.kingdom", dat.kingdom, envir = .GlobalEnv)
  ## Function - Pathways/Enzymes/Genes
  load("files/Phyloseq_objects/Pathways_PhyloseqObj.RData")
  cat("Loading: MetaCyc Pathways: Stratified \n"); assign("dat.path", dat.path, envir = .GlobalEnv)
  load("files/Phyloseq_objects/Pathways.slim_PhyloseqObj.RData")
  cat("Loading: MetaCyc Pathways \n"); assign("dat.path.slim", dat.path.slim, envir = .GlobalEnv)
  load("files/Phyloseq_objects/Enzymes_PhyloseqObj.RData")
  cat("Loading: Enzymes: Stratified \n"); assign("dat.ec", dat.ec, envir = .GlobalEnv)
  load("files/Phyloseq_objects/Enzymes.slim_PhyloseqObj.RData")
  cat("Loading: Enzymes \n"); assign("dat.ec.slim", dat.ec.slim, envir = .GlobalEnv)
  load("files/Phyloseq_objects/KOs_PhyloseqObj.RData")
  cat("Loading: Kegg Orthology: Stratified \n"); assign("dat.KOs", dat.KOs, envir = .GlobalEnv)
  load("files/Phyloseq_objects/KOs.slim_PhyloseqObj.RData")
  cat("Loading: Kegg Orthology \n"); assign("dat.KOs.slim", dat.KOs.slim, envir = .GlobalEnv)
  load("files/Phyloseq_objects/GOs_PhyloseqObj.RData")
  cat("Loading: Gene Ontology: Stratified \n"); assign("dat.GOs", dat.GOs, envir = .GlobalEnv)
  load("files/Phyloseq_objects/GOs.slim_PhyloseqObj.RData")
  cat("Loading: Gene Ontology \n"); assign("dat.GOs.slim", dat.GOs.slim, envir = .GlobalEnv)
  load("files/Phyloseq_objects/PFAMs_PhyloseqObj.RData")
  cat("Loading: PFAMs: Stratified \n"); assign("dat.PFAMs", dat.PFAMs, envir = .GlobalEnv)
  load("files/Phyloseq_objects/PFAMs.slim_PhyloseqObj.RData")
  cat("Loading: PFAMs \n"); assign("dat.PFAMs.slim", dat.PFAMs.slim, envir = .GlobalEnv)
  load("files/Phyloseq_objects/EGGNOGs_PhyloseqObj.RData")
  cat("Loading: EggNOGs: Stratified \n"); assign("dat.EGGNOGs", dat.EGGNOGs, envir = .GlobalEnv)
  load("files/Phyloseq_objects/EGGNOGs.slim_PhyloseqObj.RData")
  cat("Loading: EggNOGs \n"); assign("dat.EGGNOGs.slim", dat.EGGNOGs.slim, envir = .GlobalEnv)
}




