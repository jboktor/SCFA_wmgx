### SCFA Analysis ### 


######## Load Data & functions
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/DAF_Functions.R")
# Will take some time
source("src/omixer-rmpR_setup.R")

# ArcSinSqrt - Transformation
df.gbm[-1] <- asin(sqrt(df.gbm[-1]))

df.gbm.plt <- df.gbm %>% dplyr::select(-module) %>% 
  t() %>% melt()
df.gbm.plt <- df.gbm.plt %>% group_col_from_ids(df.gbm.plt$Var1)

# Adjust feature names to match MaAsLin 
df.gbm.plt$Var2 <- gsub(" ", ".", df.gbm.plt$Var2)
df.gbm.plt$Var2 <- gsub("-", ".", df.gbm.plt$Var2)
df.gbm.plt$Var2 <- gsub("[[:punct:][:blank:]]",".", df.gbm.plt$Var2)



# Prep plotting vars
df.gbm.plt$group <- factor(df.gbm.plt$group, levels = c("PC", "PD","HC"))
cols <- c("PC"= "#ed7d31", "PD" = "#bfbfbf", "HC" = "#5b9bd5")
ylabel <- "Normalized Abundance"

# Pull stats from MaAsLin2
Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/GBM_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Population Control")
Maas.pd.hc <- read_tsv(paste0("data/MaAsLin2_Analysis/GBM_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Household Control")

#####################################################################################
PlotGBMs <- function(df.gbm.plt, df.targets, Maas.pd.pc, Maas.pd.hc) {
  for (i in unique(df.targets$Var2)) {
    
    cat("Plotting Feature: ", i, "\n")
    Maas.pd.pc.Q <- Maas.pd.pc %>% filter(feature == i)
    cat("MaAsLin2 PDvPC GLM Q-value: ", Maas.pd.pc.Q$qval, "\n")
    Maas.pd.hc.Q <- Maas.pd.hc %>% filter(feature == i)
    cat("MaAsLin2 PDvHC GLM Q-value: ", Maas.pd.hc.Q$qval, "\n")
    
    df.gbm.plt1 <- df.gbm.plt %>% 
      filter(Var2 == i)
    
    p <- boxplot_all(df.gbm.plt1, x=df.gbm.plt1$group, y=df.gbm.plt1$value,
                     cols=cols, ylabel=ylabel, title= gsub("\\.", " ",  df.gbm.plt1$Var2))
    p <- p + geom_signif(comparisons = list(c("PD", "HC")), #tip_length = 0.02,
                         annotations = sig_mapper(Maas.pd.hc.Q$qval, porq = "q", symbols = F)) +
      geom_signif(comparisons = list(c("PC", "PD")), #tip_length = 0.02,
                  annotations = sig_mapper(Maas.pd.pc.Q$qval, porq = "q", symbols = F)) +
      theme(legend.position = "none",
            panel.grid.major.x = element_blank(),
            plot.title = element_text(size = 14),
            axis.title.y = element_text(size = 12),
            axis.line.x = element_line(colour="grey"),
            axis.line.y = element_line(colour="grey"))
    
    plot(p)
    ggsave(p, filename = paste0("data/GBM_SCFA_Analysis/Boxplot_",i ,".svg"), height = 4, width =3)
    p <- NULL  
  }
}

#####################################################################################

##### Acetate Plots ##### 
## Explore acetate variables
df.acetate <- df.gbm.plt %>% 
  filter(grepl("Acetate", Var2, ignore.case = T))
unique(df.acetate$Var2)

PlotGBMs(df.gbm.plt = df.gbm.plt,
         df.targets = df.acetate,
         Maas.pd.pc = Maas.pd.pc, 
         Maas.pd.hc = Maas.pd.hc)

##### Propionate 
## Explore Propionate variables
df.propionate <- df.gbm.plt %>% 
  filter(grepl("Propionate", Var2, ignore.case = T))
unique(df.propionate$Var2)

PlotGBMs(df.gbm.plt = df.gbm.plt,
         df.targets = df.propionate,
         Maas.pd.pc = Maas.pd.pc, 
         Maas.pd.hc = Maas.pd.hc)

##### Butyrate
## Explore Butyrate variables
df.butyrate <- df.gbm.plt %>% 
  filter(grepl("Butyrate", Var2, ignore.case = T))
unique(df.butyrate$Var2)

PlotGBMs(df.gbm.plt = df.gbm.plt,
         df.targets = df.butyrate,
         Maas.pd.pc = Maas.pd.pc, 
         Maas.pd.hc = Maas.pd.hc)

##### Valerate
## Explore Butyrate variables
df.valerate <- df.gbm.plt %>% 
  filter(grepl("Valer", Var2, ignore.case = T))
unique(df.valerate$Var2)

PlotGBMs(df.gbm.plt = df.gbm.plt,
         df.targets = df.valerate,
         Maas.pd.pc = Maas.pd.pc, 
         Maas.pd.hc = Maas.pd.hc)



##### Filter all significant plots from GBM Model 

Maas.pd.hc.sig <- Maas.pd.hc %>% 
  dplyr::filter(qval < 0.25) %>% 
  dplyr::select(feature)

df.pdhc.sig <- df.gbm.plt %>% 
  filter(Var2 %in% Maas.pd.hc.sig$feature)
unique(df.pdhc.sig$Var2)

PlotGBMs(df.gbm.plt = df.gbm.plt,
         df.targets = df.pdhc.sig,
         Maas.pd.pc = Maas.pd.pc, 
         Maas.pd.hc = Maas.pd.hc)

