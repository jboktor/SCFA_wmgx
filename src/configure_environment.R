## Configure enviornment

#### Install CRAN Packages
install.packages(
  c(
    "ggplot2",
    "tidyverse",
    "readxl",
    "plyr",
    "dplyr",
    "ggrepel",
    "gridExtra",
    "reshape2",
    "devtools",
    "RColorBrewer",
    "ggfortify",
    "vegan",
    "MASS",
    "compositions",
    "zCompositions",
    "gplots",
    "viridis",
    "lme4",
    "jtools",
    "phangorn",
    "plotly",
    "VennDiagram",
    "viridis",
    "foreach",
    "doParallel" ,
    "parallel",
    "ggbeeswarm",
    "FSA",
    "ggpubr",
    "ggsci",
    "ggridges",
    "future",
    "svglite",
    "cowplot",
    "coin",
    "EnvStats",
    "sjlabelled",
    "sjmisc",
    "sjPlot",
    "nlme",
    "eulerr",
    "ggthemes",
    "ggforce",
    "huge",
    "Matrix",
    "magrittr",
    "randomForest",
    "pROC",
    "plotROC",
    "SGL",
    "rstatix",
    "mlbench",
    "caret",
    "MLeval",
    "lazyeval",
    "ppcor",
    "ggdendro",
    "tidymodels"
  )
)
# 
#### Install Bioconductor Packages
install.packages("BiocManager")
BiocManager::install(c(
  "phyloseq",
  "microbiome",
  "Biobase",
  "Maaslin2"
))

#### Install Other Github available Packages
install.packages("remotes")
remotes::install_github("gmteunisse/Fantaxtic")
remotes::install_github("schuyler-smith/phyloschuyler")




