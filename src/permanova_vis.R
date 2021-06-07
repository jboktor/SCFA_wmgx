# PERMANOVA PLOTS

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")

permdf <- 
  read.csv(file = "files/permanova_analysis_2021-04-12.csv", header= TRUE) 

#-------------------------------------------------------------------------------
#####                           Wrangling                                 ##### 
#-------------------------------------------------------------------------------

permdf_aitch <- permdf %>% filter(metacat != "AA_notmatched") %>% 
  dplyr::filter(distance == "Aitchisons") %>% 
  mutate(FDR.symbol = sig.symbol.generator(FDR, shh = T))
permdf_bray <- permdf %>% filter(metacat != "AA_notmatched") %>% 
  dplyr::filter(distance == "BrayCurtis") %>% 
  mutate(FDR.symbol = sig.symbol.generator(FDR, shh = T))

#-------------------------------------------------------------------------------
#####                           Plotting                                 ##### 
#-------------------------------------------------------------------------------

perm.aitchisons <- 
  permdf_aitch %>% 
  dplyr::mutate(var_labs = paste0(vars, " (n=", n_meta, ")") ) %>% 
  dplyr::filter(vars %ni% c("donor_id", "description")) %>% 
  ggplot(aes(x=-log10(FDR), y=reorder(var_labs, -log10(FDR)), fill=metacat)) +
  geom_point(aes(size = R2), shape=21, stroke = 0.2) +
  theme_bw() +
  geom_vline(xintercept = -log10(.05), linetype = 1, color = "grey") +
  labs(x = expression('-log' [10]*'[Adjusted P-value]' ), y = NULL, fill = "Metadata", 
       size = expression(R^"2"~"(%)") ) +
  scale_y_discrete(position = "right") +
  scale_fill_nejm() +
  scale_size(breaks = c(5, 10, 20, 60), labels =  c(5, 10, 20, 60)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position =  "left",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(perm.aitchisons, filename = "data/PERMANOVA/PERMANOVA_Aitchisons_summary.svg",
       width = 7, height = 4.5)



perm.Bray <- 
  permdf_bray %>% 
  dplyr::mutate(var_labs = paste0(vars, " (n=", n_meta, ")") ) %>% 
  dplyr::filter(vars %ni% c("donor_id", "description")) %>% 
  ggplot(aes(x=-log10(FDR), y=reorder(var_labs, -log10(FDR)), fill=metacat)) +
  geom_point(aes(size = R2), shape=21, stroke = 0.2) +
  theme_bw() +
  geom_vline(xintercept = -log10(.05), linetype = 1, color = "grey") +
  labs(x = expression('-log' [10]*'[Adjusted P-value]' ), y = NULL, fill = "Metadata", 
       size = expression(R^"2"~"(%)") ) +
  scale_y_discrete(position = "right") +
  scale_fill_nejm() +
  scale_size(breaks = c(5, 10, 20, 60), labels =  c(5, 10, 20, 60)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position =  "left",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(perm.Bray, filename = "data/PERMANOVA/PERMANOVA_BrayCurtis_summary.svg",
       width = 7, height = 4.5)


