


### Read-in MaAsLin2 output
df.mas <-
  read_tsv(
    paste0(
      "data/MaAsLin2_Analysis/Species_genotype_diet_maaslin2_output/all_results.tsv"),
    col_names = T
  )
df.mas.sig <- df.mas %>% filter(qval < 0.25)


df.plot.nonsig <- df.mas %>% 
  filter(qval > 0.25) %>% 
  pivot_wider(
    !value, 
    values_from = c("coef", "stderr", "pval", "qval"), 
    names_from = c("metadata"))
df.plot.sig <- df.mas.sig %>% 
  pivot_wider(
    !value, 
    values_from = c("coef", "stderr", "pval", "qval"), 
    names_from = c("metadata"))



df.plot %>% 
  ggplot() +
  geom_point(data = obj.plot.sig,
             aes(x = coef_diet, y = coef_gen), 
             alpha = 0.2, color = "grey") +
  geom_point(data = obj.plot.sig,
             aes(x = hc_mean_l2fc, y = pc_mean_l2fc, color = group), 
             alpha = 0.7) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "darkgrey") +
  geom_abline(intercept = 0, slope = -1, linetype = 3, color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 1, color = "grey") +
  geom_hline(yintercept = 0, linetype = 1, color = "grey") +
  labs(x = expression(paste("log"[10] * '(PD/HC)')),
       y = expression(paste("log"[10] * '(PD/PC)')),
       title = obj.name,
       color = "Feature \nAssociation") +
  scale_color_manual(values = pal.color) +
  geom_text_repel(data = obj.plot.label,
                  aes(x = hc_mean_l2fc, y = pc_mean_l2fc, label = fullnames), 
                  segment.alpha = 0.5,
                  segment.size = 0.2, 
                  size = 1.75,
                  force = 1,
                  max.time = 1,
                  max.iter = Inf,
                  nudge_y = ifelse(obj.plot.label$pc_mean_l2fc > 0, repelyup, 
                                   ifelse(obj.plot.label$pc_mean_l2fc < 0, -repelydn, 0)),
                  max.overlaps = Inf) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))