## Statistical test functions

#-------------------------------------------------------------------------
#                       Two-way ANOVA + Tukey
#-------------------------------------------------------------------------

anova_tukey <- function(test_col, df){
  
  # # TROUBLE
  # test_col <- "Observed"
  # df = env
  
  formula <- as.formula(paste(test_col, "~ treatment_group"))
  aov_oneway <- aov(formula, data = df)
  print(summary(aov_oneway))
  tukey <- TukeyHSD(aov_oneway, which = "treatment_group")
  return(tukey)
}

#-------------------------------------------------------------------------