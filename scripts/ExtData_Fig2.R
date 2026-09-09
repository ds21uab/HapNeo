#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================================
## Extended Data Figure 2: Multivariate Cox regression (Haplotype + TMB + PD-L1)
##   a) Discovery set (n=117, PD-L1 subset)
##   b) POPLAR cohort (n=49)
## Input:  source_data/Source_Data_Extended_Data_Fig_2.xlsx
## Output: figures/ExtData_Fig2a-b.pdf
## ============================================================================

library(readxl)
library(forestploter)

input_file <- "source_data/Source Data Extended Data Fig.2.xlsx"
dir.create("figures", showWarnings = FALSE)

plot_forest <- function(sheet_name, output_file) {
  cat("Generating", output_file, "...\n")

  df <- as.data.frame(read_excel(input_file, sheet = sheet_name))
  df$`Hazard ratio (95% CI)` <- sprintf("%.2f (%.2f, %.2f)", df$HR, df$CI_lower, df$CI_upper)
  df$P_display <- ifelse(df$P < 0.001, "<0.001", sprintf("%.3f", df$P))

  tm <- forest_theme(base_size = 11, ci_pch = 15, ci_col = "red",
                     ci_lty = 1, ci_lwd = 1.5, refline_col = "grey70", refline_lty = 2)

  p <- forest(
    df[, c("Variable", "N", "Hazard ratio (95% CI)", "P_display")],
    est = df$HR, lower = df$CI_lower, upper = df$CI_upper,
    ci_column = 3, ref_line = 1, theme = tm
  )

  cairo_pdf(output_file, width = 10, height = 6)
  print(p)
  dev.off()
}

plot_forest("Extended Data Fig.1a", "figures/ExtData_Fig2a_Discovery_TMB_PDL1.pdf")
plot_forest("Extended Data Fig.1b", "figures/ExtData_Fig2b_POPLAR_TMB_PDL1.pdf")

cat("Extended Data Figure 2 complete.\n")
