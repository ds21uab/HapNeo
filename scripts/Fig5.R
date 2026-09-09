#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================================
## Figure 5: rs7698680 association with OS and predictive accuracy
## Panel a: Forest plot — HR across ICI-treated and ICI-naive cohorts
## Panel b: Time-dependent ROC at 1 year (8cov+TMB vs 8cov+TMB+SNP)
## Panel c: Time-dependent ROC at 2 years (8cov+TMB vs 8cov+TMB+SNP)
##
## Input:  source_data/Source_Data_Fig_5.xlsx
## Output: figures/Fig5a_Forest_multicohort.pdf
##         figures/Fig5b_ROC_1yr.pdf
##         figures/Fig5c_ROC_2yr.pdf
##
## NOTE: ROC curves were generated using the timeROC package from
## controlled-access clinical data (ORIEN+SU2C discovery set, n=309).
## The source data contains pre-computed FPR/TPR curve coordinates.
## To regenerate from raw data, see Methods and the timeROC analysis
## described in the manuscript.
## ============================================================================

## ---- Load libraries ---------------------------------------------------------
library(readxl)
library(ggplot2)
library(ggtext)
library(forestploter)

## ---- Setup ------------------------------------------------------------------
input_file <- "source_data/Source Data Fig.5.xlsx"
dir.create("figures", showWarnings = FALSE)

## ---- Helper functions -------------------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  paste0("P = ", signif(p, 3))
}

p_color <- function(p) {
  if (is.na(p)) return("azure4")
  if (p < 0.05) return("red")
  return("azure4")
}

roc_theme <- theme_classic(base_size = 20) +
  theme(
    axis.text         = element_text(color = "black", size = 20),
    axis.title        = element_text(size = 20, face = "bold"),
    legend.title      = element_blank(),
    legend.position   = c(0.60, 0.12),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(size = 15),
    legend.key.width  = unit(1.5, "cm"),
    plot.margin       = margin(10, 10, 10, 10)
  )


## ============================================================================
## PANEL B: ROC curve at 1 year
## ============================================================================

cat("Generating Fig 5b: ROC 1-year...\n")

df_5b <- as.data.frame(read_excel(input_file, sheet = "Fig.5b"))
df_5b$Model <- factor(df_5b$Model, levels = unique(df_5b$Model))

## Extract AUCs and P-value from model labels
models_5b <- levels(df_5b$Model)
auc_snp_5b <- as.numeric(sub(".*AUC=([0-9.]+).*", "\\1", models_5b[1]))
auc_cov_5b <- as.numeric(sub(".*AUC=([0-9.]+).*", "\\1", models_5b[2]))

## P-value from the manuscript (timeROC::compare, not recomputable from curve coordinates)
pval_5b <- 0.0778

p5b <- ggplot(df_5b, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 1.5) +
  scale_color_manual(values = c("#B2182B", "#FF7F00")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  annotate("text", x = 0.05, y = 0.80, label = fmt_p(pval_5b),
           hjust = 0, size = 7, fontface = "bold.italic",
           color = p_color(pval_5b)) +
  labs(x = "False positive rate", y = "True positive rate",
       title = expression(italic("Prediction at 1 year"))) +
  roc_theme + coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold.italic"))

ggsave("figures/Fig5b_ROC_1yr.pdf", p5b,
       width = 7, height = 6.5, device = cairo_pdf)
cat("  Saved figures/Fig5b_ROC_1yr.pdf\n")


## ============================================================================
## PANEL C: ROC curve at 2 years
## ============================================================================

cat("Generating Fig 5c: ROC 2-year...\n")

df_5c <- as.data.frame(read_excel(input_file, sheet = "Fig.5c"))
df_5c$Model <- factor(df_5c$Model, levels = unique(df_5c$Model))

models_5c <- levels(df_5c$Model)
auc_snp_5c <- as.numeric(sub(".*AUC=([0-9.]+).*", "\\1", models_5c[1]))
auc_cov_5c <- as.numeric(sub(".*AUC=([0-9.]+).*", "\\1", models_5c[2]))

## P-value from the manuscript
pval_5c <- 0.0246

p5c <- ggplot(df_5c, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 1.5) +
  scale_color_manual(values = c("#B2182B", "#FF7F00")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  annotate("text", x = 0.05, y = 0.80, label = fmt_p(pval_5c),
           hjust = 0, size = 7, fontface = "bold.italic",
           color = p_color(pval_5c)) +
  labs(x = "False positive rate", y = "True positive rate",
       title = expression(italic("Prediction at 2 year"))) +
  roc_theme + coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold.italic"))

ggsave("figures/Fig5c_ROC_2yr.pdf", p5c,
       width = 7, height = 6.5, device = cairo_pdf)
cat("  Saved figures/Fig5c_ROC_2yr.pdf\n")

cat("\nFigure 5 complete. Panels saved to figures/\n")