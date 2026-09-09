#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================================
## Figure 4: Kaplan-Meier survival + Forest plots by CCDC110 haplotype
## Panel a: KM — Discovery set (n=309)
## Panel b: KM — Hold-out validation set (n=78)
## Panel c: KM — POPLAR cohort (n=63)
## Panel d-f: Forest plots from Cox regression (read from source data)
## Input:  source_data/Source_Data_Fig_4.xlsx
## Output: figures/Fig4a_KM_Discovery.pdf
##         figures/Fig4b_KM_Holdout.pdf
##         figures/Fig4c_KM_POPLAR.pdf
##         figures/Fig4d_Forest_Discovery.pdf
##         figures/Fig4e_Forest_Holdout.pdf
##         figures/Fig4f_Forest_POPLAR.pdf
## ============================================================================

## ---- Load libraries ---------------------------------------------------------
library(readxl)
library(survival)
library(survminer)
library(ggplot2)
library(ggtext)
library(forestploter)

## ---- Setup ------------------------------------------------------------------
input_file <- "source_data/Source Data Fig.4.xlsx"
dir.create("figures", showWarnings = FALSE)

## ---- Helper functions -------------------------------------------------------

## Format p-value to 4 decimal places
fmt_p <- function(p) sprintf("%.4f", p)

## Color p-value: red if significant, grey otherwise
p_col <- function(p) ifelse(p <= 0.05, "red", "grey50")

## Reusable KM plot function
plot_km <- function(fit, data, legend_labels, annot_colored,
                    annot_x = 20, annot_y = 0.95,
                    legend_pos = c(0.35, 0.15)) {
  p <- ggsurvplot(
    fit,
    data         = data,
    pval         = FALSE,
    risk.table   = FALSE,
    legend.labs  = legend_labels,
    palette      = c("red", "royalblue1", "springgreen4"),
    xlab         = "Time (months)",
    ylab         = "Overall survival probability",
    ggtheme      = theme_classic(base_size = 20),
    size         = 1.5,
    censor.size  = 1
  )

  p$plot <- p$plot +
    ggtext::geom_richtext(
      data        = data.frame(x = annot_x, y = annot_y, label = annot_colored),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust       = 0,
      vjust       = 1,
      size        = 5.8,
      fill        = NA,
      label.color = NA,
      lineheight  = 1.2
    ) +
    theme(
      axis.title         = element_text(size = 20),
      axis.text          = element_text(size = 20),
      legend.title       = element_blank(),
      legend.text        = ggtext::element_markdown(size = 20),
      legend.position    = legend_pos,
      legend.key.height  = unit(0.5, "cm"),
      legend.spacing.y   = unit(0.02, "cm"),
      legend.background  = element_blank(),
      plot.margin        = margin(10, 15, 10, 10)
    )

  return(p)
}

## Reusable pairwise log-rank comparison
run_pairwise <- function(data) {
  pw <- pairwise_survdiff(
    Surv(OS_month, OS_status) ~ Haplotype,
    data            = data,
    p.adjust.method = "none"
  )
  comparisons <- data.frame(
    Comparison = c("Hap2/Hap2 vs Hap1/Hap1",
                   "Hap2/Hap2 vs Hap1/Hap2",
                   "Hap1/Hap2 vs Hap1/Hap1"),
    Log_Rank_P = c(pw$p.value["Hap2/Hap2", "Hap1/Hap1"],
                   pw$p.value["Hap2/Hap2", "Hap1/Hap2"],
                   pw$p.value["Hap1/Hap2", "Hap1/Hap1"])
  )
  comparisons$BH_Adjusted_P <- p.adjust(comparisons$Log_Rank_P, method = "BH")
  return(comparisons)
}


## ============================================================================
## PANEL A: KM — Discovery set
## ============================================================================

cat("Generating Fig 4a: KM Discovery set...\n")

surv_disc <- as.data.frame(read_excel(input_file, sheet = "fig.4a"))
surv_disc$Haplotype <- factor(surv_disc$Haplotype,
                               levels = c("Hap1/Hap1", "Hap1/Hap2", "Hap2/Hap2"))

fit_disc <- survfit(Surv(OS_month, OS_status) ~ Haplotype, data = surv_disc)
sd_disc  <- survdiff(Surv(OS_month, OS_status) ~ Haplotype, data = surv_disc)
p_disc   <- 1 - pchisq(sd_disc$chisq, df = length(sd_disc$n) - 1)

n_disc <- summary(fit_disc)$table[, "records"]

legend_disc <- c(
  sprintf("<span style='color:red'>Hap1/Hap1 (n=%d)</span>", n_disc[1]),
  sprintf("<span style='color:royalblue1'>Hap1/Hap2 (n=%d)</span>", n_disc[2]),
  sprintf("<span style='color:springgreen4'>Hap2/Hap2 (n=%d)</span>", n_disc[3])
)

annot_disc <- sprintf("Log-rank p = <span style='color:%s'>%s</span>",
                      p_col(p_disc), fmt_p(p_disc))

cairo_pdf("figures/Fig4a_KM_Discovery.pdf", width = 7, height = 5.5)
print(plot_km(fit_disc, surv_disc, legend_disc, annot_disc,
              annot_x = 20, annot_y = 0.95), newpage = FALSE)
dev.off()

cat("  Pairwise log-rank:\n")
print(run_pairwise(surv_disc))


## ============================================================================
## PANEL B: KM — Hold-out validation set
## ============================================================================

cat("\nGenerating Fig 4b: KM Hold-out validation...\n")

surv_hold <- as.data.frame(read_excel(input_file, sheet = "fig.4b"))
surv_hold$Haplotype <- factor(surv_hold$Haplotype,
                               levels = c("Hap1/Hap1", "Hap1/Hap2", "Hap2/Hap2"))

fit_hold <- survfit(Surv(OS_month, OS_status) ~ Haplotype, data = surv_hold)
sd_hold  <- survdiff(Surv(OS_month, OS_status) ~ Haplotype, data = surv_hold)
p_hold   <- 1 - pchisq(sd_hold$chisq, df = length(sd_hold$n) - 1)

n_hold <- summary(fit_hold)$table[, "records"]

legend_hold <- c(
  sprintf("<span style='color:red'>Hap1/Hap1 (n=%d)</span>", n_hold[1]),
  sprintf("<span style='color:royalblue1'>Hap1/Hap2 (n=%d)</span>", n_hold[2]),
  sprintf("<span style='color:springgreen4'>Hap2/Hap2 (n=%d)</span>", n_hold[3])
)

annot_hold <- sprintf("Log-rank p = <span style='color:%s'>%s</span>",
                      p_col(p_hold), fmt_p(p_hold))

cairo_pdf("figures/Fig4b_KM_Holdout.pdf", width = 7, height = 5.5)
print(plot_km(fit_hold, surv_hold, legend_hold, annot_hold,
              annot_x = 20, annot_y = 0.95), newpage = FALSE)
dev.off()

cat("  Pairwise log-rank:\n")
print(run_pairwise(surv_hold))


## ============================================================================
## PANEL C: KM — POPLAR cohort
## ============================================================================

cat("\nGenerating Fig 4c: KM POPLAR...\n")

surv_pop <- as.data.frame(read_excel(input_file, sheet = "fig.4c"))

## Remove rows with "NA" haplotype (missing genotype)
surv_pop <- surv_pop[surv_pop$Haplotype != "NA" & !is.na(surv_pop$Haplotype), ]
surv_pop$Haplotype <- factor(surv_pop$Haplotype,
                              levels = c("Hap1/Hap1", "Hap1/Hap2", "Hap2/Hap2"))

fit_pop <- survfit(Surv(OS_month, OS_status) ~ Haplotype, data = surv_pop)
sd_pop  <- survdiff(Surv(OS_month, OS_status) ~ Haplotype, data = surv_pop)
p_pop   <- 1 - pchisq(sd_pop$chisq, df = length(sd_pop$n) - 1)

n_pop <- summary(fit_pop)$table[, "records"]

legend_pop <- c(
  sprintf("<span style='color:red'>Hap1/Hap1 (n=%d)</span>", n_pop[1]),
  sprintf("<span style='color:royalblue1'>Hap1/Hap2 (n=%d)</span>", n_pop[2]),
  sprintf("<span style='color:springgreen4'>Hap2/Hap2 (n=%d)</span>", n_pop[3])
)

annot_pop <- sprintf("Log-rank p = <span style='color:%s'>%s</span>",
                     p_col(p_pop), fmt_p(p_pop))

cairo_pdf("figures/Fig4c_KM_POPLAR.pdf", width = 7, height = 5.5)
print(plot_km(fit_pop, surv_pop, legend_pop, annot_pop,
              annot_x = 5, annot_y = 0.75), newpage = FALSE)
dev.off()

cat("  Pairwise log-rank:\n")
print(run_pairwise(surv_pop))


## ============================================================================
## PANELS D-F: Forest plots from pre-computed Cox regression results
## ============================================================================

cat("\nGenerating Fig 4d-f: Forest plots...\n")

## Reusable forest plot function
plot_forest <- function(sheet_name, title, output_file) {
  df <- as.data.frame(read_excel(input_file, sheet = sheet_name))

  ## Format for display
  df$`Hazard ratio (95% CI)` <- sprintf("%.2f (%.2f, %.2f)", df$HR, df$CI_lower, df$CI_upper)
  df$P_display <- ifelse(df$P < 0.001, "<0.001", sprintf("%.3f", df$P))

  ## Clean variable names for display
  df$Variable <- gsub("^Haplotype", "  ", df$Variable)
  df$Variable <- gsub("^Sex", "  ", df$Variable)
  df$Variable <- gsub("^Batch", "  ", df$Variable)

  ## Add category headers
  header_rows <- data.frame(
    Variable = c("Haplotype", "  Hap1/Hap1"),
    N = c(NA, df$N[1]),
    HR = c(NA, NA),
    CI_lower = c(NA, NA),
    CI_upper = c(NA, NA),
    P = c(NA, NA),
    `Hazard ratio (95% CI)` = c("", "Reference"),
    P_display = c("", ""),
    check.names = FALSE
  )

  ## Create forest plot
  cairo_pdf(output_file, width = 10, height = 6)

  tm <- forest_theme(
    base_size = 11,
    ci_pch = 15,
    ci_col = "red",
    ci_lty = 1,
    ci_lwd = 1.5,
    refline_col = "grey70",
    refline_lty = 2
  )

  p <- forest(
    df[, c("Variable", "N", "Hazard ratio (95% CI)", "P_display")],
    est   = df$HR,
    lower = df$CI_lower,
    upper = df$CI_upper,
    ci_column = 3,
    ref_line = 1,
    theme = tm
  )
  print(p)
  dev.off()
}

plot_forest("fig.4d", "Discovery set",         "figures/Fig4d_Forest_Discovery.pdf")
plot_forest("fig.4e", "Hold-out validation set", "figures/Fig4e_Forest_Holdout.pdf")
plot_forest("fig.4f", "POPLAR cohort",          "figures/Fig4f_Forest_POPLAR.pdf")

cat("\nFigure 4 complete. All panels saved to figures/\n")