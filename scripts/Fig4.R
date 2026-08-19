#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================
## Figure 4: Kaplan-Meier survival analysis by CCDC110 haplotype
## Panel a: Discovery cohort (ORIEN, n=309)
## Panel b: Hold-out external validation (SU2C, n=78)
## Panel c: POPLAR clinical trial (n=63)
## CCDC110 germline haplotype manuscript
## Input:  Source_Data_Fig4.xlsx
## Output: Fig4a_KM_Discovery_set.pdf
##         Fig4b_KM_Holdout_external_set.pdf
##         Fig4c_KM_POPLAR.pdf
## R version: 4.5.2
## Key packages: survival, survminer, ggplot2, ggtext, readxl
## ============================================================

## ---- Load libraries -----------------------------------------------------
library(readxl)
library(survival)
library(survminer)
library(ggplot2)
library(ggtext)

## ---- Helper functions ---------------------------------------------------

## Format p-value to 4 decimal places
fmt_p <- function(p) {
  sprintf("%.4f", p)
}

## Color p-value: red if significant, grey otherwise
p_col <- function(p) {
  ifelse(p <= 0.05, "red", "grey50")
}

## Columns to coerce from character to numeric
cols_to_numeric <- c(
  "PC1", "PC2", "PC3", "PC4", "PC5", "TMB", "PD-L1",
  "chr4_185458745_G_C", "chr4_185459089_A_C",
  "chr4_185459361_G_A", "chr4_185459692_A_T",
  "chr4_185459961_G_T", "chr4_185460011_G_A",
  "chr4_185461052_A_G"
)

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
    ggtheme      = theme_classic(base_size = 20, base_family = "Arial"),
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
      legend.title        = element_blank(),
      legend.text         = ggtext::element_markdown(size = 20),
      legend.position     = legend_pos,
      legend.key.height   = unit(0.5, "cm"),
      legend.spacing.y    = unit(0.02, "cm"),
      legend.background   = element_blank(),
      plot.margin         = margin(10, 15, 10, 10)
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


## =====================================================================
## PANEL A: Discovery cohort (ORIEN)
## =====================================================================

## ---- Read and prepare data ----------------------------------------------
surv_data_disc <- read_excel("Source_Data_Fig4.xlsx", sheet = "surv_data_disc")
surv_data_disc[cols_to_numeric] <- lapply(surv_data_disc[cols_to_numeric], as.numeric)

## ---- KM plot ------------------------------------------------------------
fit_disc <- survfit(Surv(OS_month, OS_status) ~ Haplotype, data = surv_data_disc)
summary(fit_disc)$table[, "records"]

sd_disc <- survdiff(Surv(OS_month, OS_status) ~ Haplotype, data = surv_data_disc)
p_disc  <- 1 - pchisq(sd_disc$chisq, df = length(sd_disc$n) - 1)

legend_disc <- c(
  "<span style='color:red'>Hap1/Hap1 (n=70)</span>",
  "<span style='color:royalblue1'>Hap1/Hap2 (n=150)</span>",
  "<span style='color:springgreen4'>Hap2/Hap2 (n=89)</span>"
)

annot_disc <- sprintf(
  "Log-rank p = <span style='color:%s'>%s</span>",
  p_col(p_disc), fmt_p(p_disc)
)

cairo_pdf("Fig4a_KM_Discovery_set.pdf", width = 7, height = 5.5)
p_disc_plot <- plot_km(fit_disc, surv_data_disc, legend_disc, annot_disc,
                       annot_x = 20, annot_y = 0.95)
print(p_disc_plot, newpage = FALSE)
dev.off()

## ---- Pairwise log-rank -------------------------------------------------
cat("\n=== Discovery cohort: Pairwise log-rank ===\n")
print(run_pairwise(surv_data_disc))

## ---- Cox regression: covariates -----------------------------------------
cat("\n=== Discovery: Cox ~ Haplotype + Age + Sex + Batch + PC1-5 ===\n")
fit_cox1 <- coxph(
  Surv(OS_month, OS_status) ~ Haplotype + Age + Sex + Batch +
    PC1 + PC2 + PC3 + PC4 + PC5,
  data  = surv_data_disc,
  model = TRUE
)
print(summary(fit_cox1))



## =====================================================================
## PANEL B: Hold-out external validation (SU2C)
## =====================================================================

## ---- Read and prepare data ----------------------------------------------
surv_data_holdout <- read_excel("Source_Data_Fig4.xlsx", sheet = "surv_data_holdout")
surv_data_holdout[cols_to_numeric] <- lapply(surv_data_holdout[cols_to_numeric], as.numeric)

## ---- KM plot ------------------------------------------------------------
fit_hold <- survfit(Surv(OS_month, OS_status) ~ Haplotype, data = surv_data_holdout)
summary(fit_hold)$table[, "records"]

sd_hold <- survdiff(Surv(OS_month, OS_status) ~ Haplotype, data = surv_data_holdout)
p_hold  <- 1 - pchisq(sd_hold$chisq, df = length(sd_hold$n) - 1)

legend_hold <- c(
  "<span style='color:red'>Hap1/Hap1 (n=16)</span>",
  "<span style='color:royalblue1'>Hap1/Hap2 (n=40)</span>",
  "<span style='color:springgreen4'>Hap2/Hap2 (n=22)</span>"
)

annot_hold <- sprintf(
  "Log-rank p = <span style='color:%s'>%s</span>",
  p_col(p_hold), fmt_p(p_hold)
)

cairo_pdf("Fig4b_KM_Holdout_external_set.pdf", width = 7, height = 5.5)
p_hold_plot <- plot_km(fit_hold, surv_data_holdout, legend_hold, annot_hold,
                       annot_x = 20, annot_y = 0.95)
print(p_hold_plot, newpage = FALSE)
dev.off()

## ---- Pairwise log-rank -------------------------------------------------
cat("\n=== Hold-out cohort: Pairwise log-rank ===\n")
print(run_pairwise(surv_data_holdout))

## ---- Cox regression: covariates -----------------------------------------
cat("\n=== Hold-out: Cox ~ Haplotype + Age + Sex + Batch + PC1-5 ===\n")
fit_cox2 <- coxph(
  Surv(OS_month, OS_status) ~ Haplotype + Age + Sex + Batch +
    PC1 + PC2 + PC3 + PC4 + PC5,
  data  = surv_data_holdout,
  model = TRUE
)
print(summary(fit_cox2))


## =====================================================================
## PANEL C: POPLAR clinical trial
## =====================================================================

## ---- Read and prepare data ----------------------------------------------
surv_data_poplar <- read_excel("Source_Data_Fig4.xlsx", sheet = "surv_data_poplar")

## Remove string "NA" rows and drop unused factor levels
surv_data_poplar <- surv_data_poplar[surv_data_poplar$Haplotype != "NA", ]
surv_data_poplar$Haplotype <- droplevels(factor(surv_data_poplar$Haplotype))

surv_data_poplar[cols_to_numeric] <- lapply(surv_data_poplar[cols_to_numeric], as.numeric)

## ---- KM plot ------------------------------------------------------------
fit_pop <- survfit(Surv(OS_month, OS_status) ~ Haplotype, data = surv_data_poplar)
summary(fit_pop)$table[, "records"]

sd_pop <- survdiff(Surv(OS_month, OS_status) ~ Haplotype, data = surv_data_poplar)
p_pop  <- 1 - pchisq(sd_pop$chisq, df = length(sd_pop$n) - 1)

legend_pop <- c(
  "<span style='color:red'>Hap1/Hap1 (n=11)</span>",
  "<span style='color:royalblue1'>Hap1/Hap2 (n=29)</span>",
  "<span style='color:springgreen4'>Hap2/Hap2 (n=23)</span>"
)

annot_pop <- sprintf(
  "Log-rank p = <span style='color:%s'>%s</span>",
  p_col(p_pop), fmt_p(p_pop)
)

cairo_pdf("Fig4c_KM_POPLAR.pdf", width = 7, height = 5.5)
p_pop_plot <- plot_km(fit_pop, surv_data_poplar, legend_pop, annot_pop,
                      annot_x = 5, annot_y = 0.75)
print(p_pop_plot, newpage = FALSE)
dev.off()

## ---- Pairwise log-rank -------------------------------------------------
cat("\n=== POPLAR cohort: Pairwise log-rank ===\n")
print(run_pairwise(surv_data_poplar))

## ---- Cox regression: covariates -----------------------------------------
cat("\n=== POPLAR: Cox ~ Haplotype + Age + Sex + PC1-5 ===\n")
fit_cox3 <- coxph(
  Surv(OS_month, OS_status) ~ Haplotype + Age + Sex +
    PC1 + PC2 + PC3 + PC4 + PC5,
  data  = surv_data_poplar,
  model = TRUE
)
print(summary(fit_cox3))


## ---- Export Source Data for Fig 4 ---------------------------------------
library(writexl)

## Helper: extract forest plot table from coxph object
extract_forest <- function(fit) {
  s <- summary(fit)
  data.frame(
    Variable = rownames(s$coefficients),
    N        = fit$n,
    HR       = round(s$conf.int[, "exp(coef)"], 2),
    CI_lower = round(s$conf.int[, "lower .95"], 2),
    CI_upper = round(s$conf.int[, "upper .95"], 2),
    P        = signif(s$coefficients[, "Pr(>|z|)"], 3),
    row.names = NULL
  )
}

write_xlsx(
  list(
    "fig.4a" = surv_data_disc[, c("OS_month", "OS_status", "Haplotype")],
    "fig.4b" = surv_data_holdout[, c("OS_month", "OS_status", "Haplotype")],
    "fig.4c" = surv_data_poplar[, c("OS_month", "OS_status", "Haplotype")],
    "fig.4d" = extract_forest(fit_cox1),
    "fig.4e" = extract_forest(fit_cox2),
    "fig.4f" = extract_forest(fit_cox3)
  ),
  "Raw_Data_Fig4.xlsx"
)
