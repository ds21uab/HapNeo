#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================================
## Extended Data Figure 3: CCDC110 haplotype does not affect mRNA or protein
## expression but is associated with overall survival across NSCLC
##   a) CCDC110 mRNA expression: TCGA tumor vs normal (LUAD, LUSC)
##   b) CCDC110 protein expression: CPTAC (LUAD, LUSC)
##   c) CCDC110 mRNA by haplotype (ORIEN+SU2C, violin)
##   d) CCDC110 protein by haplotype (CPTAC, violin)
##   e) Forest plot: LUAD discovery set
##   f) Forest plot: LUSC discovery set
## Input:  source_data/Source_Data_Extended_Data_Fig_3.xlsx
## Output: figures/ExtData_Fig3a-f.pdf
## ============================================================================

library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(forestploter)

input_file <- "source_data/Source Data Extended Data Fig.3.xlsx"
dir.create("figures", showWarnings = FALSE)


## ============================================================================
## PANEL A: CCDC110 mRNA — TCGA tumor vs normal (LUAD and LUSC)
## ============================================================================

cat("Generating ExtData Fig 3a: mRNA tumor vs normal...\n")

f3a <- as.data.frame(read_excel(input_file, sheet = "Extended Data Fig.3a"))
f3a$Group <- factor(f3a$Group, levels = unique(f3a$Group))

## Wilcoxon comparisons
comparisons <- list(
  c("LUAD_Normal", "LUAD_LUAD"),
  c("LUSC_Normal", "LUSC_LUSC")
)

p3a <- ggplot(f3a, aes(x = Group, y = log2TPM, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                     label = "p.format", size = 5) +
  scale_fill_manual(values = c("Normal" = "#4DAF4A", "Tumor" = "#E41A1C")) +
  labs(title = "CCDC110 (TCGA RNA-seq)",
       x = NULL, y = expression("Expression — log"[2]*"(TPM + 1)"), fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 16),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )

ggsave("figures/ExtData_Fig3a_mRNA_tumor_normal.pdf", p3a, width = 7, height = 5)


## ============================================================================
## PANEL B: CCDC110 protein — CPTAC LUAD vs LUSC
## ============================================================================

cat("Generating ExtData Fig 3b: Protein CPTAC...\n")

f3b <- as.data.frame(read_excel(input_file, sheet = "Extended Data Fig.3b"))
f3b$Histology <- factor(f3b$Histology, levels = c("LUAD", "LUSC"))

wt <- wilcox.test(protein_exp ~ Histology, data = f3b)

p3b <- ggplot(f3b, aes(x = Histology, y = protein_exp, fill = Histology)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.4) +
  annotate("text", x = 1.5, y = max(f3b$protein_exp, na.rm = TRUE) + 1,
           label = paste0("Wilcoxon p = ", signif(wt$p.value, 3)),
           size = 5, fontface = "bold") +
  scale_fill_manual(values = c("LUAD" = "#E41A1C", "LUSC" = "#377EB8")) +
  labs(title = "CCDC110 Protein (CPTAC)",
       x = NULL, y = expression("Protein Expression — log"[2]*"(TMT ratio)"), fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 16),
    legend.position = "none"
  )

ggsave("figures/ExtData_Fig3b_protein_CPTAC.pdf", p3b, width = 5, height = 5)


## ============================================================================
## PANEL C: CCDC110 mRNA by haplotype (ORIEN+SU2C, violin)
## ============================================================================

cat("Generating ExtData Fig 3c: mRNA by haplotype...\n")

f3c <- as.data.frame(read_excel(input_file, sheet = "Extended Data Fig.3c"))
f3c$haplotype <- factor(f3c$haplotype,
                         levels = c("Hap1/Hap1", "Hap1/Hap2", "Hap2/Hap2"))

## Kruskal-Wallis
kw <- kruskal.test(CCDC110_expr ~ haplotype, data = f3c)

## Sample sizes for labels
n_per_hap <- f3c %>% count(haplotype)
x_labels <- paste0(n_per_hap$haplotype, "\n(n=", n_per_hap$n, ")")

p3c <- ggplot(f3c, aes(x = haplotype, y = CCDC110_expr, fill = haplotype)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.3, color = "red") +
  annotate("text", x = 2, y = max(f3c$CCDC110_expr, na.rm = TRUE) + 0.3,
           label = paste0("Kruskal-Wallis p = ", signif(kw$p.value, 3)),
           size = 5, fontface = "bold") +
  scale_x_discrete(labels = x_labels) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  labs(title = "CCDC110 (ORIEN + SU2C RNA-seq)",
       x = "Haplotype",
       y = expression("Expression — log"[2]*"(TPM + 1)"), fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 16),
    legend.position = "none"
  )

ggsave("figures/ExtData_Fig3c_mRNA_by_haplotype.pdf", p3c, width = 6, height = 5)


## ============================================================================
## PANEL D: CCDC110 protein by haplotype (CPTAC, violin)
## ============================================================================

cat("Generating ExtData Fig 3d: Protein by haplotype...\n")

f3d <- as.data.frame(read_excel(input_file, sheet = "Extended Data Fig.3d"))
f3d$haplotype <- factor(f3d$haplotype,
                         levels = c("Hap1/Hap1", "Hap1/Hap2", "Hap2/Hap2"))

kw_d <- kruskal.test(protein_exp ~ haplotype, data = f3d)

n_per_hap_d <- f3d %>% count(haplotype)
x_labels_d <- paste0(n_per_hap_d$haplotype, "\n(n=", n_per_hap_d$n, ")")

p3d <- ggplot(f3d, aes(x = haplotype, y = protein_exp, fill = haplotype)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.3, color = "red") +
  annotate("text", x = 2, y = max(f3d$protein_exp, na.rm = TRUE) + 1,
           label = paste0("Kruskal-Wallis p = ", signif(kw_d$p.value, 3)),
           size = 5, fontface = "bold") +
  scale_x_discrete(labels = x_labels_d) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  labs(title = "CCDC110 (CPTAC)",
       x = "Haplotype",
       y = expression("Protein Expression — log"[2]*"(TMT ratio)"), fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 16),
    legend.position = "none"
  )

ggsave("figures/ExtData_Fig3d_protein_by_haplotype.pdf", p3d, width = 6, height = 5)


## ============================================================================
## PANELS E-F: Forest plots — LUAD and LUSC discovery sets
## ============================================================================

plot_forest <- function(sheet_name, title_text, output_file) {
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

plot_forest("Extended Data Fig.3e", "LUAD (Discovery set)",
            "figures/ExtData_Fig3e_Forest_LUAD.pdf")
plot_forest("Extended Data Fig.3f", "LUSC (Discovery set)",
            "figures/ExtData_Fig3f_Forest_LUSC.pdf")

cat("\nExtended Data Figure 3 complete.\n")
