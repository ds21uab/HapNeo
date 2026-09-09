#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================================
## Extended Data Figure 4: Immune cell infiltration (quanTIseq) by CCDC110 haplotype
##   a) ORIEN+SU2C (ICI-treated, pre-treatment RNA-seq)
##   b) TCGA NSCLC (ICI-naive)
## Input:  source_data/Source_Data_Extended_Data_Fig_4.xlsx
## Output: figures/ExtData_Fig4a_immune_ORIEN_SU2C.pdf
##         figures/ExtData_Fig4b_immune_TCGA.pdf
## ============================================================================

library(readxl)
library(ggplot2)
library(dplyr)

input_file <- "source_data/Source Data Extended Data Fig.4.xlsx"
dir.create("figures", showWarnings = FALSE)

## ---- Reusable immune boxplot function ---------------------------------------
plot_immune <- function(sheet_name, title_text, output_file) {
  cat("Generating", output_file, "...\n")

  df <- as.data.frame(read_excel(input_file, sheet = sheet_name))

  ## Simplify haplotype labels for x-axis
  df$haplotype <- factor(df$haplotype,
                          levels = c("Hap1/Hap1", "Hap1/Hap2", "Hap2/Hap2"))

  ## Sample sizes
  n_per_hap <- df %>% distinct(sample, haplotype) %>% count(haplotype)

  ## Kruskal-Wallis per cell type
  kw_pvals <- df %>%
    group_by(cell_type) %>%
    summarise(
      p = kruskal.test(fraction ~ haplotype)$p.value,
      .groups = "drop"
    ) %>%
    mutate(p_label = paste0("p = ", signif(p, 2)))

  ## Merge p-values for facet labels
  df <- df %>%
    left_join(kw_pvals, by = "cell_type") %>%
    mutate(facet_label = paste0(cell_type, "\n", p_label))

  ## Order facets by cell type name
  cell_order <- sort(unique(df$cell_type))
  df$facet_label <- factor(df$facet_label,
    levels = df %>%
      distinct(cell_type, facet_label) %>%
      arrange(match(cell_type, cell_order)) %>%
      pull(facet_label)
  )

  p <- ggplot(df, aes(x = haplotype, y = fraction, fill = haplotype)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.3) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c("Hap1/Hap1" = "#E41A1C",
                                  "Hap1/Hap2" = "#377EB8",
                                  "Hap2/Hap2" = "#4DAF4A")) +
    scale_x_discrete(labels = c("Hap1/\nHap1", "Hap1/\nHap2", "Hap2/\nHap2")) +
    labs(title = title_text,
         x = "Haplotype", y = "Cell Fraction", fill = "Haplotype") +
    theme_classic(base_size = 12) +
    theme(
      plot.title     = element_text(face = "bold.italic", hjust = 0.5, size = 16),
      strip.text     = element_text(size = 10, face = "bold"),
      axis.text.x    = element_text(size = 9),
      legend.position = "bottom",
      legend.text    = element_text(size = 12)
    )

  ggsave(output_file, p, width = 14, height = 8, device = cairo_pdf)
}

plot_immune(
  "Extended Data Fig.4a",
  "Immune cell infiltration from pre-ICIs RNA-seq (CCDC110 haplotype from ORIEN and SU2C)",
  "figures/ExtData_Fig4a_immune_ORIEN_SU2C.pdf"
)

plot_immune(
  "Extended Data Fig.4b",
  "Immune cell infiltration from non-ICI RNA-seq (CCDC110 haplotype from TCGA NSCLC)",
  "figures/ExtData_Fig4b_immune_TCGA.pdf"
)

cat("\nExtended Data Figure 4 complete.\n")
