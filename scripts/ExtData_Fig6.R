#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================================
## Extended Data Figure 6: MHC-I neoantigen presentation of additional S409F
## predicted peptide (10-mer QGSIIFENEK vs QGSIISENEK)
##   a) ORIEN+SU2C cohort
##   b) POPLAR cohort
## Input:  source_data/Source_Data_Extended_Data_Fig_6.xlsx
## Output: figures/ExtData_Fig6a_10mer_ORIEN_SU2C.pdf
##         figures/ExtData_Fig6b_10mer_POPLAR.pdf
## ============================================================================

library(readxl)
library(ggplot2)

input_file <- "source_data/Source Data Extended Data Fig.6.xlsx"
dir.create("figures", showWarnings = FALSE)

make_mhc_barplot <- function(sheet_name, title_text, output_file) {
  cat("Generating", output_file, "...\n")

  df <- as.data.frame(read_excel(input_file, sheet = sheet_name))

  allele_order <- df %>%
    dplyr::distinct(Allele, diff) %>%
    dplyr::arrange(diff) %>%
    dplyr::pull(Allele)

  df$Allele <- factor(df$Allele, levels = allele_order)
  df$Type   <- factor(df$Type, levels = c("QGSIIFENEK (F409)", "QGSIISENEK (S409)"))

  p <- ggplot(df, aes(x = Rank_EL, y = Allele, fill = Type)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "darkgreen",
               linewidth = 0.8) +
    geom_vline(xintercept = 2.0, linetype = "dashed", color = "orange",
               linewidth = 0.8) +
    annotate("text", x = 0.25, y = Inf, label = "SB", vjust = 1.5,
             color = "darkgreen", size = 5, fontface = "bold") +
    annotate("text", x = 1.25, y = Inf, label = "WB", vjust = 1.5,
             color = "orange", size = 5, fontface = "bold") +
    annotate("text", x = 2.5, y = Inf, label = "NB", vjust = 1.5,
             color = "gray40", size = 5, fontface = "bold") +
    scale_fill_manual(
      values = c("QGSIIFENEK (F409)" = "#CC3333", "QGSIISENEK (S409)" = "#4477AA")
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = title_text,
      x = "pctl. EL Rank (lower = stronger MHC-I binding)",
      y = NULL, fill = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title    = element_text(face = "bold.italic", size = 16, hjust = 0.5),
      legend.position = "bottom",
      legend.text   = element_text(size = 12),
      axis.text.y   = element_text(size = 14, face = "bold"),
      axis.text.x   = element_text(size = 12),
      axis.title.x  = element_text(size = 14)
    )

  ggsave(output_file, p, width = 8, height = 4, device = cairo_pdf)
}

make_mhc_barplot(
  "Extended Data Fig.6a",
  "MHC-I binding of S409F peptides (ORIEN and SU2C)",
  "figures/ExtData_Fig6a_10mer_ORIEN_SU2C.pdf"
)

make_mhc_barplot(
  "Extended Data Fig.6b",
  "MHC-I binding of S409F peptides (POPLAR)",
  "figures/ExtData_Fig6b_10mer_POPLAR.pdf"
)

cat("\nExtended Data Figure 6 complete.\n")