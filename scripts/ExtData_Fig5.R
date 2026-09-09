#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================================
## Extended Data Figure 5: MHC-I strong binder analysis of CCDC110 variants
##   a) ORIEN+SU2C cohort
##   b) POPLAR cohort
## Stacked bar: variant-only SB (red), both SB (yellow), reference-only SB (blue)
## Input:  source_data/Source_Data_Extended_Data_Fig_5.xlsx
## Output: figures/ExtData_Fig5a_SB_ORIEN_SU2C.pdf
##         figures/ExtData_Fig5b_SB_POPLAR.pdf
## ============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

input_file <- "source_data/Source Data Extended Data Fig.5.xlsx"
dir.create("figures", showWarnings = FALSE)

## ---- Reusable function ------------------------------------------------------
plot_sb_barplot <- function(sheet_name, title_text, output_file) {
  cat("Generating", output_file, "...\n")

  raw <- as.data.frame(read_excel(input_file, sheet = sheet_name))

  ## Compute SB counts per variant
  sb_counts <- raw %>%
    filter(Variant_ID %in% c("S409F", "P209Q", "I614M", "L299M", "Y500D")) %>%
    group_by(Variant_ID) %>%
    summarise(
      VAR_only = sum(Rank_EL_VAR < 0.5 & Rank_EL_WT >= 0.5, na.rm = TRUE),
      Both     = sum(Rank_EL_VAR < 0.5 & Rank_EL_WT < 0.5, na.rm = TRUE),
      REF_only = sum(Rank_EL_WT < 0.5 & Rank_EL_VAR >= 0.5, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Variant_ID = factor(Variant_ID,
                                levels = c("S409F", "P209Q", "I614M", "L299M", "Y500D")))

  ## Minimum display height for visibility
  min_height <- 4

  df_display <- sb_counts %>%
    mutate(
      VAR_only_d = ifelse(VAR_only > 0, pmax(VAR_only, min_height), 0),
      Both_d     = ifelse(Both > 0, pmax(Both, min_height), 0),
      REF_only_d = ifelse(REF_only > 0, pmax(REF_only, min_height), 0)
    )

  ## Reshape for stacking
  df_long <- df_display %>%
    select(Variant = Variant_ID,
           VAR_only = VAR_only_d, Both = Both_d, REF_only = REF_only_d) %>%
    pivot_longer(-Variant, names_to = "Category", values_to = "Count") %>%
    mutate(
      Count = ifelse(Category == "REF_only", -Count, Count),
      Category = factor(Category,
        levels = c("VAR_only", "Both", "REF_only"),
        labels = c(
          "Variant peptide\u2013MHC-I allele pair SB, Reference not SB",
          "Both peptide\u2013MHC-I allele pairs SB",
          "Reference peptide\u2013MHC-I allele pair SB, Variant not SB"
        )
      )
    )

  ## Label positions (actual counts as labels)
  label_df <- df_display %>%
    rename(Variant = Variant_ID) %>%
    mutate(
      y_var_only = Both_d + VAR_only_d / 2,
      y_both     = Both_d / 2,
      y_ref_only = -REF_only_d / 2
    )

  p <- ggplot(df_long, aes(x = Variant, y = Count, fill = Category)) +
    geom_col(width = 0.65) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_text(data = label_df %>% filter(VAR_only > 0),
              aes(x = Variant, y = y_var_only, label = VAR_only),
              inherit.aes = FALSE, size = 7, fontface = "bold", color = "white") +
    geom_text(data = label_df %>% filter(Both > 0),
              aes(x = Variant, y = y_both, label = Both),
              inherit.aes = FALSE, size = 7, fontface = "bold") +
    geom_text(data = label_df %>% filter(REF_only > 0),
              aes(x = Variant, y = y_ref_only, label = REF_only),
              inherit.aes = FALSE, size = 7, fontface = "bold", color = "white") +
    scale_fill_manual(values = c(
      "Variant peptide\u2013MHC-I allele pair SB, Reference not SB" = "red1",
      "Both peptide\u2013MHC-I allele pairs SB"                      = "#F0C75E",
      "Reference peptide\u2013MHC-I allele pair SB, Variant not SB"  = "royalblue3"
    )) +
    labs(title = title_text,
         y = "Number of Strong Binders (pctl. EL Rank < 0.5)",
         fill = NULL) +
    theme_classic(base_size = 11) +
    theme(
      legend.position  = "bottom",
      legend.direction = "vertical",
      legend.text      = element_text(size = 18),
      plot.title       = element_text(hjust = 0.5, face = "bold.italic", size = 18),
      axis.text.x      = element_text(face = "bold", size = 18),
      axis.text.y      = element_text(face = "bold", size = 18),
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(face = "bold", size = 18)
    )

  ggsave(output_file, p, width = 9, height = 7, device = cairo_pdf)
}

plot_sb_barplot(
  "Extended Data Fig.5a",
  "MHC-I Strong Binder Count by CCDC110 Variant (ORIEN and SU2C)",
  "figures/ExtData_Fig5a_SB_ORIEN_SU2C.pdf"
)

plot_sb_barplot(
  "Extended Data Fig.5b",
  "MHC-I Strong Binder Count by CCDC110 Variant (POPLAR)",
  "figures/ExtData_Fig5b_SB_POPLAR.pdf"
)

cat("\nExtended Data Figure 5 complete.\n")