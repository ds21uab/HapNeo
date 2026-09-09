#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================================
## Figure 6: MHC-I neoantigen presentation analysis of S409F in CCDC110
## Panel a: NetChop proteasomal cleavage (F409 vs S409)
## Panel b: MHC-I binding bar plot (ORIEN+SU2C)
## Panel c: MHC-I binding bar plot (POPLAR)
## Panel d: PHBR density plot (ORIEN+SU2C)
## Panel e: PHBR density plot (POPLAR)
## Panel f: Presenter/Non-presenter bar plot
## Panel g: KM survival by presentation status
## Input:  source_data/Source_Data_Fig_6.xlsx
## Output: figures/Fig6a-g PDFs
## ============================================================================

## ---- Load libraries ---------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
  library(ggrepel)
  library(survival)
  library(survminer)
})

## ---- Setup ------------------------------------------------------------------
input_file <- "source_data/Source Data Fig.6.xlsx"
dir.create("figures", showWarnings = FALSE)


## ============================================================================
## PANEL A: NetChop proteasomal cleavage barplot (F409 vs S409)
## ============================================================================

cat("Generating Fig 6a: NetChop cleavage...\n")

## Read source data (full protein, 834 rows)
## Columns have duplicate names; use .name_repair to make unique
f6a <- as.data.frame(read_excel(input_file, sheet = "fig.6a", .name_repair = "unique"))

## Extract F409 and S409 data for positions 395-415
positions <- 395:415
idx <- positions  # row indices matching position numbers

df_var <- data.frame(
  pos   = positions,
  aa    = as.character(f6a[idx, 2]),   # amino_acid (F409)
  score = as.numeric(f6a[idx, 3])      # prediction_score (F409)
)

df_wt <- data.frame(
  pos   = positions,
  aa    = as.character(f6a[idx, 5]),   # amino_acid (S409)
  score = as.numeric(f6a[idx, 6])      # prediction_score (S409)
)

df_var$is_mut <- df_var$pos == 409
df_wt$is_mut  <- FALSE

show_pos <- c(395, 409, 415)
df_var$label <- ifelse(df_var$pos %in% show_pos, paste0(df_var$aa, "\n", df_var$pos), df_var$aa)
df_wt$label  <- ifelse(df_wt$pos %in% show_pos, paste0(df_wt$aa, "\n", df_wt$pos), df_wt$aa)

df_var$score_label <- ifelse(df_var$score > 0.45, sprintf("%.2f", df_var$score), "")
df_wt$score_label  <- ifelse(df_wt$score > 0.45, sprintf("%.2f", df_wt$score), "")

bar_high  <- "#D85A30"
bar_low   <- "#B4B2A9"
mut_color <- "#A32D2D"

make_netchop <- function(df, title_text, pair1_end, rect_fill, rect_border) {
  label_colors <- ifelse(df$is_mut, mut_color, "black")
  arrow_pos <- df$pos[df$score > 0.45]
  has_arrows <- length(arrow_pos) > 0
  if (has_arrows) df_arrows <- data.frame(x = arrow_pos, y = -0.05)

  p <- ggplot(df, aes(x = pos, y = score)) +
    annotate("rect", xmin = 400.5, xmax = pair1_end + 0.5,
             ymin = -0.08, ymax = 1.10,
             fill = rect_fill, alpha = 0.15, color = rect_border,
             linewidth = 0.6, linetype = "solid") +
    annotate("rect", xmin = 403.5, xmax = 413.5,
             ymin = -0.08, ymax = 1.10,
             fill = rect_fill, alpha = 0.10, color = rect_border,
             linewidth = 0.6, linetype = "dotted") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_col(aes(fill = score > 0.45), width = 0.7, color = NA) +
    scale_fill_manual(values = c("TRUE" = bar_high, "FALSE" = bar_low), guide = "none") +
    geom_text(aes(label = score_label), vjust = -0.5, size = 2.5,
              fontface = "bold", color = "darkblue") +
    scale_x_continuous(breaks = positions, labels = df$label) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    coord_cartesian(ylim = c(-0.08, 1.18), clip = "off") +
    ggtitle(title_text) +
    labs(x = NULL, y = "Prediction score") +
    theme_classic(base_size = 9) +
    theme(
      plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 6, lineheight = 0.85,
                                 color = label_colors,
                                 face = ifelse(df$is_mut, "bold", "plain")),
      axis.text.y  = element_text(size = 6),
      axis.title.y = element_text(size = 7),
      axis.line    = element_line(linewidth = 0.3),
      plot.margin  = margin(6, 6, 6, 2)
    )

  if (has_arrows) {
    p <- p + geom_point(data = df_arrows, aes(x = x, y = y),
                        shape = 24, size = 2.5, fill = "#2E8B57", color = "#2E8B57")
  }
  return(p)
}

p6a_var <- make_netchop(df_var, "F409", 409, "#FF6B6B", "#CC3333")
p6a_wt  <- make_netchop(df_wt, "S409", 408, "#6BB5FF", "#2266CC")

p6a <- (p6a_var | p6a_wt) +
  plot_annotation(
    caption = "CCDC110 amino acid sequence (positions 395\u2013415)",
    theme = theme(plot.caption = element_text(size = 7, face = "bold", hjust = 0.5))
  )

ggsave("figures/Fig6a_NetChop.pdf", p6a, width = 180/25.4, height = 60/25.4)
cat("  Saved Fig6a\n")


## ============================================================================
## PANELS B-C: MHC-I binding bar plots (ORIEN+SU2C and POPLAR)
## ============================================================================

make_mhc_barplot <- function(sheet_name, title_text, output_file) {
  cat("Generating", output_file, "...\n")

  df <- as.data.frame(read_excel(input_file, sheet = sheet_name))

  ## Order alleles by diff (strongest binding difference first)
  allele_order <- df %>%
    distinct(Allele, diff) %>%
    arrange(diff) %>%
    pull(Allele)

  df$Allele <- factor(df$Allele, levels = allele_order)
  df$Type   <- factor(df$Type, levels = c("LVKQGSIIF (F409)", "LVKQGSII (S409)"))

  ## SB threshold
  sb_thresh <- 0.5

  p <- ggplot(df, aes(x = Rank_EL, y = Allele, fill = Type)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_vline(xintercept = sb_thresh, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    scale_fill_manual(
      values = c("LVKQGSIIF (F409)" = "#CC3333", "LVKQGSII (S409)" = "#4477AA")
    ) +
    labs(
      title = title_text,
      subtitle = paste0("SB threshold (pctl. EL Rank < ", sb_thresh, ")"),
      x = "pctl. EL Rank (lower = stronger MHC-I binding)",
      y = NULL, fill = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title    = element_text(face = "bold.italic", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom",
      legend.text   = element_text(size = 12),
      axis.text.y   = element_text(size = 12, face = "bold"),
      axis.text.x   = element_text(size = 12)
    )

  ggsave(output_file, p, width = 8, height = 5, device = cairo_pdf)
}

make_mhc_barplot("fig.6b", "MHC-I binding of S409F peptides (ORIEN and SU2C)",
                 "figures/Fig6b_MHC_binding_ORIEN_SU2C.pdf")
make_mhc_barplot("fig.6c", "MHC-I binding of S409F peptides (POPLAR)",
                 "figures/Fig6c_MHC_binding_POPLAR.pdf")


## ============================================================================
## PANELS D-E: PHBR density plots (ORIEN+SU2C and POPLAR)
## ============================================================================

make_density <- function(sheet_name, title_text, p_value, output_file) {
  cat("Generating", output_file, "...\n")

  df <- as.data.frame(read_excel(input_file, sheet = sheet_name))
  df$Group <- factor(df$Group, levels = unique(df$Group))

  ## Sample sizes
  n_var <- sum(grepl("variant|Hap2", df$Group, ignore.case = TRUE))
  n_ref <- nrow(df) - n_var

  ## Format p-value
  p_label <- if (p_value < 2.2e-16) {
    "p < 2.2e-16"
  } else {
    paste0("p = ", formatC(p_value, format = "e", digits = 2))
  }

  x_range <- range(df$log10_PHBR, na.rm = TRUE)

  p <- ggplot(df, aes(x = log10_PHBR, fill = Group, color = Group)) +
    geom_density(alpha = 0.45, linewidth = 1) +
    geom_vline(xintercept = log10(0.5), linetype = "dashed",
               color = "red", linewidth = 1.5) +
    geom_vline(xintercept = log10(2.0), linetype = "dashed",
               color = "orange", linewidth = 1.5) +
    annotate("text", x = log10(0.5) - 0.15, y = Inf,
             label = "SB", vjust = 1.5, color = "red", size = 15, fontface = "bold") +
    annotate("text", x = (log10(0.5) + log10(2.0)) / 2, y = Inf,
             label = "WB", vjust = 1.5, color = "orange", size = 15, fontface = "bold") +
    annotate("text", x = log10(2.0) + 0.15, y = Inf,
             label = "NB", vjust = 1.5, color = "gray40", size = 15, fontface = "bold") +
    annotate("text", x = x_range[1], y = 1, label = p_label,
             hjust = 0, vjust = 2, size = 10, fontface = "bold.italic") +
    scale_fill_manual(values = c("#C0392B", "#4477AA")) +
    scale_color_manual(values = c("#C0392B", "#4477AA")) +
    coord_cartesian(clip = "off") +
    labs(
      title    = title_text,
      subtitle = paste0("Variant (Hap2 carriers, n = ", n_var,
                        ") | Reference (Hap1/Hap1, n = ", n_ref, ")"),
      x = "log10(PHBR Score)", y = "Density", fill = "", color = ""
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title    = element_text(size = 20, face = "bold.italic", hjust = 0.5),
      legend.position = "bottom",
      legend.text   = element_text(size = 20),
      axis.text     = element_text(size = 20),
      axis.title.x  = element_text(size = 20),
      axis.title.y  = element_text(size = 20),
      axis.line     = element_line(linewidth = 0.5)
    )

  ggsave(output_file, p, width = 11, height = 7, device = cairo_pdf)
}

## P-values from Fisher's exact test (computed from patient-level data, not recomputable from source)
make_density("fig.6d",
             "S409F Variant vs Reference Allele Peptide (ORIEN and SU2C)",
             p_value = 1.59e-14,
             "figures/Fig6d_PHBR_density_ORIEN_SU2C.pdf")

make_density("fig.6e",
             "S409F Variant vs Reference Allele Peptide (POPLAR)",
             p_value = 1.50e-02,
             "figures/Fig6e_PHBR_density_POPLAR.pdf")


## ============================================================================
## PANEL F: Presenter / Non-presenter barplot
## ============================================================================

cat("Generating Fig 6f: Presenter barplot...\n")

f6f <- as.data.frame(read_excel(input_file, sheet = "fig.6f"))
f6f$Peptide      <- factor(f6f$Peptide, levels = c("F409", "S409"))
f6f$Presentation <- factor(f6f$Presentation, levels = c("Presenter", "Non-presenter"))
f6f$Label <- paste0(f6f$Percent, "%\n(", f6f$N, ")")

## Fisher's exact test P-value (from manuscript)
fisher_p <- 2.411e-13

subtitle_text <- paste0(
  "Threshold = 1.7<br><b><i>Fisher exact p = ",
  format(fisher_p, scientific = TRUE, digits = 4), "</i></b>"
)

p6f <- ggplot(f6f, aes(x = Peptide, y = Percent, fill = Presentation)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black") +
  geom_text(aes(label = Label), position = position_dodge(width = 0.75),
            vjust = -0.25, size = 10) +
  scale_fill_manual(values = c("Presenter" = "#2EAD6B", "Non-presenter" = "#D95F5F")) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 25, 50, 75, 100),
                     labels = c("0%", "25%", "50%", "75%", "100%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "PHBR Score Presentation Status (ORIEN and SU2C)",
       subtitle = subtitle_text, x = "", y = "% of patients", fill = "") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    plot.title    = element_text(face = "bold.italic", size = 24, hjust = 0.5),
    plot.subtitle = element_markdown(size = 23, hjust = 0.5),
    axis.title.y  = element_text(size = 24),
    axis.text     = element_text(size = 24),
    legend.text   = element_text(size = 24)
  )

ggsave("figures/Fig6f_Presenter_barplot.pdf", p6f, width = 10, height = 10, device = cairo_pdf)
cat("  Saved Fig6f\n")


## ============================================================================
## PANEL G: KM survival by presentation status (4 groups)
## ============================================================================

cat("Generating Fig 6g: KM presentation status...\n")

f6g <- as.data.frame(read_excel(input_file, sheet = "fig.6g"))
f6g$Group <- factor(f6g$Group, levels = c(
  "F409 Presenter", "F409 Non-presenter",
  "S409 Presenter", "S409 Non-presenter"
))

## KM fit
km_fit <- survfit(Surv(OS_month, OS_status) ~ Group, data = f6g)

## Log-rank test
logrank <- survdiff(Surv(OS_month, OS_status) ~ Group, data = f6g)
logrank_p <- pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
logrank_label <- paste0("p = ", format(logrank_p, scientific = TRUE, digits = 3))

## Sample sizes for legend
group_n <- table(f6g$Group)
legend_labs <- paste0(names(group_n), " (n = ", group_n, ")")

km_plot <- ggsurvplot(
  km_fit,
  data        = f6g,
  size        = 2,
  risk.table  = FALSE,
  conf.int    = FALSE,
  legend.title = "",
  legend.labs = legend_labs,
  legend      = c(0.68, 0.82),
  palette     = c("#1B9E77", "#D95F02", "#377EB8", "#E41A1C"),
  xlab        = "Time in months",
  ylab        = "Overall survival probability",
  title       = "PHBR Score Presentation Status (ORIEN and SU2C)",
  ggtheme     = theme_classic(base_size = 20)
)

km_plot$plot <- km_plot$plot +
  annotate("text", x = 2, y = 0.16, label = logrank_label,
           size = 12, fontface = "bold.italic", hjust = 0) +
  labs(subtitle = "Threshold = 1.7") +
  theme(
    legend.position    = c(0.68, 0.82),
    legend.direction   = "vertical",
    legend.background  = element_blank(),
    plot.title    = element_text(face = "bold.italic", size = 24, hjust = 0.5),
    plot.subtitle = element_text(size = 22, hjust = 0.5),
    axis.title    = element_text(size = 23),
    axis.text     = element_text(size = 23),
    legend.text   = element_text(size = 23),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(1.0, "cm")
  )

ggsave("figures/Fig6g_KM_presentation.pdf", km_plot$plot,
       width = 10, height = 7, device = cairo_pdf)
cat("  Saved Fig6g\n")

cat("\nFigure 6 complete. All panels saved to figures/\n")