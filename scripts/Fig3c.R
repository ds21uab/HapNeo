#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================================
## Figure 3c: LD regional plot — Haploblock structure around rs7698680
## Panel A: Zoomed view (R² >= 0.8 SNPs with rsID labels)
## Panel B: Full range (all variants)
## Input:  source_data/Source_Data_Fig_3.xlsx (sheet: fig.3c)
## Output: figures/Fig3c.tiff, figures/Fig3c.pdf
## R version: 4.5.2
## ============================================================================

## ---- Load libraries ---------------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(readxl)

## ---- Setup ------------------------------------------------------------------
input_file <- "source_data/Source Data Fig.3.xlsx"
dir.create("figures", showWarnings = FALSE)

## ---- Read source data -------------------------------------------------------
df_raw <- read_excel(input_file, sheet = "fig.3c")

## ---- Define lead SNP --------------------------------------------------------
LEAD_RS <- "rs7698680"
LEAD_BP <- 185459692

## ---- Parse coordinates and assign LD tiers ----------------------------------
df <- df_raw %>%
  mutate(
    RS_Number = trimws(RS_Number),
    pos_bp    = as.numeric(sub(".*:", "", Coord)),
    R2        = as.numeric(R2),
    ld_tier   = case_when(
      R2 >= 0.8 ~ "0.8-1.0",
      R2 >= 0.6 ~ "0.6-0.8",
      R2 >= 0.4 ~ "0.4-0.6",
      R2 >= 0.2 ~ "0.2-0.4",
      TRUE      ~ "< 0.2"
    ),
    ld_tier = factor(ld_tier, levels = c(
      "0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "< 0.2"
    ))
  )

## ---- LD color palette (LocusZoom convention) --------------------------------
ld_colors <- c(
  "0.8-1.0" = "#E24B4A",
  "0.6-0.8" = "#EF9F27",
  "0.4-0.6" = "#E8B020",
  "0.2-0.4" = "#378ADD",
  "< 0.2"   = "#B5D4F4"
)

## ---- Panel A: Zoomed view (R² >= 0.8 only) ----------------------------------
df_zoom <- df %>%
  filter(R2 >= 0.8) %>%
  mutate(label = RS_Number)

xA_range <- range(df_zoom$pos_bp, na.rm = TRUE)
xA_pad   <- 500

pA <- ggplot(df_zoom, aes(x = pos_bp, y = R2)) +
  geom_point(aes(fill = ld_tier),
             shape = 21, size = 6, color = "white", stroke = 0.3) +
  geom_label_repel(
    aes(label = label),
    size = 9, label.padding = unit(0.15, "lines"), label.size = 0.2,
    box.padding = 0.5, min.segment.length = 0.1, max.overlaps = 20,
    na.rm = TRUE, fill = "white", color = "grey30",
    segment.color = "grey60", segment.size = 0.3
  ) +
  scale_fill_manual(values = ld_colors, name = expression(r^2), drop = FALSE) +
  scale_x_continuous(
    name   = "Chromosomal position (chr4, bp)",
    limits = c(xA_range[1] - xA_pad, xA_range[2] + xA_pad),
    breaks = seq(185458500, 185462500, 500),
    labels = function(x) formatC(x, format = "d", big.mark = ",")
  ) +
  scale_y_continuous(
    name   = expression(r^2),
    limits = c(0.8, 1.05),
    breaks = seq(0.8, 1.0, 0.05),
    expand = c(0, 0)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 10, shape = 21))) +
  ggtitle("Haploblock length in European population from 1000G\nchr4:185,458,745-185,462,261 (3,516 bases)") +
  theme_classic(base_size = 14) +
  theme(
    axis.title         = element_text(size = 20),
    axis.text          = element_text(size = 20, color = "grey20"),
    axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position    = "none",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(size = 20, face = "bold"),
    plot.margin        = margin(10, 10, 10, 10)
  )

## ---- Panel B: Full range (all variants) -------------------------------------
x_range <- range(df$pos_bp, na.rm = TRUE)

pB <- ggplot(df, aes(x = pos_bp, y = R2)) +
  geom_point(aes(fill = ld_tier),
             shape = 21, size = 2.5, color = "white", stroke = 0.2, alpha = 0.85) +
  scale_fill_manual(values = ld_colors, name = expression(r^2), drop = FALSE) +
  scale_x_continuous(
    name   = "Chromosomal position (chr4, bp)",
    limits = c(x_range[1] - 1000, x_range[2] + 1000),
    labels = function(x) formatC(x, format = "d", big.mark = ",")
  ) +
  scale_y_continuous(
    name   = expression(r^2),
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5, shape = 21))) +
  theme_classic(base_size = 14) +
  theme(
    axis.title         = element_text(size = 14),
    axis.text          = element_text(size = 12, color = "grey20"),
    legend.title       = element_text(size = 12, face = "bold"),
    legend.text        = element_text(size = 11),
    legend.key.size    = unit(0.5, "cm"),
    legend.position    = "right",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(size = 16, face = "bold"),
    plot.margin        = margin(10, 10, 10, 10)
  )

## ---- Combine panels ---------------------------------------------------------
combined <- pA / pB + plot_layout(guides = "collect") &
  theme(legend.position = "right")

## ---- Save -------------------------------------------------------------------
ggsave("figures/Fig3c.tiff", plot = combined,
       width = 10, height = 8, units = "in",
       dpi = 800, compression = "lzw")

ggsave("figures/Fig3c.pdf", plot = combined,
       width = 10, height = 8, units = "in")

cat("Fig 3c saved to figures/\n")