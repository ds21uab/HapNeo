#!/usr/bin/Rscript
## Copyright Divya Sahu, 2026
## ============================================================
## Figure 2d: Bootstrap lollipop plot + CCDC110 gene structure
## CCDC110 germline haplotype manuscript
## Input:  Source_Data_Fig2d.xlsx (Source Data Fig. 2d)
## Output: Fig2d.tiff, Fig2d.pdf
## R version: 4.5.2
## Key packages: ggplot2, patchwork, ggrepel, scales, readxl, dplyr, stringr
## ============================================================

## ---- Load libraries -----------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
  library(readxl)
})

## ---- Read source data ---------------------------------------------------
stab     <- read_excel("Source_Data_Fig2d.xlsx", sheet = "Fig 2d Lollipop")
exons_df <- read_excel("Source_Data_Fig2d.xlsx", sheet = "Fig 2d Gene structure")

## Ensure SNP factor order matches genomic position
stab <- stab %>%
  arrange(pos) %>%
  mutate(SNP = factor(SNP, levels = SNP))

## ---- Shared x-axis limits (genomic window around CCDC110 locus) ---------
x_min <- 185457700
x_max <- 185461400

## ---- Filter exons visible within plotting window ------------------------
exons_vis <- exons_df %>%
  filter(exon_end >= x_min - 500, exon_start <= x_max + 500)

## ---- Label nudges for non-overlapping rsID placement --------------------
label_df <- stab %>%
  mutate(
    nudge_x = case_when(
      rsid == "rs35596415" ~ -420,
      rsid == "rs59319722" ~    0,
      rsid == "rs11132306" ~ -380,
      rsid == "rs7698680"  ~  360,
      rsid == "rs7699687"  ~ -380,
      rsid == "rs7699724"  ~  420,
      rsid == "rs11132309" ~ -420,
      TRUE ~ 0
    ),
    nudge_y = case_when(
      rsid == "rs35596415" ~  10,
      rsid == "rs59319722" ~  12,
      rsid == "rs11132306" ~ -10,
      rsid == "rs7698680"  ~  12,
      rsid == "rs7699687"  ~   9,
      rsid == "rs7699724"  ~  14,
      rsid == "rs11132309" ~   9,
      TRUE ~ 8
    )
  )

## ---- Panel A: Lollipop plot ---------------------------------------------
## Y-axis: training set selection frequency (out of 100 bootstrap splits)
## Point size: test set selection frequency
gg_lollipop <- ggplot(stab, aes(x = pos, y = freq_train)) +

  ## Lollipop stems
  geom_segment(aes(xend = pos, yend = 0),
               colour = "grey60", linewidth = 1.5, alpha = 1) +

  ## Lollipop heads (size encodes test-set frequency)
  geom_point(aes(size = freq_test),
             colour = "#0072B2", alpha = 1) +

  ## rsID labels with leader lines (ggrepel avoids overlap)
  geom_text_repel(
    data = label_df,
    aes(x = pos, y = freq_train, label = rsid),
    nudge_x = label_df$nudge_x,
    nudge_y = label_df$nudge_y,
    size = 10,
    family = "Helvetica",
    fontface = "plain",
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.5,
    max.overlaps = Inf
  ) +

  scale_size_continuous(
    name  = "Test-set selection\nfrequency (out of 100)",
    range = c(3, 7)
  ) +
  scale_y_continuous(
    "Selection frequency in train sets\n(out of 100 splits)",
    limits = c(0, 100), breaks = seq(0, 100, 20)
  ) +
  coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(
    axis.title.x         = element_blank(),
    axis.text.x          = element_blank(),
    axis.ticks.x         = element_blank(),
    axis.title.y         = element_text(size = 22, face = "bold"),
    axis.text.y          = element_text(size = 22, colour = "black"),
    legend.position      = "right",
    legend.justification = c(0, 0.9),
    legend.box.margin    = margin(0, 0, 0, 2),
    legend.margin        = margin(5, 5, 5, 5),
    legend.background    = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
    legend.title         = element_text(size = 20),
    legend.text          = element_text(size = 20),
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    panel.grid.major.y   = element_line(colour = "#AAAAAA", linewidth = 0.5, linetype = "dashed"),
    panel.border         = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.margin          = margin(5, 40, 0, 5)
  )

## ---- Panel B: CCDC110 gene structure ------------------------------------
## Intron line spanning the locus
gene_line <- data.frame(x = x_min - 300, xend = x_max + 300, y = 0.5)

## Intron segments (gaps between exons) for directional arrows
intron_segs <- data.frame(
  start = c(x_min,     185460238, 185461159),
  end   = c(185458125, 185461048, x_max)
)

## Place leftward arrows along introns (CCDC110 is on minus strand)
arrow_pos <- do.call(rbind, lapply(1:nrow(intron_segs), function(i) {
  s <- intron_segs$start[i] + 150
  e <- intron_segs$end[i] - 150
  if (e - s < 200) return(NULL)
  n <- max(2, round((e - s) / 400))
  data.frame(x = seq(s, e, length.out = n))
}))

gg_gene <- ggplot() +

  ## Intron line
  geom_segment(data = gene_line,
               aes(x = x, xend = xend, y = y, yend = y),
               colour = "black", linewidth = 1.2) +

  ## Leftward arrows indicating minus-strand transcription direction
  geom_segment(data = arrow_pos,
               aes(x = x + 120, xend = x - 120, y = 0.5, yend = 0.5),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               colour = "black", linewidth = 1.0) +

  ## Exon boxes (filled black)
  geom_rect(data = exons_vis,
            aes(xmin = pmax(exon_start, x_min - 400),
                xmax = pmin(exon_end, x_max + 400),
                ymin = 0.15, ymax = 0.85),
            fill = "black", colour = "black", linewidth = 0.6) +

  ## Exon number labels below boxes
  geom_text(data = exons_vis,
            aes(x = (pmax(exon_start, x_min) + pmin(exon_end, x_max)) / 2,
                y = -0.1, label = paste0("Exon ", exon_num)),
            size = 7, fontface = "bold", colour = "black",
            family = "Helvetica") +

  ## 5-prime label (right side, minus strand)
  annotate("text", x = x_max + 200, y = 1.0, label = "5 prime",
           size = 8, fontface = "bold", hjust = 0, family = "Helvetica") +

  ## 3-prime label (left side, minus strand)
  annotate("text", x = x_min - 200, y = 1.0, label = "3 prime",
           size = 8, fontface = "bold", hjust = 1, family = "Helvetica") +

  ## Gene name centered above
  annotate("text",
           x = (x_min + x_max) / 2, y = 1.05,
           label = "CCDC110", fontface = "bold.italic",
           size = 8, family = "Helvetica") +

  coord_cartesian(xlim = c(x_min, x_max), ylim = c(-0.5, 1.5), clip = "off") +
  scale_x_continuous("Chromosome 4 position (bp)", labels = comma) +
  theme_bw(base_size = 22, base_family = "Helvetica") +
  theme(
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.text.x  = element_text(size = 22, colour = "black"),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    axis.line.x  = element_line(colour = "black", linewidth = 0.4),
    plot.margin  = margin(0, 5, 5, 5)
  )

## ---- Combine panels with patchwork --------------------------------------
combined <- gg_lollipop / gg_gene +
  plot_layout(heights = c(4, 1.5))

## Display plot
combined

## ---- Save high-resolution outputs ---------------------------------------
ggsave("Fig2d.tiff", plot = combined,
       width = 18, height = 7.5, units = "in",
       dpi = 800, compression = "lzw")

ggsave("Fig2d.pdf", plot = combined,
       width = 14, height = 7.5, units = "in")