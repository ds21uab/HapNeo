#!/usr/bin/Rscript
## ============================================================================
## Fig2.R — Figure 2: EWAS analysis identifies germline variants in CCDC110
##
## Panels:
##   a) Manhattan plot — Discovery cohort GWAS
##   b) Venn diagram  — Discovery, Hold-out validation, POPLAR
##   c) Venn diagram  — Testing set, Hold-out validation, POPLAR
##   d) Lollipop plot + CCDC110 gene structure
##
## Input:  source_data/Source_Data_Fig_2.xlsx (sheets: fig.2a, fig.2b, fig.2c, fig.2d)
## Output: figures/Fig2a.tiff, Fig2b.pdf, Fig2c.pdf, Fig2d.tiff, Fig2d.pdf
##
## Copyright Divya Sahu, 2026
## ============================================================================

## ---- Load libraries ---------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(ggVennDiagram)
  library(sf)
  library(patchwork)
  library(scales)
  library(grid)
  library(writexl)
})

## ---- Setup ------------------------------------------------------------------
input_file <- "source_data/Source Data Fig.2.xlsx"
dir.create("figures", showWarnings = FALSE)


## ============================================================================
## PANEL A: Manhattan plot — Discovery cohort GWAS
## ============================================================================

cat("Generating Fig 2a: Manhattan plot...\n")

f2a <- read_excel(input_file, sheet = "fig.2a")

## Lead SNP and LD SNPs
lead_snp   <- "chr4_185459692_A_T"
lead_label <- "CCDC110 (rs7698680)"

ccdc110_ld_snps <- c(
  "chr4_185458745_G_C", "chr4_185461052_A_G", "chr4_185459089_A_C",
  "chr4_185459361_G_A", "chr4_185459961_G_T", "chr4_185460011_G_A"
)

## Prepare data
gwas2 <- f2a %>%
  mutate(P = ifelse(is.na(P) | P == 0, .Machine$double.xmin, P)) %>%
  filter(!is.na(P)) %>%
  mutate(
    CHR  = as.integer(str_match(SNP, "^chr([0-9]+)_")[, 2]),
    BP   = as.integer(str_match(SNP, "^chr[0-9]+_([0-9]+)_")[, 2]),
    LOGP = -log10(P)
  ) %>%
  filter(!is.na(CHR), !is.na(BP)) %>%
  arrange(CHR, BP)

## Cumulative chromosome positions
chr_pos <- gwas2 %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(
    chr_len   = as.numeric(chr_len),
    chr_start = dplyr::lag(cumsum(chr_len), default = 0)
  )

gwas2 <- gwas2 %>%
  left_join(chr_pos, by = "CHR") %>%
  mutate(
    BP_cum     = as.numeric(BP) + chr_start,
    chr_parity = CHR %% 2
  )

axis_df <- chr_pos %>% mutate(center = chr_start + chr_len / 2)

## Thresholds
thresh_suggestive <- 5e-4
thresh_gw         <- 5e-8

## Plot
p2a <- ggplot(gwas2, aes(x = BP_cum, y = LOGP)) +
  geom_point(aes(color = factor(chr_parity)), size = 0.45, alpha = 0.85) +
  scale_color_manual(values = c("grey60", "#0072B2"), guide = "none") +
  geom_point(
    data = subset(gwas2, SNP %in% ccdc110_ld_snps),
    aes(shape = "SNPs in LD with rs7698680"),
    color = "black", size = 2.5, stroke = 1
  ) +
  geom_point(
    data = subset(gwas2, SNP == lead_snp),
    aes(shape = "Lead SNP rs7698680"),
    color = "red", size = 2.5
  ) +
  geom_hline(aes(yintercept = -log10(thresh_suggestive), linetype = "thresh_line"),
             color = "deeppink", linewidth = 0.8) +
  geom_hline(yintercept = -log10(thresh_gw), color = "grey40",
             linetype = "dotted", linewidth = 0.8) +
  scale_shape_manual(
    name = NULL,
    values = c("Lead SNP rs7698680" = 16, "SNPs in LD with rs7698680" = 4),
    guide = guide_legend(order = 1,
      override.aes = list(color = c("red", "black"), size = c(2.5, 2.5), stroke = c(0, 1)))
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("thresh_line" = "solid"),
    labels = "*P* < 5 \u00d7 10<sup>\u22124</sup>",
    guide = guide_legend(order = 2, override.aes = list(color = "deeppink"))
  ) +
  geom_text_repel(
    data = subset(gwas2, SNP == lead_snp),
    aes(label = lead_label),
    color = "red", size = 6, fontface = "bold",
    nudge_y = 0.6, min.segment.length = 0
  ) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR,
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(bquote(bold(-log[10](P))), breaks = seq(0, 8, 2),
                     expand = expansion(mult = c(0, 0.05))) +
  xlab("Chromosome") +
  ggtitle("Discovery set") +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold.italic", hjust = 0.5, size = 20),
    axis.title         = element_text(face = "bold"),
    axis.title.x       = element_text(size = 17, face = "bold"),
    axis.title.y       = element_text(size = 17, face = "bold"),
    axis.text.x        = element_text(size = 15, angle = 0, hjust = 0.5),
    axis.text.y        = element_text(size = 15),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border       = element_rect(colour = "black", fill = NA),
    legend.position    = c(0.85, 0.85),
    legend.background  = element_blank(),
    legend.box.background = element_blank(),
    legend.text        = element_markdown(size = 12, face = "bold"),
    legend.key         = element_blank(),
    legend.key.height  = unit(0.4, "cm"),
    legend.key.width   = unit(0.86, "cm"),
    legend.margin      = margin(0, 0, 0, 0)
  )

ggsave("figures/Fig2a.tiff", p2a, width = 12, height = 5,
       units = "in", dpi = 800, compression = "lzw")
cat("  Saved figures/Fig2a.tiff\n")


## ============================================================================
## PANEL B: Venn diagram — Discovery, Hold-out validation, POPLAR
## ============================================================================

cat("Generating Fig 2b: Venn diagram (Discovery approach)...\n")

f2b <- read_excel(input_file, sheet = "fig.2b")

x_b <- list(
  Discovery             = na.omit(f2b$Discovery),
  POPLAR                = na.omit(f2b$POPLAR),
  `Hold-out validation` = na.omit(f2b$`Hold-out validation`)
)

p2b <- ggVennDiagram(x_b, label = "count", label_alpha = 0,
                     label_size = 7, set_size = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_color_manual(values = rep("grey40", 3)) +
  theme(plot.margin = margin(60, 40, 60, 40), legend.position = "none") +
  coord_sf(clip = "off") +
  annotate("text", x = -4.5, y = 6.5, label = "Discovery",
           fontface = "bold.italic", size = 10) +
  annotate("text", x = 8.5, y = 6.5, label = "POPLAR",
           fontface = "bold.italic", size = 10) +
  annotate("text", x = 2.0, y = -10.5, label = "Hold-out validation",
           fontface = "bold.italic", size = 10)

ggsave("figures/Fig2b.pdf", p2b, width = 7, height = 5, dpi = 800)
cat("  Saved figures/Fig2b.pdf\n")


## ============================================================================
## PANEL C: Venn diagram — Testing set, Hold-out validation, POPLAR
## ============================================================================

cat("Generating Fig 2c: Venn diagram (Iterative resampling approach)...\n")

f2c <- read_excel(input_file, sheet = "fig.2c")

x_c <- list(
  `Testing set`         = na.omit(f2c$`Testing set`),
  POPLAR                = na.omit(f2c$POPLAR),
  `Hold-out validation` = na.omit(f2c$`Hold-out validation`)
)

p2c <- ggVennDiagram(x_c, label = "count", label_alpha = 0,
                     label_size = 7, set_size = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_color_manual(values = rep("grey40", 3)) +
  theme(plot.margin = margin(60, 40, 60, 40), legend.position = "none") +
  coord_sf(clip = "off") +
  annotate("text", x = -4.5, y = 6.5, label = "Testing set",
           fontface = "bold.italic", size = 10) +
  annotate("text", x = 8.5, y = 6.5, label = "POPLAR",
           fontface = "bold.italic", size = 10) +
  annotate("text", x = 2.0, y = -10.5, label = "Hold-out validation",
           fontface = "bold.italic", size = 10)

ggsave("figures/Fig2c.pdf", p2c, width = 7, height = 5, dpi = 800)
cat("  Saved figures/Fig2c.pdf\n")


## ============================================================================
## PANEL D: Lollipop plot + CCDC110 gene structure
## ============================================================================

cat("Generating Fig 2d: Lollipop + gene structure...\n")

stab     <- read_excel(input_file, sheet = "fig.2d")
exons_df <- read_excel(input_file, sheet = "fig.2d_gene_structure")

stab <- stab %>% arrange(pos) %>% mutate(SNP = factor(SNP, levels = SNP))

## Genomic window
x_min <- 185457700
x_max <- 185461400

exons_vis <- exons_df %>%
  filter(exon_end >= x_min - 500, exon_start <= x_max + 500)

## Label nudges for non-overlapping rsID placement
label_df <- stab %>%
  mutate(
    nudge_x = case_when(
      rsid == "rs35596415" ~ -420, rsid == "rs59319722" ~    0,
      rsid == "rs11132306" ~ -380, rsid == "rs7698680"  ~  360,
      rsid == "rs7699687"  ~ -380, rsid == "rs7699724"  ~  420,
      rsid == "rs11132309" ~ -420, TRUE ~ 0
    ),
    nudge_y = case_when(
      rsid == "rs35596415" ~  10, rsid == "rs59319722" ~  12,
      rsid == "rs11132306" ~ -10, rsid == "rs7698680"  ~  12,
      rsid == "rs7699687"  ~   9, rsid == "rs7699724"  ~  14,
      rsid == "rs11132309" ~   9, TRUE ~ 8
    )
  )

## Panel D top: Lollipop plot
gg_lollipop <- ggplot(stab, aes(x = pos, y = freq_train)) +
  geom_segment(aes(xend = pos, yend = 0), colour = "grey60", linewidth = 1.5) +
  geom_point(aes(size = freq_test), colour = "#0072B2") +
  geom_text_repel(
    data = label_df,
    aes(x = pos, y = freq_train, label = rsid),
    nudge_x = label_df$nudge_x, nudge_y = label_df$nudge_y,
    size = 10, fontface = "plain", box.padding = 0.35,
    point.padding = 0.25, min.segment.length = 0,
    segment.color = "black", segment.size = 0.5, max.overlaps = Inf
  ) +
  scale_size_continuous(name = "Test set validation frequency", range = c(3, 7),
                       breaks = c(3, 4, 5, 14), limits = c(3, 14)) +
  scale_y_continuous("Selection frequency in training sets\n(out of 100 iterations)",
                     limits = c(0, 100), breaks = seq(0, 100, 20)) +
  coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
  theme_bw(base_size = 20) +
  theme(
    axis.title.x       = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.title.y       = element_text(size = 22, face = "bold"),
    axis.text.y        = element_text(size = 22, colour = "black"),
    legend.position    = "right",
    legend.justification = c(0, 0.9),
    legend.background  = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
    legend.title       = element_text(size = 20),
    legend.text        = element_text(size = 20),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#AAAAAA", linewidth = 0.5, linetype = "dashed"),
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.margin        = margin(5, 40, 0, 5)
  )

## Panel D bottom: Gene structure
gene_line  <- data.frame(x = x_min - 300, xend = x_max + 300, y = 0.5)
intron_segs <- data.frame(
  start = c(x_min, 185460238, 185461159),
  end   = c(185458125, 185461048, x_max)
)

arrow_pos <- do.call(rbind, lapply(1:nrow(intron_segs), function(i) {
  s <- intron_segs$start[i] + 150
  e <- intron_segs$end[i] - 150
  if (e - s < 200) return(NULL)
  n <- max(2, round((e - s) / 400))
  data.frame(x = seq(s, e, length.out = n))
}))

gg_gene <- ggplot() +
  geom_segment(data = gene_line,
               aes(x = x, xend = xend, y = y, yend = y),
               colour = "black", linewidth = 1.2) +
  geom_segment(data = arrow_pos,
               aes(x = x + 120, xend = x - 120, y = 0.5, yend = 0.5),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               colour = "black", linewidth = 1.0) +
  geom_rect(data = exons_vis,
            aes(xmin = pmax(exon_start, x_min - 400),
                xmax = pmin(exon_end, x_max + 400),
                ymin = 0.15, ymax = 0.85),
            fill = "black", colour = "black", linewidth = 0.6) +
  geom_text(data = exons_vis,
            aes(x = (pmax(exon_start, x_min) + pmin(exon_end, x_max)) / 2,
                y = -0.1, label = paste0("Exon ", exon_num)),
            size = 7, fontface = "bold", colour = "black") +
  annotate("text", x = x_max + 200, y = 1.0, label = "5 prime",
           size = 8, fontface = "bold", hjust = 0) +
  annotate("text", x = x_min - 200, y = 1.0, label = "3 prime",
           size = 8, fontface = "bold", hjust = 1) +
  annotate("text", x = (x_min + x_max) / 2, y = 1.05,
           label = "CCDC110", fontface = "bold.italic", size = 8) +
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(-0.5, 1.5), clip = "off") +
  scale_x_continuous("Chromosome 4 position (bp)", labels = comma) +
  theme_bw(base_size = 22) +
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

## Combine lollipop + gene structure
p2d <- gg_lollipop / gg_gene + plot_layout(heights = c(4, 1.5))

ggsave("figures/Fig2d.tiff", p2d, width = 18, height = 7.5,
       units = "in", dpi = 800, compression = "lzw")
ggsave("figures/Fig2d.pdf", p2d, width = 14, height = 8, units = "in")
cat("  Saved figures/Fig2d.tiff and Fig2d.pdf\n")

cat("\nFigure 2 complete. All panels saved to figures/\n")