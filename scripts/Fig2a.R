#!/usr/bin/Rscript

## Copyright Divya Sahu, 2026
## ============================================================
## Figure 2a: Manhattan plot — Discovery cohort GWAS
## CCDC110 germline haplotype manuscript
## Input:  Fig2a.csv (Source Data Fig. 2a)
## Output: Fig2a.tiff
## R version: 4.5.2

## ============================================================

## ---- Load libraries -----------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(ggtext)

## ---- Read GWAS summary statistics ---------------------------------------
f1 <- read.csv("Fig2a.csv", sep = ",")

## ---- Define lead and LD SNPs -------------------------------------------
## Lead SNP identified from bootstrap cross-validated GWAS (rs7698680)
lead_snp   <- "chr4_185459692_A_T"
lead_label <- "CCDC110 (rs7698680)"

## SNPs in linkage disequilibrium with rs7698680 (LDhap/PLINK2)
ccdc110_ld_snps <- c(
  "chr4_185458745_G_C",
  "chr4_185461052_A_G",
  "chr4_185459089_A_C",
  "chr4_185459361_G_A",
  "chr4_185459961_G_T",
  "chr4_185460011_G_A"
)

## ---- Prepare GWAS data for plotting -------------------------------------
## Replace missing or zero P-values with machine minimum to avoid -log10 errors
gwas2 <- f1 %>%
  mutate(P = ifelse(is.na(P) | P == 0, .Machine$double.xmin, P)) %>%
  filter(!is.na(P)) %>%
  mutate(
    CHR  = as.integer(str_match(SNP, "^chr([0-9]+)_")[, 2]),
    BP   = as.integer(str_match(SNP, "^chr[0-9]+_([0-9]+)_")[, 2]),
    LOGP = -log10(P)
  ) %>%
  filter(!is.na(CHR), !is.na(BP)) %>%
  arrange(CHR, BP)

## ---- Compute cumulative chromosome positions for x-axis -----------------
chr_pos <- gwas2 %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(
    chr_len   = as.numeric(chr_len),
    chr_start = dplyr::lag(cumsum(chr_len), default = 0)
  )

## Add cumulative base-pair position and chromosome parity (for alternating colors)
gwas2 <- gwas2 %>%
  left_join(chr_pos, by = "CHR") %>%
  mutate(
    BP_cum     = as.numeric(BP) + chr_start,
    chr_parity = CHR %% 2
  )

## Chromosome label positions (centered on each chromosome)
axis_df <- chr_pos %>%
  mutate(center = chr_start + chr_len / 2)

## ---- Define significance thresholds -------------------------------------
## Suggestive threshold justified by power analysis (powerSurvEpi; see Methods)
thresh_suggestive <- 5e-4
## Standard genome-wide significance threshold
thresh_gw         <- 5e-8

## ---- Generate Manhattan plot --------------------------------------------
p_manh_alt <- ggplot(gwas2, aes(x = BP_cum, y = LOGP)) +

  ## All SNPs, colored by chromosome parity (grey/blue alternation)
  geom_point(
    aes(color = factor(chr_parity)),
    size = 0.45, alpha = 0.85
  ) +
  scale_color_manual(
    values = c("grey60", "#0072B2"),
    guide  = "none"
  ) +

  ## Highlight SNPs in LD with lead SNP (black crosses)
  geom_point(
    data   = subset(gwas2, SNP %in% ccdc110_ld_snps),
    aes(shape = "SNPs in LD with rs7698680"),
    color  = "black",
    size   = 2.5,
    stroke = 1
  ) +

  ## Highlight lead SNP rs7698680 (red filled circle)
  geom_point(
    data  = subset(gwas2, SNP == lead_snp),
    aes(shape = "Lead SNP rs7698680"),
    color = "red",
    size  = 2.5
  ) +

  ## Suggestive significance threshold line (P < 5e-4)
  geom_hline(
    aes(yintercept = -log10(thresh_suggestive),
        linetype   = "thresh_line"),
    color     = "deeppink",
    linewidth = 0.8
  ) +

  ## Genome-wide significance threshold line (P < 5e-8, no legend entry)
  geom_hline(
    yintercept = -log10(thresh_gw),
    color      = "grey40",
    linetype   = "dotted",
    linewidth  = 0.8
  ) +

  ## Shape legend: lead SNP (red circle) and LD SNPs (black cross)
  scale_shape_manual(
    name   = NULL,
    values = c("Lead SNP rs7698680"       = 16,
               "SNPs in LD with rs7698680" = 4),
    guide  = guide_legend(order = 1,
                          override.aes = list(
                            color  = c("red", "black"),
                            size   = c(2.5, 2.5),
                            stroke = c(0, 1)
                          ))
  ) +

  ## Linetype legend: suggestive threshold line
  scale_linetype_manual(
    name   = NULL,
    values = c("thresh_line" = "solid"),
    labels = "*P* < 5 × 10<sup>−4</sup>",
    guide  = guide_legend(order = 2,
                          override.aes = list(color = "deeppink"))
  ) +

  ## Label lead SNP with gene name
  geom_text_repel(
    data = subset(gwas2, SNP == lead_snp),
    aes(label = lead_label),
    color = "red", size = 6, fontface = "bold",
    nudge_y = 0.6, min.segment.length = 0
  ) +

  ## X-axis: chromosome numbers centered on each chromosome
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR,
    expand = expansion(mult = c(0.01, 0.01))
  ) +

  ## Y-axis: -log10(P)
  scale_y_continuous(
    bquote(bold(-log[10](P))),
    breaks = seq(0, 8, 2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  xlab("Chromosome") +
  ggtitle("Discovery cohort") +

  ## Theme and formatting
  theme_bw(base_size = 11) +
  theme(
    plot.title            = element_text(face = "bold.italic", hjust = 0.5, size = 20),
    axis.title            = element_text(face = "bold"),
    axis.title.x          = element_text(size = 17, face = "bold"),
    axis.title.y          = element_text(size = 17, face = "bold"),
    axis.text.x           = element_text(size = 15, angle = 0, hjust = 0.5),
    axis.text.y           = element_text(size = 15),
    panel.grid.minor      = element_blank(),
    panel.grid.major.x    = element_blank(),
    panel.grid.major.y    = element_blank(),
    panel.border          = element_rect(colour = "black", fill = NA),
    legend.position       = c(0.85, 0.85),
    legend.background     = element_blank(),
    legend.box.background = element_blank(),
    legend.text           = element_markdown(size = 12, face = "bold"),
    legend.key            = element_blank(),
    legend.key.height     = unit(0.4, "cm"),
    legend.key.width      = unit(0.86, "cm"),
    legend.margin         = margin(0, 0, 0, 0)
  )

## Display plot
#p_manh_alt

## ---- Save high-resolution TIFF (Nature Medicine specifications) ---------
ggsave(
  "Fig2a.tiff",
  p_manh_alt,
  width = 12, height = 5,
  units = "in", dpi = 800, compression = "lzw"
)