#!/usr/bin/Rscript
## Copyright Divya Sahu, 2025
## ============================================================
## Figure 2b: Venn Diagram — Discovery, Hold-out external, and POPLAR
## CCDC110 germline haplotype manuscript
## Input:  Source_Data_Fig2b.xlsx (Source Data Fig. 2b)
## Output: Fig2b.pdf
## R version: 4.5.2
## Key packages: ggVennDiagram, readxl, sf, ggplot2
## ============================================================

## ---- Load libraries -----------------------------------------------------
library(readxl)
library(ggVennDiagram)
library(sf)
library(ggplot2)

## ---- Read source data (one tab per cohort) ------------------------------
discovery_snps <- read_excel("Source_Data_Fig2b.xlsx", sheet = "Discovery")
ext_snps       <- read_excel("Source_Data_Fig2b.xlsx", sheet = "Hold-out external")
poplar         <- read_excel("Source_Data_Fig2b.xlsx", sheet = "POPLAR")

## ---- Apply cohort-specific significance filters -------------------------
## Discovery: P < 5e-4 (suggestive threshold; see Methods)
discovery_snps <- discovery_snps[discovery_snps$P <= 5e-4, ]
## Hold-out external: same effect direction and one-tailed P < 0.05
ext_snps       <- ext_snps[ext_snps$direction_match == TRUE & ext_snps$p_one_tailed_ext <= 0.05, ]
## POPLAR: same effect direction and one-tailed P < 0.05
poplar         <- poplar[poplar$Direction == "same" & poplar$P_one_tailed <= 0.05, ]

## ---- Build Venn diagram -------------------------------------------------
x <- list(
  Discovery           = discovery_snps$SNP,
  POPLAR              = poplar$SNP,
  `Hold-out external` = ext_snps$SNP
)

## ---- Export Source Data for Fig 2b --------------------------------------
library(writexl)

write_xlsx(
  list(
    "Discovery"          = data.frame(SNP = x[["Discovery"]]),
    "Hold-out external"  = data.frame(SNP = x[["Hold-out external"]]),
    "POPLAR"             = data.frame(SNP = x[["POPLAR"]])
  ),
  "Source_Data_Fig2b.xlsx"
)

p_venn <- ggVennDiagram(x,
  label       = "count",
  label_alpha = 0,
  label_size  = 7,
  set_size    = 0
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_color_manual(values = c("grey40", "grey40", "grey40")) +
  theme(
    plot.margin     = margin(60, 40, 60, 40),
    legend.position = "none"
  ) +
  coord_sf(clip = "off") +
  annotate("text", x = -4.5, y = 6.5, label = "Discovery",
           fontface = "bold.italic", size = 10) +
  annotate("text", x = 8.5, y = 6.5, label = "POPLAR",
           fontface = "bold.italic", size = 10) +
  annotate("text", x = 2.0, y = -10.5, label = "Hold-out external",
           fontface = "bold.italic", size = 10)

## Display plot
p_venn

## ---- Save high-resolution PDF -------------------------------------------
ggsave("Fig2b.pdf", plot = p_venn, width = 7, height = 5, dpi = 800)
