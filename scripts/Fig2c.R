#!/usr/bin/Rscript
## Copyright Divya Sahu, 2025
## ============================================================
## Figure 2c: Venn Diagram — Testing set, Hold-out external, and POPLAR
## CCDC110 germline haplotype manuscript
## Input:  Source_Data_Fig2c.xlsx (Source Data Fig. 2c)
## Output: Fig2c.pdf
## R version: 4.5.2
## Key packages: ggVennDiagram, readxl, sf
## ============================================================

## ---- Load libraries -----------------------------------------------------
library(readxl)
library(ggVennDiagram)
library(sf)

## ---- Read source data (one tab per cohort) ------------------------------
test_snps <- read_excel("Source_Data_Fig2c.xlsx", sheet = "Testing set")
ext_snps  <- read_excel("Source_Data_Fig2c.xlsx", sheet = "Hold-out external")
poplar    <- read_excel("Source_Data_Fig2c.xlsx", sheet = "POPLAR")

## ---- Apply cohort-specific significance filters -------------------------
## Testing set: already filtered via bootstrap cross-validation (see Methods)
## Hold-out external: same effect direction and one-tailed P < 0.05
ext_snps <- ext_snps[ext_snps$direction_match == TRUE & ext_snps$p_one_tailed_ext <= 0.05, ]
## POPLAR: same effect direction and one-tailed P < 0.06
poplar   <- poplar[poplar$Direction == "same" & poplar$P_one_tailed <= 0.06, ]

## ---- Build Venn diagram -------------------------------------------------
x <- list(
  Test_set            = test_snps$SNP,
  POPLAR              = poplar$SNP,
  `Hold-out external` = ext_snps$SNP
)

## ---- Export Source Data for Fig 2b --------------------------------------
library(writexl)

write_xlsx(
  list(
    "Testing set"          = data.frame(SNP = x[["Test_set"]]),
    "Hold-out external"  = data.frame(SNP = x[["Hold-out external"]]),
    "POPLAR"             = data.frame(SNP = x[["POPLAR"]])
  ),
  "Raw_Data_Fig2c.xlsx"
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
  annotate("text", x = -4.5, y = 6.5, label = "Testing set",
           fontface = "bold.italic", size = 10) +
  annotate("text", x = 8.5, y = 6.5, label = "POPLAR",
           fontface = "bold.italic", size = 10) +
  annotate("text", x = 2.0, y = -10.5, label = "Hold-out external",
           fontface = "bold.italic", size = 10)

## Display plot
p_venn

## ---- Save high-resolution PDF -------------------------------------------
ggsave("Fig2c.pdf", plot = p_venn, width = 7, height = 5, dpi = 800)