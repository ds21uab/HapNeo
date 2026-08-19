# HapNeo

Predict ICI response from germline haplotype driven neoantigen presentation

## REQUIREMENTS

- R (≥ 4.0)
- R packages: survival, survminer, ggplot2, data.table, forestploter, writexl, ggpubr

## SUMMARY

HapNeo identifies a germline CCDC110 haplotype that predicts improved overall survival in immune checkpoint inhibitor (ICI)-treated non-small cell lung cancer (NSCLC).

The ORIEN and SU2C cohorts (n=387) were split into a Discovery set and a Hold-out external validation set. Seven germline variants in linkage disequilibrium (lead SNP rs7698680) were identified through a discovery and validation approach, with additional bootstrap cross-validation (100 iterations, 80:20 splits). The haplotype was further validated in an independent POPLAR cohort.

The mechanistic basis involves the S409F missense variant generating a proteasomal cleavage site that produces a novel MHC-I-binding 9-mer peptide (LVKQGSIIF), enabling neoantigen presentation across multiple HLA alleles.

## REPOSITORY STRUCTURE

- `scripts/` — R scripts for figure generation
- `source_data/` — Source data files for manuscript figures

## CITATION

*Manuscript under review*

## AUTHORS

Divya Sahu, Aakrosh Ratan, Anindya Dutta
