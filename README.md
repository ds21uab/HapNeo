# HapNeo

Germline haplotype predicts survival after immunotherapy in NSCLC via neoantigen presentation

## REQUIREMENTS

- R (≥ 4.0)
- R packages: dplyr, survival, survminer, ggplot2, data.table, forestploter, writexl, ggpubr, pROC

## SUMMARY

HapNeo identifies a germline haplotype in CCDC110 associated with improved overall survival in immune checkpoint inhibitor (ICI)-treated non-small cell lung cancer (NSCLC).

The ORIEN and SU2C cohorts (n=387) were merged and split into a Discovery set and a Hold-out validation set. Germline variants were tested for association with overall survival using multivariate Cox proportional hazards regression, adjusting for age, sex, batch effects, and population ancestry (PC1–PC5). A two-step iterative resampling approach (100 iterations, stratified 80:20 splits) was used to ensure robustness of variant discovery. Within each iteration, training and test samples were non-overlapping. Variants were considered validated if they showed concordant hazard ratio direction and one-tailed 
significance (P < 0.05) in independent test sets. This identified seven germline variants in linkage disequilibrium, all mapping to CCDC110. The haplotype was further validated in two independent cohorts including Hold-out validation set and POPLAR clinical trial cohort.

The mechanistic basis involves the S409F missense variant generating a proteasomal cleavage site that produces a novel MHC-I-binding 9-mer peptide (LVKQGSIIF), enabling neoantigen presentation across multiple HLA alleles.

## REPOSITORY STRUCTURE

- `scripts/` — R scripts for figure generation
- `source_data/` — Source data files for manuscript figures

## CITATION

*Manuscript under review*

## AUTHORS

Divya Sahu, Aakrosh Ratan, Anindya Dutta
