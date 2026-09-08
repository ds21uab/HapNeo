## ============================================================================
## 01_discovery_holdout_split.R
## Split ORIEN+SU2C cohort into Discovery (80%) and Hold-out validation (20%)
## Stratified by OS_status to maintain outcome balance
## ============================================================================

library(data.table)
library(dplyr)

# --- Load data ----------------------------------------------------------------
clinical <- read.csv("source_data/ORIEN_SU2C_CLINICAL.csv")
genotype <- fread("source_data/ORIEN_SU2C_GT.csv",
                  stringsAsFactors = FALSE,
                  strip.white = TRUE,
                  check.names = FALSE,
                  header = TRUE,
                  fill = TRUE)

# --- Split 80:20 stratified by OS_status -------------------------------------
set.seed(123)

dis_indices <- clinical %>%
  group_by(OS_status) %>%
  sample_frac(0.8) %>%
  pull(patient_id)

holdout_indices <- setdiff(clinical$patient_id, dis_indices)

# --- Subset clinical data -----------------------------------------------------
clinical_dis     <- clinical[clinical$patient_id %in% dis_indices, ]
clinical_holdout <- clinical[clinical$patient_id %in% holdout_indices, ]

cat("Discovery set: ", nrow(clinical_dis), "patients\n")
cat("Hold-out set:  ", nrow(clinical_holdout), "patients\n")

# Check OS_status balance
cat("\n--- OS_status counts ---\n")
cat("Discovery:\n"); print(table(clinical_dis$OS_status))
cat("Hold-out:\n");  print(table(clinical_holdout$OS_status))

cat("\n--- OS_status proportions ---\n")
cat("Discovery:\n"); print(prop.table(table(clinical_dis$OS_status)))
cat("Hold-out:\n");  print(prop.table(table(clinical_holdout$OS_status)))

# --- Subset and reorder genotype data -----------------------------------------
# Discovery genotype
geno_dis <- genotype[, c("CHROM", "POS", "REF", "ALT", "var_id", dis_indices), with = FALSE]
geno_dis_info    <- geno_dis[, 1:5]
geno_dis_samples <- geno_dis[, clinical_dis$patient_id, with = FALSE]
geno_dis_ordered <- cbind(geno_dis_info, geno_dis_samples)
stopifnot(all(clinical_dis$patient_id == colnames(geno_dis_ordered)[6:ncol(geno_dis_ordered)]))

# Hold-out genotype
geno_holdout <- genotype[, c("CHROM", "POS", "REF", "ALT", "var_id", holdout_indices), with = FALSE]
geno_holdout_info    <- geno_holdout[, 1:5]
geno_holdout_samples <- geno_holdout[, clinical_holdout$patient_id, with = FALSE]
geno_holdout_ordered <- cbind(geno_holdout_info, geno_holdout_samples)
stopifnot(all(clinical_holdout$patient_id == colnames(geno_holdout_ordered)[6:ncol(geno_holdout_ordered)]))

# --- Save outputs -------------------------------------------------------------
fwrite(clinical_dis, file = "source_data/clinical_discovery.csv",
       sep = ",", quote = FALSE, col.names = TRUE)
fwrite(clinical_holdout, file = "source_data/clinical_holdout.csv",
       sep = ",", quote = FALSE, col.names = TRUE)
fwrite(geno_dis_ordered, file = "source_data/geno_discovery.csv",
       sep = ",", quote = FALSE, col.names = TRUE)
fwrite(geno_holdout_ordered, file = "source_data/geno_holdout.csv",
       sep = ",", quote = FALSE, col.names = TRUE)

cat("\nDone. Files saved to source_data/\n")