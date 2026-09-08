## ============================================================================
## generate_resampling_splits.R
## Generate 100 stratified 80:20 train/test splits within the Discovery set
## for the two-step iterative resampling and validation approach
## ============================================================================
##
## INPUT FILES (controlled access — not included in repo):
##   - source_data/clinical_discovery.csv
##     Columns: patient_id, OS_month, OS_status, age, sex, batch, PC1-PC5
##   - source_data/geno_discovery.csv
##     Columns: CHROM, POS, REF, ALT, var_id, <patient_id columns>
##
## OUTPUT:
##   - splits/train_split_001.clinical.csv ... splits/train_split_100.clinical.csv
##   - splits/test_split_001.clinical.csv  ... splits/test_split_100.clinical.csv
##   - splits/train_split_001.geno.csv     ... splits/train_split_100.geno.csv
##   - splits/test_split_001.geno.csv      ... splits/test_split_100.geno.csv
##   (400 files total)
##
## Each iteration uses a different seed (123 + i) for reproducibility.
## Splits are stratified by OS_status to maintain outcome balance.
## Training and test samples are non-overlapping within each iteration.
## ============================================================================

library(data.table)
library(dplyr)

# --- Load data ----------------------------------------------------------------
clinical <- fread("source_data/clinical_discovery.csv", check.names = FALSE)
geno     <- fread("source_data/geno_discovery.csv", check.names = FALSE)

stopifnot(all(clinical$patient_id == colnames(geno)[6:ncol(geno)]))
cat("Discovery set:", nrow(clinical), "patients,", nrow(geno), "variants\n\n")

# --- Generate 100 stratified splits -------------------------------------------
n_splits <- 100
dir.create("splits", showWarnings = FALSE)

for (i in 1:n_splits) {
  set.seed(123 + i)

  # Stratified 80:20 split by OS_status
  train_ids <- clinical %>%
    group_by(OS_status) %>%
    sample_frac(0.8) %>%
    ungroup() %>%
    pull(patient_id)

  test_ids <- setdiff(clinical$patient_id, train_ids)

  # Verify no overlap
  stopifnot(length(intersect(train_ids, test_ids)) == 0)

  # Subset clinical
  clinical_train <- clinical[patient_id %in% train_ids]
  clinical_test  <- clinical[patient_id %in% test_ids]

  # Subset genotype (columns 1:5 = variant info, rest = patient genotypes)
  geno_train <- geno[, c("CHROM", "POS", "REF", "ALT", "var_id", train_ids), with = FALSE]
  geno_test  <- geno[, c("CHROM", "POS", "REF", "ALT", "var_id", test_ids), with = FALSE]

  # Save
  fwrite(clinical_train, sprintf("splits/train_split_%03d.clinical.csv", i))
  fwrite(clinical_test,  sprintf("splits/test_split_%03d.clinical.csv", i))
  fwrite(geno_train,     sprintf("splits/train_split_%03d.geno.csv", i))
  fwrite(geno_test,      sprintf("splits/test_split_%03d.geno.csv", i))

  if (i %% 10 == 0) {
    cat(sprintf("Split %3d/%d — train: %d, test: %d\n",
                i, n_splits, length(train_ids), length(test_ids)))
  }
}

cat("\nDone.", n_splits, "splits saved to splits/\n")