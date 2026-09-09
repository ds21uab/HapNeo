## ============================================================================
## 03_validate_independent_cohorts.R
## Validate discovery-significant variants in two independent cohorts:
##   1. Hold-out validation set (20% of ORIEN+SU2C)
##   2. POPLAR clinical trial cohort (fully external)
## Validation criteria: concordant HR direction + one-tailed P < 0.05
## ============================================================================
##
## INPUT FILES:
##   Hold-out (controlled access — not included in repo):
##     - clinical_holdout.csv
##       Columns: patient_id, OS_month, OS_status, age, sex, batch, PC1-PC5
##     - geno_holdout.csv
##       Columns: CHROM, POS, REF, ALT, var_id, <patient_id columns>
##
##   POPLAR (controlled access — available from EGA: EGAS00001002460):
##     - clinical_poplar.csv
##       Columns: ANON_PATIENT_ID, AGE, GENDER, PC1-PC5, OS_MONTHS, OS_status
##     - geno_poplar.csv
##       Columns: CHROM, POS, REF, ALT, var_info, <patient_id columns>
##
##   Discovery results (from discovery_cox_gwas.R):
##     - results/ewas_discovery_results.csv
##
## OUTPUT:
##   - results/validation_holdout.csv
##   - results/validation_poplar.csv
##   - results/validation_summary.csv   (intersection across cohorts)
##
## See README.md for data access instructions.
## ============================================================================

library(data.table)
library(dplyr)
library(survival)
library(future)
library(future.apply)

# --- Load discovery significant variants --------------------------------------
lead_snps     <- fread("results/ewas_discovery_results.csv")
sig_threshold <- 5e-4
sig_snps      <- lead_snps[!is.na(P) & P < sig_threshold]
selected_snps <- sig_snps$SNP

cat("Discovery significant variants (P <", sig_threshold, "):", nrow(sig_snps), "\n\n")

# --- Cox regression function --------------------------------------------------
run_cox <- function(i, snps, geno_mat, merged_df, covariates) {
  out <- data.frame(
    SNP = snps[i], Beta = NA, HR = NA, SE = NA,
    Z = NA, P = NA, Death = NA
  )
  tryCatch({
    g <- geno_mat[, i]
    valid <- !is.na(g)
    g_val <- g[valid]
    df <- merged_df[valid, ]
    df$GENO <- as.numeric(g_val)

    n_death <- sum(df$OS_status == 1, na.rm = TRUE)
    out$Death <- n_death

    if (n_death <= 10) return(out)
    if (length(unique(g_val)) <= 1) return(out)

    formula <- as.formula(paste(
      "Surv(OS_month, OS_status) ~ GENO +",
      paste(covariates, collapse = " + ")
    ))

    fit <- coxph(formula, data = df)
    ss  <- summary(fit)$coefficients["GENO", ]

    out$Beta <- ss["coef"]
    out$HR   <- ss["exp(coef)"]
    out$SE   <- ss["se(coef)"]
    out$Z    <- ss["z"]
    out$P    <- ss["Pr(>|z|)"]
    return(out)
  }, error = function(e) {
    return(out)
  })
}

# --- Helper: prepare genotype matrix and run validation -----------------------
validate_cohort <- function(geno_file, clinical_file, selected_snps,
                            covariates, cohort_name) {

  cat("=== Validating in", cohort_name, "===\n")

  clinical <- fread(clinical_file, check.names = FALSE)
  geno     <- fread(geno_file, check.names = FALSE)

  cat("Patients:", nrow(clinical), "\n")

  # Standardize column names across cohorts
  if ("var_info" %in% colnames(geno))            setnames(geno, "var_info", "var_id")
 if ("ID" %in% colnames(geno))            setnames(geno, "ID", "var_id") 
 if ("ANON_PATIENT_ID" %in% colnames(clinical)) setnames(clinical, "ANON_PATIENT_ID", "patient_id")
  if ("AGE" %in% colnames(clinical))             setnames(clinical, "AGE", "age")
  if ("GENDER" %in% colnames(clinical))          setnames(clinical, "GENDER", "sex")
  if ("OS_MONTHS" %in% colnames(clinical))       setnames(clinical, "OS_MONTHS", "OS_month")
  if ("OS_months" %in% colnames(clinical))       setnames(clinical, "OS_months", "OS_month")

  # Align patient order
  geno_ids     <- colnames(geno)[6:ncol(geno)]
  clinical_ids <- clinical$patient_id

  reorder_idx <- match(clinical_ids, geno_ids)
  if (any(is.na(reorder_idx))) {
    # Keep only overlapping patients
    common_ids  <- intersect(clinical_ids, geno_ids)
    clinical    <- clinical[patient_id %in% common_ids]
    reorder_idx <- match(clinical$patient_id, geno_ids)
    cat("Patients after overlap:", nrow(clinical), "\n")
  }

  geno <- geno[, c(1:5, (5 + reorder_idx)), with = FALSE]
  stopifnot(all(clinical$patient_id == colnames(geno)[6:ncol(geno)]))

  # Subset to discovery-significant variants
  geno <- geno[geno$var_id %in% selected_snps, ]
  cat("Variants to validate:", nrow(geno), "\n")

  if (nrow(geno) == 0) {
    cat("No matching variants found in", cohort_name, "\n\n")
    return(NULL)
  }

  # Prepare matrix
  geno_mat <- t(as.matrix(geno[, 6:ncol(geno)]))
  colnames(geno_mat) <- geno$var_id
  geno_df <- as.data.frame(geno_mat)
  geno_df$patient_id <- rownames(geno_df)

  merged <- merge(clinical, geno_df, by = "patient_id")

  # Scale PCs for numerical stability
  for (pc in covariates[grep("^PC", covariates)]) {
    if (pc %in% colnames(merged)) {
      merged[[pc]] <- scale(merged[[pc]])
    }
  }

  snp_cols <- intersect(geno$var_id, colnames(merged))
  geno_matrix <- as.matrix(merged[, ..snp_cols, drop = FALSE])
  rownames(geno_matrix) <- merged$patient_id
  snp_list <- snp_cols

  # Run Cox regression
  plan(multisession, workers = 10)

  res <- future_lapply(
    seq_along(snp_list),
    run_cox,
    snps = snp_list,
    geno_mat = geno_matrix,
    merged_df = merged,
    covariates = covariates
  )

  res <- rbindlist(res)
  cat("Variants tested:", nrow(res), "\n\n")
  return(res)
}

# --- Helper: merge with discovery and assess concordance ----------------------
assess_validation <- function(res_val, sig_snps, suffix) {

  colnames(res_val) <- paste0(c("SNP", "Beta", "HR", "SE", "Z", "P", "Death"),
                               "_", suffix)
  setnames(res_val, paste0("SNP_", suffix), "SNP")

  disc_set <- copy(sig_snps)
  colnames(disc_set) <- c("SNP", "Beta_disc", "HR_disc", "SE_disc",
                           "Z_disc", "P_disc", "Death_disc")

  val <- merge(disc_set, res_val, by = "SNP", all = FALSE)

  # Concordant HR direction
  beta_val <- val[[paste0("Beta_", suffix)]]
  z_val    <- val[[paste0("Z_", suffix)]]

  val[, direction_match := sign(Beta_disc) == sign(beta_val)]

  # One-tailed P (only meaningful if direction matches)
  val[, p_one_tailed := fifelse(
    direction_match & !is.na(z_val),
    pnorm(-abs(z_val)),
    NA_real_
  )]

  val[, validated := direction_match == TRUE & p_one_tailed < 0.05]

  cat("Direction match:   ", sum(val$direction_match, na.rm = TRUE), "/", nrow(val), "\n")
  cat("Validated (P<0.05):", sum(val$validated, na.rm = TRUE), "/", nrow(val), "\n")

  return(val)
}


## ============================================================================
## PART 1: Hold-out validation
## ============================================================================

res_holdout <- validate_cohort(
  geno_file     = "geno_holdout.csv",
  clinical_file = "clinical_holdout.csv",
  selected_snps = selected_snps,
  covariates    = c("age", "sex", "batch", "PC1", "PC2", "PC3", "PC4", "PC5"),
  cohort_name   = "Hold-out set"
)

if (!is.null(res_holdout) && nrow(res_holdout) > 0) {
  val_holdout <- assess_validation(res_holdout, sig_snps, suffix = "holdout")
} else {
  cat("Hold-out validation skipped — no matching variants.\n")
  val_holdout <- NULL
}

## ============================================================================
## PART 2: POPLAR validation
## ============================================================================
## NOTE: POPLAR data not included in repo (EGA controlled access).
## To run this section, obtain data from EGA (EGAS00001002460) and place
## clinical_poplar.csv and geno_poplar.csv in the working directory.
## POPLAR model does not include batch (single cohort).
## ============================================================================

res_poplar <- validate_cohort(
  geno_file     = "geno_poplar.csv",
  clinical_file = "clinical_poplar.csv",
  selected_snps = selected_snps,
  covariates    = c("age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5"),
  cohort_name   = "POPLAR"
)

if (!is.null(res_poplar) && nrow(res_poplar) > 0) {
  val_poplar <- assess_validation(res_poplar, sig_snps, suffix = "poplar")
} else {
  cat("POPLAR validation skipped — no matching variants.\n")
  cat("Check if var_id format in geno_poplar.csv matches discovery results.\n")
  val_poplar <- NULL
}

## ============================================================================
## PART 3: Cross-cohort summary
## ============================================================================

snps_disc    <- sig_snps$SNP
snps_holdout <- if (!is.null(val_holdout)) val_holdout[validated == TRUE]$SNP else character(0)
snps_poplar  <- if (!is.null(val_poplar))  val_poplar[validated == TRUE]$SNP  else character(0)

cat("\n========================================\n")
cat("CROSS-COHORT VALIDATION SUMMARY\n")
cat("========================================\n")
cat("Discovery significant (P <", sig_threshold, "):", length(snps_disc), "\n")
cat("Validated in Hold-out:  ", length(snps_holdout), "\n")
cat("Validated in POPLAR:    ", length(snps_poplar), "\n")

# Intersection: validated in BOTH independent cohorts
snps_both <- intersect(snps_holdout, snps_poplar)
cat("Validated in BOTH:      ", length(snps_both), "\n")

if (length(snps_both) > 0) {
  cat("\nVariants validated in both cohorts:\n")
  print(snps_both)

  summary_table <- merge(
    val_holdout[SNP %in% snps_both, .(SNP, Beta_disc, HR_disc, P_disc,
                                       Beta_holdout, HR_holdout, p_one_tailed)],
    val_poplar[SNP %in% snps_both, .(SNP, Beta_poplar, HR_poplar,
                                      p_one_tailed_poplar = p_one_tailed)],
    by = "SNP"
  )
  cat("\n")
  print(summary_table)
}

# --- Save results -------------------------------------------------------------
dir.create("results", showWarnings = FALSE)

if (!is.null(val_holdout)) {
  fwrite(val_holdout, file = "results/validation_holdout.csv",
         sep = ",", quote = FALSE, col.names = TRUE)
}

if (!is.null(val_poplar)) {
  fwrite(val_poplar, file = "results/validation_poplar.csv",
         sep = ",", quote = FALSE, col.names = TRUE)
}

# Combined summary
if (!is.null(val_holdout) && !is.null(val_poplar)) {
  summary_all <- merge(
    val_holdout[, .(SNP, Beta_disc, HR_disc, P_disc,
                    Beta_holdout, HR_holdout, p_one_tailed_holdout = p_one_tailed,
                    validated_holdout = validated)],
    val_poplar[, .(SNP, Beta_poplar, HR_poplar, p_one_tailed_poplar = p_one_tailed,
                   validated_poplar = validated)],
    by = "SNP", all = TRUE
  )
  summary_all[, validated_both := validated_holdout == TRUE & validated_poplar == TRUE]
  fwrite(summary_all, file = "results/validation_summary.csv",
         sep = ",", quote = FALSE, col.names = TRUE)
} else if (!is.null(val_holdout)) {
  fwrite(val_holdout, file = "results/validation_summary.csv",
         sep = ",", quote = FALSE, col.names = TRUE)
} else if (!is.null(val_poplar)) {
  fwrite(val_poplar, file = "results/validation_summary.csv",
         sep = ",", quote = FALSE, col.names = TRUE)
}

cat("\nDone. Results saved to results/\n")
