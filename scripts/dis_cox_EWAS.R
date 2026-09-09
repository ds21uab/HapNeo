#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(survival)
library(future)
library(future.apply)

# --- Load data ----------------------------------------------------------------
clinical <- fread("source_data/clinical_discovery.csv", check.names = FALSE)
geno     <- fread("source_data/geno_discovery.csv", check.names = FALSE)
if ("ID" %in% colnames(geno)) setnames(geno, "ID", "var_id")

stopifnot(all(clinical$patient_id == colnames(geno)[6:ncol(geno)]))
cat("Loaded", nrow(clinical), "patients and", nrow(geno), "variants\n")

# --- Setup parallel -----------------------------------------------------------
plan(multisession, workers = 10)
sig_threshold <- 5e-4

# --- Cox regression function --------------------------------------------------
# For each variant, fits: Surv(OS_month, OS_status) ~ GENO + covariates
# Returns beta, HR, SE, Z, P, and number of events

run_cox <- function(i, snps, geno_mat, merged_df) {
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

    # Skip if too few events or no genotype variation
    if (n_death <= 10) return(out)
    if (length(unique(g_val)) <= 1) return(out)

    fit <- coxph(
      Surv(OS_month, OS_status) ~ GENO + age + sex +
        PC1 + PC2 + PC3 + PC4 + PC5 + batch,
      data = df
    )

    ss <- summary(fit)$coefficients["GENO", ]
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

# --- Prepare genotype matrix --------------------------------------------------
geno_mat <- t(as.matrix(geno[, 6:ncol(geno)]))
colnames(geno_mat) <- geno$var_id
geno_df <- as.data.frame(geno_mat)
geno_df$patient_id <- rownames(geno_df)

merged <- merge(clinical, geno_df, by = "patient_id")

snp_cols <- setdiff(
  colnames(merged),
  c("patient_id", "age", "sex", "capture_panel", "capture_panel_short",
    "batch", paste0("PC", 1:10), "histology", "tumor_stage",
    "OS_status", "OS_month")
)

geno_matrix <- as.matrix(merged[, ..snp_cols, drop = FALSE])
rownames(geno_matrix) <- merged$patient_id

# --- Filter variants: require ≥10 ref and ≥10 alt carriers -------------------
n_ref <- colSums(geno_matrix == 0, na.rm = TRUE)
n_alt <- colSums((geno_matrix == 1) | (geno_matrix == 2), na.rm = TRUE)
keep_snps <- which(n_ref >= 10 & n_alt >= 10)

cat("Variants before filtering:", length(snp_cols), "\n")
cat("Variants after filtering: ", length(keep_snps), "\n")

if (length(keep_snps) == 0) {
  cat("No valid variants after filtering\n")
  quit(save = "no", status = 0)
}

geno_matrix <- geno_matrix[, keep_snps, drop = FALSE]
snp_list <- snp_cols[keep_snps]

# --- Run Cox regression across all variants -----------------------------------
cat("Running Cox regression on", length(snp_list), "variants...\n")

res_discovery <- future_lapply(
  seq_along(snp_list),
  run_cox,
  snps = snp_list,
  geno_mat = geno_matrix,
  merged_df = merged
)

res_discovery <- rbindlist(res_discovery)

# --- Summary ------------------------------------------------------------------
sig_snps <- res_discovery[!is.na(P) & P < sig_threshold]
cat("\nVariants tested:", nrow(res_discovery), "\n")
cat("Variants with P <", sig_threshold, ":", nrow(sig_snps), "\n")

# --- Save results -------------------------------------------------------------
dir.create("results", showWarnings = FALSE)
fwrite(res_discovery, file = "results/ewas_discovery_results.csv",
       sep = ",", quote = FALSE, col.names = TRUE)
fwrite(sig_snps, file = "results/ewas_discovery_significant.csv",
       sep = ",", quote = FALSE, col.names = TRUE)

cat("\nDone. Results saved to results/\n")
