#!/usr/bin/env Rscript
## ============================================================================
## iterative_resampling_ewas.R
## Two-step iterative resampling and validation approach
##
## For each of 100 iterations:
##   Step 1: Cox GWAS on training set (80%) — identify variants with P < 5e-4
##   Step 2: Validate in test set (20%) — concordant HR + one-tailed P < 0.05
## Then aggregate across all iterations to identify robustly replicated variants.
## ============================================================================
##
## USAGE:
##   Option A — SLURM array (recommended for HPC):
##     sbatch 05_iterative_resampling_gwas.sh     # runs 100 iterations
##     Rscript 05_iterative_resampling_gwas.R aggregate  # then aggregate
##
##   Option B — Single iteration:
##     Rscript 05_iterative_resampling_gwas.R 42
##
##   Option C — All iterations sequentially (slow, for testing/reproducibility):
##     Rscript 05_iterative_resampling_gwas.R all
##
## INPUT FILES (from 04_generate_resampling_splits.R):
##   - splits/train_split_XXX.clinical.csv
##   - splits/train_split_XXX.geno.csv
##   - splits/test_split_XXX.clinical.csv
##   - splits/test_split_XXX.geno.csv
##
## OUTPUT:
##   Per iteration:
##     - results/resampling/train_iter_XXX.csv
##     - results/resampling/test_iter_XXX.csv
##   Aggregated:
##     - results/resampling_all_replicated.csv
##     - results/resampling_summary.csv
##     - results/resampling_frequency_plot.png
##
## See README.md for data access instructions.
## ============================================================================

library(data.table)
library(dplyr)
library(survival)
library(future)
library(future.apply)

# --- Settings -----------------------------------------------------------------
n_iters       <- 100
sig_threshold <- 5e-4
val_threshold <- 0.05

dir.create("results/resampling", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

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

align_geno_to_clinical <- function(geno, clinical, label) {
  geno_ids     <- colnames(geno)[6:ncol(geno)]
  clinical_ids <- clinical$patient_id
  reorder_idx  <- match(clinical_ids, geno_ids)

  if (any(is.na(reorder_idx))) {
    stop(paste("Mismatch: patient_ids in", label, "not found in genotype."))
  }

  geno <- geno[, c(1:5, (5 + reorder_idx)), with = FALSE]
  stopifnot(all(clinical$patient_id == colnames(geno)[6:ncol(geno)]))
  return(geno)
}

prepare_geno_matrix <- function(geno, clinical) {
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

  list(merged = merged, snp_cols = snp_cols)
}

# ==============================================================================
# RUN ONE ITERATION: train GWAS → test validation → directional concordance
# ==============================================================================
run_one_iteration <- function(iter) {

  cat("\n========== ITERATION:", iter, "==========\n")
  iter_label <- sprintf("%03d", iter)
  set.seed(123 + iter)

  # Load split data
  clinical_train <- fread(sprintf("splits/train_split_%s.clinical.csv", iter_label), check.names = FALSE)
  clinical_test  <- fread(sprintf("splits/test_split_%s.clinical.csv", iter_label), check.names = FALSE)
  geno_train     <- fread(sprintf("splits/train_split_%s.geno.csv", iter_label), check.names = FALSE)
  geno_test      <- fread(sprintf("splits/test_split_%s.geno.csv", iter_label), check.names = FALSE)

  cat("  Train:", nrow(clinical_train), "| Test:", nrow(clinical_test), "patients\n")

  geno_train <- align_geno_to_clinical(geno_train, clinical_train, "train")
  geno_test  <- align_geno_to_clinical(geno_test, clinical_test, "test")

  # --- STEP 1: Cox GWAS on training set (80%) ---------------------------------
  train_prep   <- prepare_geno_matrix(geno_train, clinical_train)
  merged_train <- train_prep$merged
  snp_cols     <- train_prep$snp_cols

  geno_matrix <- as.matrix(merged_train[, ..snp_cols, drop = FALSE])
  rownames(geno_matrix) <- merged_train$patient_id

  n_ref <- colSums(geno_matrix == 0, na.rm = TRUE)
  n_alt <- colSums((geno_matrix == 1) | (geno_matrix == 2), na.rm = TRUE)
  keep_snps <- which(n_ref >= 10 & n_alt >= 10)

  if (length(keep_snps) == 0) {
    cat("  No valid variants after filtering\n")
    return(NULL)
  }

  geno_matrix <- geno_matrix[, keep_snps, drop = FALSE]
  snp_list <- snp_cols[keep_snps]

  cat("  Training: testing", length(snp_list), "variants...\n")

  res_train <- future_lapply(
    seq_along(snp_list), run_cox,
    snps = snp_list, geno_mat = geno_matrix, merged_df = merged_train
  )
  res_train <- rbindlist(res_train)

  sig_snps <- res_train[!is.na(P) & P < sig_threshold]
  cat("  Training significant (P <", sig_threshold, "):", nrow(sig_snps), "\n")

  fwrite(res_train, sprintf("results/resampling/train_iter_%s.csv", iter_label))

  if (nrow(sig_snps) == 0) return(NULL)

  # --- STEP 2: Validate in test set (20%) -------------------------------------
  test_prep   <- prepare_geno_matrix(geno_test, clinical_test)
  merged_test <- test_prep$merged

  snp_cols_test <- intersect(sig_snps$SNP, colnames(merged_test))

  if (length(snp_cols_test) == 0) {
    cat("  No significant variants found in test set\n")
    return(NULL)
  }

  geno_test_matrix <- as.matrix(merged_test[, ..snp_cols_test, drop = FALSE])
  rownames(geno_test_matrix) <- merged_test$patient_id

  cat("  Testing:", length(snp_cols_test), "variants in test set...\n")

  res_test <- future_lapply(
    seq_along(snp_cols_test), run_cox,
    snps = snp_cols_test, geno_mat = geno_test_matrix, merged_df = merged_test
  )
  res_test <- rbindlist(res_test)

  fwrite(res_test, sprintf("results/resampling/test_iter_%s.csv", iter_label))

  # --- Directional concordance + one-tailed test ------------------------------
  train_sig <- sig_snps[, .(SNP, Beta_train = Beta, HR_train = HR, P_train = P)]
  dt <- merge(res_test, train_sig, by = "SNP", all = FALSE)

  dt[, direction_match := sign(Beta) == sign(Beta_train)]
  dt[, p_one_tailed := fifelse(
    direction_match & !is.na(Z),
    pnorm(-abs(Z)),
    NA_real_
  )]

  replicated <- dt[direction_match == TRUE &
                   !is.na(p_one_tailed) &
                   p_one_tailed < val_threshold]

  cat("  Replicated:", nrow(replicated), "variants\n")

  if (nrow(replicated) > 0) {
    replicated[, iteration := iter]
    return(replicated)
  }

  return(NULL)
}

# ==============================================================================
# AGGREGATE: combine results from all 100 iterations
# ==============================================================================
aggregate_results <- function() {

  library(ggplot2)

  cat("\n========================================\n")
  cat("AGGREGATING RESULTS\n")
  cat("========================================\n")

  all_replicated <- list()

  for (i in 1:n_iters) {
    iter_label <- sprintf("%03d", i)
    train_file <- sprintf("results/resampling/train_iter_%s.csv", iter_label)
    test_file  <- sprintf("results/resampling/test_iter_%s.csv", iter_label)

    if (!file.exists(train_file) || !file.exists(test_file)) next

    train_dt <- fread(train_file)
    test_dt  <- fread(test_file)

    train_sig <- train_dt[!is.na(P) & P < sig_threshold, .(
      SNP, Beta_train = Beta, HR_train = HR, P_train = P
    )]

    if (nrow(train_sig) == 0) next

    dt <- merge(test_dt, train_sig, by = "SNP", all = FALSE)
    if (nrow(dt) == 0) next

    dt[, direction_match := sign(Beta) == sign(Beta_train)]
    dt[, p_one_tailed := fifelse(
      direction_match & !is.na(Z),
      pnorm(-abs(Z)),
      NA_real_
    )]

    replicated <- dt[direction_match == TRUE &
                     !is.na(p_one_tailed) &
                     p_one_tailed < val_threshold,
                     .(SNP, Beta, HR, SE, Z, P,
                       Beta_train, HR_train, P_train,
                       p_one_tailed)]

    if (nrow(replicated) > 0) {
      replicated[, iteration := i]
      all_replicated[[length(all_replicated) + 1]] <- replicated
    }
  }

  all_replicated <- rbindlist(all_replicated, use.names = TRUE, fill = TRUE)

  cat("Total replicated signals:", nrow(all_replicated), "\n")
  cat("Unique variants:         ", uniqueN(all_replicated$SNP), "\n")

  # --- Per-variant summary ----------------------------------------------------
  resampling_summary <- all_replicated[, .(
    n_replicated = .N,
    beta_mean    = mean(Beta),
    beta_sd      = sd(Beta),
    beta_min     = min(Beta),
    beta_max     = max(Beta),
    hr_mean      = mean(HR),
    p_median     = median(p_one_tailed),
    sign_flip    = length(unique(sign(Beta))) > 1
  ), by = SNP]

  setorder(resampling_summary, -n_replicated)

  cat("\nTop replicated variants:\n")
  print(head(resampling_summary, 20))

  # --- Frequency plot ---------------------------------------------------------
  plot_df <- resampling_summary[n_replicated >= 10]

  if (nrow(plot_df) > 0) {
    plot_df[, SNP := factor(SNP, levels = SNP[order(n_replicated)])]

    p <- ggplot(plot_df, aes(x = n_replicated, y = SNP)) +
      geom_bar(stat = "identity", fill = "#2C7FB8") +
      labs(
        title = expression("Variants replicated in" >= "10 of 100 iterations"),
        x = "Number of iterations replicated (concordant HR, one-tailed P < 0.05)",
        y = "Variant"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 11)
      )

    ggsave("results/resampling_frequency_plot.png",
           plot = p, width = 8, height = 6, dpi = 300)
    cat("\nPlot saved: results/resampling_frequency_plot.png\n")
  } else {
    cat("\nNo variants replicated in >= 10 iterations; skipping plot.\n")
  }

  # --- Save -------------------------------------------------------------------
  fwrite(all_replicated, "results/resampling_all_replicated.csv",
         sep = ",", quote = FALSE, col.names = TRUE)
  fwrite(resampling_summary, "results/resampling_summary.csv",
         sep = ",", quote = FALSE, col.names = TRUE)

  cat("\nDone. Results saved to results/\n")
}

# ==============================================================================
# MAIN: determine run mode
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # SLURM array mode
  slurm_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
  if (!is.na(slurm_id)) {
    plan(multisession, workers = 10)
    run_one_iteration(as.integer(slurm_id))
  } else {
    stop("No argument provided. Use: 'all', 'aggregate', or an iteration number (1-100).")
  }

} else if (args[1] == "aggregate") {
  aggregate_results()

} else if (args[1] == "all") {
  plan(multisession, workers = 10)
  all_results <- list()
  for (iter in 1:n_iters) {
    result <- run_one_iteration(iter)
    if (!is.null(result) && nrow(result) > 0) {
      all_results[[length(all_results) + 1]] <- result
    }
  }
  aggregate_results()

} else {
  plan(multisession, workers = 10)
  iter <- as.integer(args[1])
  if (is.na(iter) || iter < 1 || iter > n_iters) {
    stop("Invalid argument. Use: 'all', 'aggregate', or a number 1-100.")
  }
  run_one_iteration(iter)
}
