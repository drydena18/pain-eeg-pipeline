# =============================================================================
# run_robustness.R
# -----------------------------------------------------------------------------
# Four robustness tests, all evaluated on the held-out set:
#
#   1. PERMUTATION TESTING
#      Shuffles held-out pain ratings, evaluates prediction quality on the
#      shuffled labels, repeats n_perm times → null distribution. Compares
#      real held-out performance (from cv_external) against null. A
#      permutation p-value < 0.05 means the model predicts pain ratings
#      better than chance. The model is NOT refit — only the evaluation
#      labels are shuffled.
#
#   2. BOOTSTRAP CONFIDENCE INTERVALS
#      Resamples held-out subjects with replacement, evaluates all five
#      metrics each time → 95% CI on held-out performance. Supplements the
#      point estimates in cv_external with uncertainty bounds. CI width
#      also flags whether the holdout N is too small for stable estimates.
#
#   3. SENSITIVITY ANALYSES
#      Drops one experiment at a time from the TRAINING set, refits bam()
#      on the remaining training subjects, evaluates on the held-out set
#      (which stays fixed). If dropping one experiment causes large
#      performance swings, that experiment is driving the result. Also
#      re-evaluates at alternative FOOOF r² QC thresholds (0.75, 0.85)
#      to check whether the QC floor meaningfully changes conclusions.
#
#   4. SPLIT-HALF STABILITY
#      Randomly splits the TRAINING set in half (stratified by cap_size),
#      fits bam() independently on each half, generates predictions on
#      held-out subjects from both half-models, compares the predictions
#      (not just their performance metrics) via correlation. High
#      prediction-to-prediction correlation means the model's functional
#      form is stable regardless of which subjects happened to be in
#      training.
#
# All tests write to cv_robustness table in results.duckdb.
# Permutation and bootstrap are the most time-consuming — adjust
# n_perm and n_boot in robustness_cfg to trade speed for precision.
# =============================================================================

library(dplyr)
library(purrr)
library(mgcv)
library(stringr)
library(tibble)

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("cv/cv_helpers.R")
source("cv/cv_split.R")

# =============================================================================
# CONFIG
# =============================================================================
robustness_cfg <- list(
  p_sig          = 0.05,
  n_perm         = 1000L,   # permutation test iterations
  n_boot         = 500L,    # bootstrap resamples
  n_split_half   = 50L,     # split-half repetitions (each uses a new random split)
  fooof_thresholds = c(0.75, 0.85),  # alternative QC thresholds for sensitivity
  seed           = 42L
)

# =============================================================================
# HELPERS
# =============================================================================

# Refit a bam() from a reference model on a new data subset
refit_bam <- function(mod_ref, df_new, label = "") {
  tryCatch(
    bam(formula(mod_ref), data = df_new,
        method = "fREML", select = TRUE, discrete = TRUE, nthreads = 4L),
    error = function(e) {
      if (nchar(label) > 0)
        message("  [WARN] refit failed (", label, "): ", conditionMessage(e))
      NULL
    }
  )
}

# Stratified-by-cap_size half split of training subjects
split_train_half <- function(df_train, seed_offset) {
  set.seed(robustness_cfg$seed + seed_offset)
  subj_df <- df_train %>%
    distinct(global_subjid, cap_size) %>%
    group_by(cap_size) %>%
    mutate(n_grp = n()) %>%
    mutate(
      in_first_half = sample(c(rep(TRUE,  ceiling(n_grp[1] / 2)),
                                rep(FALSE, floor(n_grp[1]   / 2))))
    ) %>%
    ungroup()
  list(
    half1 = df_train %>%
      filter(global_subjid %in%
               subj_df$global_subjid[subj_df$in_first_half]),
    half2 = df_train %>%
      filter(global_subjid %in%
               subj_df$global_subjid[!subj_df$in_first_half])
  )
}

# =============================================================================
# CONNECT
# =============================================================================
data_con    <- connect_data_db(gamm_cfg$data_db_path, read_only = TRUE)
results_con <- connect_results_db(gamm_cfg$results_db_path)
init_cv_table(results_con, "cv_robustness")

models <- select_cv_models(results_con, p_sig = robustness_cfg$p_sig)
if (nrow(models) == 0) stop("No significant models — run source_core.R first.")
message("Robustness tests: ", nrow(models), " model(s).\n")

# =============================================================================
# PER-MODEL LOOP
# =============================================================================
for (i in seq_len(nrow(models))) {
  mod_row     <- models[i, ]
  this_metric <- mod_row$metrics
  group_value <- mod_row$group_value
  level       <- mod_row$level
  mod_bam     <- readRDS(mod_row$rds_path)

  message(sprintf("[%d/%d] Robustness: %s @ %s (%s)",
                  i, nrow(models), this_metric, group_value, level))

  # Load data
  df_train   <- load_cv_data(data_con, level, group_value,
                              this_metric, split_group = "train")
  df_holdout <- load_holdout_data(data_con, level, group_value, this_metric)

  if (is.null(df_holdout) || nrow(df_holdout) == 0) {
    message("  Skipping — no holdout data.")
    next
  }

  df_holdout   <- predict_holdout(mod_bam, df_holdout)
  real_metrics <- eval_metrics(df_holdout$pain_rating, df_holdout$predicted)
  n_ho         <- n_distinct(df_holdout$global_subjid)

  # ===========================================================================
  # TEST 1: PERMUTATION TESTING
  # ---------------------------------------------------------------------------
  message("  [1/4] Permutation test (", robustness_cfg$n_perm, " iterations)...")
  set.seed(robustness_cfg$seed)

  perm_r <- vapply(seq_len(robustness_cfg$n_perm), function(p) {
    shuffled_actual <- sample(df_holdout$pain_rating)
    eval_metrics(shuffled_actual, df_holdout$predicted)$pearson_r
  }, numeric(1))

  perm_p_value <- mean(abs(perm_r) >= abs(real_metrics$pearson_r), na.rm = TRUE)
  perm_p_value <- max(perm_p_value, 1 / robustness_cfg$n_perm)  # floor at 1/n_perm

  message(sprintf("    Real r=%.3f  Null mean=%.3f  Null SD=%.3f  perm_p=%.4f",
                  real_metrics$pearson_r, mean(perm_r), sd(perm_r), perm_p_value))

  perm_row <- make_cv_row(
    "permutation", level, group_value, this_metric,
    "permutation_test", n_ho, real_metrics
  ) %>%
    mutate(fold_id = sprintf("perm_p=%.4f_null_r_mean=%.3f_null_r_sd=%.3f",
                              perm_p_value, mean(perm_r), sd(perm_r)))
  write_cv_rows(results_con, "cv_robustness", perm_row)

  # ===========================================================================
  # TEST 2: BOOTSTRAP CONFIDENCE INTERVALS
  # ---------------------------------------------------------------------------
  message("  [2/4] Bootstrap CIs (", robustness_cfg$n_boot, " resamples)...")
  set.seed(robustness_cfg$seed + 1L)

  holdout_subjects <- unique(df_holdout$global_subjid)

  boot_metrics <- map_dfr(seq_len(robustness_cfg$n_boot), function(b) {
    boot_subjs <- sample(holdout_subjects, replace = TRUE)
    boot_df    <- map_dfr(boot_subjs, ~filter(df_holdout, global_subjid == .x))
    m <- eval_metrics(boot_df$pain_rating, boot_df$predicted)
    tibble(pearson_r = m$pearson_r, rmse = m$rmse, r_sq = m$r_sq,
           mae = m$mae, spearman_rho = m$spearman_rho)
  })

  ci_width_r <- diff(quantile(boot_metrics$pearson_r, c(0.025, 0.975), na.rm = TRUE))
  message(sprintf("    r CI: [%.3f, %.3f]  width=%.3f  %s",
                  quantile(boot_metrics$pearson_r, 0.025, na.rm = TRUE),
                  quantile(boot_metrics$pearson_r, 0.975, na.rm = TRUE),
                  ci_width_r,
                  if (ci_width_r > 0.20) "[FLAG: wide CI]" else "OK"))

  boot_row <- make_cv_row(
    "bootstrap_ci", level, group_value, this_metric,
    "bootstrap", n_ho,
    list(pearson_r    = median(boot_metrics$pearson_r, na.rm = TRUE),
         rmse         = median(boot_metrics$rmse, na.rm = TRUE),
         r_sq         = median(boot_metrics$r_sq, na.rm = TRUE),
         mae          = median(boot_metrics$mae, na.rm = TRUE),
         spearman_rho = median(boot_metrics$spearman_rho, na.rm = TRUE),
         n_obs        = real_metrics$n_obs)
  ) %>%
    mutate(fold_id = sprintf(
      "r_ci=[%.3f,%.3f]_rmse_ci=[%.3f,%.3f]",
      quantile(boot_metrics$pearson_r, 0.025, na.rm = TRUE),
      quantile(boot_metrics$pearson_r, 0.975, na.rm = TRUE),
      quantile(boot_metrics$rmse, 0.025, na.rm = TRUE),
      quantile(boot_metrics$rmse, 0.975, na.rm = TRUE)
    ))
  write_cv_rows(results_con, "cv_robustness", boot_row)

  # ===========================================================================
  # TEST 3: SENSITIVITY ANALYSES
  # ---------------------------------------------------------------------------
  message("  [3/4] Sensitivity analyses...")

  if (!is.null(df_train)) {
    experiments <- unique(as.character(df_train$experiment_id))

    # Drop-one-experiment
    for (exp_drop in experiments) {
      df_loo_exp <- filter(df_train,
                            as.character(experiment_id) != exp_drop)
      if (n_distinct(df_loo_exp$global_subjid) < 5) next

      mod_loo <- refit_bam(mod_bam, df_loo_exp,
                            label = paste0("drop_exp_", exp_drop))
      if (is.null(mod_loo)) next

      df_ho_pred  <- predict_holdout(mod_loo, df_holdout)
      m_loo       <- eval_metrics(df_ho_pred$pain_rating, df_ho_pred$predicted)
      delta_r     <- m_loo$pearson_r - real_metrics$pearson_r

      sens_row <- make_cv_row(
        "sensitivity_drop_experiment", level, group_value, this_metric,
        paste0("drop_exp_", exp_drop),
        n_distinct(df_loo_exp$global_subjid), m_loo
      ) %>%
        mutate(fold_id = sprintf("drop_exp_%s_delta_r=%.3f", exp_drop, delta_r))
      write_cv_rows(results_con, "cv_robustness", sens_row)

      message(sprintf("    drop exp %-4s: r=%.3f  Δr=%.3f",
                      exp_drop, m_loo$pearson_r, delta_r))
    }

    # Alternative FOOOF QC thresholds
    for (threshold in robustness_cfg$fooof_thresholds) {
      df_train_alt <- load_cv_data(data_con, level, group_value,
                                    this_metric, split_group = "train")
      if (!is.null(df_train_alt) && "fooof_r2" %in% names(df_train_alt)) {
        df_train_alt <- filter(df_train_alt,
                                is.na(fooof_r2) | fooof_r2 >= threshold)
        if (n_distinct(df_train_alt$global_subjid) < 5) next

        mod_alt <- refit_bam(mod_bam, df_train_alt,
                              label = paste0("fooof_r2_", threshold))
        if (is.null(mod_alt)) next

        df_ho_alt <- predict_holdout(mod_alt, df_holdout)
        m_alt     <- eval_metrics(df_ho_alt$pain_rating, df_ho_alt$predicted)

        alt_row <- make_cv_row(
          "sensitivity_fooof_threshold", level, group_value, this_metric,
          paste0("fooof_r2_", threshold),
          n_distinct(df_train_alt$global_subjid), m_alt
        )
        write_cv_rows(results_con, "cv_robustness", alt_row)
        message(sprintf("    FOOOF r²≥%.2f: r=%.3f  n_train=%d",
                        threshold, m_alt$pearson_r,
                        n_distinct(df_train_alt$global_subjid)))
      }
    }
  } else {
    message("  Skipping sensitivity — no training data available.")
  }

  # ===========================================================================
  # TEST 4: SPLIT-HALF STABILITY
  # ---------------------------------------------------------------------------
  message("  [4/4] Split-half stability (",
          robustness_cfg$n_split_half, " repetitions)...")

  if (!is.null(df_train)) {
    split_half_rs <- vapply(seq_len(robustness_cfg$n_split_half), function(rep) {
      halves <- split_train_half(df_train, seed_offset = rep)

      if (n_distinct(halves$half1$global_subjid) < 5 ||
          n_distinct(halves$half2$global_subjid) < 5) return(NA_real_)

      mod_h1 <- refit_bam(mod_bam, halves$half1)
      mod_h2 <- refit_bam(mod_bam, halves$half2)
      if (is.null(mod_h1) || is.null(mod_h2)) return(NA_real_)

      pred_h1 <- predict_holdout(mod_h1, df_holdout)$predicted
      pred_h2 <- predict_holdout(mod_h2, df_holdout)$predicted

      # Prediction-to-prediction correlation (not performance — stability)
      if (all(is.na(pred_h1)) || all(is.na(pred_h2))) return(NA_real_)
      cor(pred_h1, pred_h2, use = "complete.obs")
    }, numeric(1))

    sh_mean <- mean(split_half_rs, na.rm = TRUE)
    sh_sd   <- sd(split_half_rs, na.rm = TRUE)
    message(sprintf("    Prediction correlation: mean=%.3f  SD=%.3f  [%.3f, %.3f]",
                    sh_mean, sh_sd,
                    quantile(split_half_rs, 0.025, na.rm = TRUE),
                    quantile(split_half_rs, 0.975, na.rm = TRUE)))

    sh_row <- make_cv_row(
      "split_half_stability", level, group_value, this_metric,
      "split_half", n_ho,
      list(pearson_r = sh_mean, rmse = NA_real_, r_sq = sh_mean^2,
           mae = NA_real_, spearman_rho = NA_real_,
           n_obs = real_metrics$n_obs)
    ) %>%
      mutate(fold_id = sprintf(
        "pred_r_mean=%.3f_sd=%.3f_ci=[%.3f,%.3f]",
        sh_mean, sh_sd,
        quantile(split_half_rs, 0.025, na.rm = TRUE),
        quantile(split_half_rs, 0.975, na.rm = TRUE)
      ))
    write_cv_rows(results_con, "cv_robustness", sh_row)
  } else {
    message("  Skipping split-half — no training data available.")
  }

  message("  Done: ", this_metric, " @ ", group_value)
}

disconnect_db(data_con)
close_results_db(results_con)
message("\nRobustness tests complete. Results: cv_robustness table in ",
        gamm_cfg$results_db_path)