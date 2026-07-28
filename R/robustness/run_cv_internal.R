# =============================================================================
# run_cv_internal.R
# ------------------------------------------------------------------------------
# Leave-One-Subject-Out (LOSO) cross-validation within the TRAINING set.
# For every significant model in gamm_fitted (source + channel), leaves
# one training subject out, refits bam() on the remaining training
# subjects, predicts on the left-out subject, evaluates.
#
# The held-out set is NEVER touched here. This is purely an internal
# generalization check within the training distribution, run BEFORE
# looking at held-out data. It answers: "does the model generalize to
# unseen subjects from the same distribution as training?"
#
# Writes to results.duckdb: cv_internal table (one row per fold +
# one summary row per metric x ROI).
# =============================================================================

library(dplyr)
library(purrr)
library(mgcv)
library(stringer)

source("gamm_defaults.R")
source("helpers/db_helpers.R")
source("cv/cv_helpers.R")
source("cv/cv_split.R")

# =============================================================================
# CONFIG
# =============================================================================
cv_int_cfg <- list(
    p_sig = 0.05,
    min_fold_obs = 5L # skip a fold if the left-out subject has < 5 folds
)

# =============================================================================
# CONNECT
# =============================================================================
data_con <- connect_data_db(gamm_cfg$data_db_path, read_only = TRUE)
if (nrow(models) == 0) stop("No significant models found - run source_core.R first.")

message("LOSO CV: ", nrow(models), " model(s) to evaluate.")
message(" (Training set only - held-out subjects not touched.)\n")

# =============================================================================
# PER-MODEL LOOP
# =============================================================================
for (i in seq_len(nrow(models))) {
    mod_row <- models[i, ]
    this_metric <- mod_row$metrics
    group_value <- mod_row$group_value
    level <- mod_row$level
    mod_bam <- readRDS(mod_row$rds_path)

    message(sprintf("[%d/%d] LOSO: %s @ %s (%s) r_sq = %.3f", i, nrow(models), this_metric, group_value, level, mod_row$r_sq))

    # Load full training data for this metric x group
    df_train <- load_cv_data(data_con, level, group_value, this_metric, split_group = "train")
    if (is.null(df_train) || nrow(df_train) == 0) {
        message(" Skipping - no training data available.")
        next
    }

    n_subjects <- n_distinct(df_train$global_subjid)
    message("Training subjects: ", n_subjects, "Trials: ", nrow(df_train))

    # Reconstrict the base formula from the fitted model
    # (avoids re-specifying it manually - takes it directly from the bam object)
    original_formula <- formula(mod_bam)

    # Generate LOSO folds
    folds <- generate_loso_folds(df_train)
    fold_rows <- list()

    for (subj_id in names(fold)) {
        fold <- folds[[subj_id]]
        df_fold_train <- fold$train
        df_fold_test <- fold$test

        if (nrow(df_fold_test) < cv_int_cfg$min_fold_obs) next

        # Refit on training-minus-one - full bam(), not an approximation
        mod_fold <- tryCatch(
            bam(original_formula, data = df_fold_train, method = "fREML", select = TRUE, discrete = TRUE, nthreads = 4L),
            error = function(e) {
                message("[WARN] refit failed (sub ", subj_id, "): ", conditionMessage(e))
                NULL
            }
        )

        if (is.null(mod_fold)) next

        # Predict on left-our subject
        df_fold_test <- predict_holdout(mod_fold, df_fold_test)
        m <- eval_metrics(df_fold_test$pain_rating, df_fold_test$predicted)

        fold_rows[[subj_id]] <- make_cv_row(
            cv_type = "loso",
            level = level,
            group_value = group_value,
            metrics = this_metric,
            fold_id = subj_id,
            n_subjects = 1L,
            metrics_list = m
        )
    }

    if (length(fold_rows) == 0) next

    all_folds_df <- bind_rows(fold_rows)

    # Write per-fold rows
    write_cv_rows(results_con, "cv_internal", all_folds_df)

    # Write summary row (mean across folds)
    fold_summary <- summaize_fold_results(all_folds_df)
    summary_row <- tibble(
        fitted_at = as.character(Sys.time()),
        cv_type = "loso_summary",
        level = level,
        group_value = group_value,
        metrics = this_metric,
        fold_id = "summary",
        n_subjects = as.integer(n_subjects),
        n_trials = as.integer(sum(all_folds_df$n_trials, na.rm = TRUE)),
        pearson_r = fold_summary$pearson_r_mean,
        rmse = fold_summary$rmse_mean,
        r_sq = fold_summary$r_sq_mean,
        mae = fold_summary$mae_mean,
        spearman_rho = fold_summary$spearman_rho_mean,
        n_obs = as.integer(sum(all_folds_df$n_obs, na.rm = TRUE))
    )
    write_cv_rows(results_con, "cv_internal", summary_row)

    message(sprintf(" LOSO summary: r = %.3f RMSE = %.3f R^2 = %.3f MAR = %.3f p = %.3f (%d/%d folds)",
                    fold_summary$pearson_r_mean, fold_summary$rmse_mean,
                    fold_summary$r_sq_mean, fold_summary$mae_mean,
                    fold_summary$spearman_rho_mean,
                    fold_summary$n_folds_complete, fold_summary$n_folds))
}

disconnect_db(data_con)
close_results_db(results_con)
message("\nLOSO CV complete. Results: cv_internal table in ", gamm_cfg$results_db_path)