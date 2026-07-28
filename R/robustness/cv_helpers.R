# =============================================================================
# cv_helpers.R
# ------------------------------------------------------------------------------
# Shared utilities for all CV and robustness scripts. No standalone logic -
# sourced by cv_split.R, run_cv_internal.R, run_cv_external.R, and
# run_robustness.R.
#
# Contents:
#   - eval_metrics()        five prediction quality metrics
#   - predict_holdout()     bam() predictions on new subjects
#   - select_cv_models()    pulls significant hits from gamm_fitted
#   - load_cv_data()        loads wide data for one level/group
#   - init_cv_tables()      creates a CV results table in results.duckdb
#   - write_cv_rows()       upserts rows into a CV results table
# =============================================================================

library(DBI)
library(duckdb)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(mgcv)

# =============================================================================
# eval_metrics()
# -----------------------------------------------------------------------------
# Computes all five prediction quality metrics between actual and predicted
# pain ratings. Returns a named list; all metrics are NA if fewer than 3
# complete pairs exist (too few to be meaningful).
# =============================================================================
eval_metrics <- function(actual, predicted) {
    idx <- !is.na(actual) & !is.na(predicted)
    a <- actual[idx]
    p <- predicted[idx]
    n <- length(a)

    if (n < 3) {
        return(list(pearson_r = NA_real_, rmse = NA_real_, r_sq = NA_real_, mae = NA_real_, spearman_rho = NA_real_, n_obs = n))
    }

    pearson_r <- cor(a, p, method = "pearson")
    spearman_rho <- cor(a, p, method = "spearman")
    rmse <- sqrt(mean((a - p)^2))
    # R^2 as squared correlation (variance in actual explained by linear
    # relationship with predicted - standard for mixed model evaluation)
    r_sq <- pearson_r^2
    
    list(pearson_r = pearson_r, rmse = rmse, r_sq = r_sq, mae = mae, spearman_rhp = spearman_rho, n_obs = n)
}

# =============================================================================
# predict_holdout()
# -----------------------------------------------------------------------------
# Generates predictions from a fitted bam() object on new subjects.
# Uses exclude = "s(global_subjid)" so held-out subjects (not in the
# training RE levels) get population-level predictions rather than NA or
# an error. Returns the data frame with a 'predicted' column appended.
# =============================================================================
predict_holdout <- function(mod, df_new) {
    re_terms <- names(mod$smooth)[
        vapply(mod$smooth, function(s) inherits(s, "random.effect"), logical(1))
    ]
    df_new$predicted <- tryCatch(
        as.numeric(predict(mod, newdata = df_new, exclude = re_terms, type = "response")),
        error = function(e) {
            message("[WARN] predict() failed: ", conditionMessage(e))
            rep(NA_real_, nrow(df_new))
        }
    )
    df_new
}

# =============================================================================
# select_cv_models()
# -----------------------------------------------------------------------------
# Pulls all significant single-metric models from gamm_fitted, for both
# 'source' and 'channel' levels. Mirrors the selection in power_analysis.R
# so CV runs on exactly the same model set.
# Returns a tibble with: level, group_value, model_name, metrics, rds_path, r_sq, p_value_num
# =============================================================================
select_cv_models <- function(results_con, p_sig = 0.05) {
    dbGetQuery(results_con, "
        SELECT level, grouo_value, model_name, metrics, n_metrics, p_value, r_sq, rds_path
        FROM gamm_fitted
        WHERE n_metrics = 1
        ") %>%
        mutate(
            p_value_num = suppressWarnings(as.numeric(p_value)),
            r_sq = suppressWarnings(as.numeric(r_sq))
        ) %>%
        filter(
            !is.na(p_value_num), p_value_num < p_sig,
            !is.na(r_sq),
            !is.na(rds_path), file.exists(rds_path)
        ) %>%
        arrange(level, group_value, r_sq)
}

# =============================================================================
# load_cv_data()
# -----------------------------------------------------------------------------
# Loads wide data for one level ('source' or 'channel') and group (ROI or
# 'all_channels'), filtered to the requested split_group. Z-scores the
# specified metric column and base covariates.
# Returns NULL if the metric is not found in the pivoted data.
# =============================================================================
load_cv_data <- function(data_con, level, group_value, this_metric, split_group = "train") {
    raw_metric <- str_remove(this_metric, "_z$")

    df <- if (level == "source") {
    if (split_group == "train") {
      load_source_metrics_wide(data_con, group_value)
    } else {
      load_holdout_source_metrics_wide(data_con, group_value)
    }
  } else {
    if (split_group == "train") {
      load_channel_metrics_wide(data_con, group_value)
    } else {
      load_holdout_channel_metrics_wide(data_con, group_value)
    }
  }

    if (is.null(df) || nrow(df) == 0) return(NULL)
    if (!raw_metric %in% names(df)) return(NULL)

    df %>%
        filter(!is.na(pain_rating), !is.na(laser_power), !is.na(trial_index), !is.na(age), !is.na(global_subjid)) %>%
        mutate(
            global_subjid = as.factor(as.character(global_subjid)),
            experiment_id = as.factor(as.character(experiment_id)),
            sex = as.factor(as.character(sex)),
            cap_size = as.factor(as.character(cap_size)),
            age_z = as.numeric(scale(age)),
            laser_power_z = as.numeric(scale(laser_power)),
            trail_index_z = as.numeric(scale(trial_index)),
            !!this_metric := as.numeric(scale(.data[[raw_metric]]))
        )
}

# =============================================================================
# DuckDB table init + write
# =============================================================================

# Schema shared by all three CV results tables = distinguished only by
# table name (cv_internal, cv_external, cv_robustness) and cv_type column.
CV_TABLE_SCHEMA <- "
    fitted_at VARCHAR,
    cv_type VARCHAR,
    level VARCHAR,
    group_value VARCHAR,
    metrics VARCHAR,
    fold_id VARCHAR,
    n_subjects INTEGER,
    n_trials INTEGER,
    person_r DOUBLE,
    rmse DOUBLE,
    r_sq DOUBLE,
    mae DOUBLE,
    spearman_rho DOUBLE,
    n_obs INTEGER
"

init_cv_table <- function(con, table_name) {
    doExecute(con, sprintf(
        "CREATE TABLE IF NOT EXISTS %s (%s)", table_name, CV_TABLE_SCHEMA
    ))
}

write_cv_rows <- function(con, table_name, df_rows) {
    if (is.null(rows) || nrow(rows) == 0) return(invisible(NULL))
    doWriteTable(con, table_name, rows, append = TRUE)
    invisible(NULL)
}

# Convenience: build one result row tibble for eval_metrics output
make_cv_row <- function(cv_type, level, group_value, metrics, fold_id, n_subjects, metrics_list) {
    tibble(
        fitted_at = as.character(Sys.time()),
        cv_type = cv_type,
        level = level,
        group_value = group_value,
        metrics = metrics,
        fold_id = as.character(fold_id),
        n_subjects = as.integer(n_subjects),
        n_trials = as.integer(metrics_list$n_obs),
        pearson_r = metrics_list$pearson_r,
        rmse = metrics_list$rmse,
        r_sq = metrics_list$r_sq,
        mae = metrics_list$mae,
        spearman_rho = metrics_list$spearman_rho,
        n_obs = as.integer(metrics_list$n_obs)
    )
}