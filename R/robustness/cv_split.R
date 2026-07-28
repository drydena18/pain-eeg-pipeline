# =============================================================================
# cv_split.R
# -----------------------------------------------------------------------------
# Manages data partitioning for both CV strategies:
#
#   LOSO (internal):    generates leave-one-subject-out fold definitions
#                       from the training set. Each fold is identified by
#                       the subject being left out. Called by
#                       run_cv_internal.R.
#
#   Holdout (external): loads held-out subjects for a given level/group.
#                       A thin wrapper around load_holdout_source_metrics_wide()
#                       and load_holdout_channel_metrics_wide() in db_helpers.R,
#                       with the metric z-scoring and QC filtering applied
#                       consistently. Called by run_cv_external.R.
#
# Neither function touches the split assignment itself - that is owned by
# merge_core.R / split_helpers.R. These functions are read-only consumers
# of subject_split.
# =============================================================================

library(dplyr)
library(purrr)

# =============================================================================
# generate_loso_folds()
# -----------------------------------------------------------------------------
# Returns a named list, one entry per unique subject in df_train. Each
# entry is a list(train = df, test = df) where test is that subject's
# rows and train is everyone else. This is the correct LOSO split: the
# held-out subject contributes NO rows to model fitting for that fold.
#
# NOTE: subject identity is determined by global_subjid, not row position.
# This ensures that multi-session subjects (where one global_subjid may
# have many rows) are left out entirely, not partially.
#
# Called once per metric x ROI combination in run_cv_internal.R, since
# the fold structure is the same regardless of which metric is being
# tested.
# =============================================================================
generate_loso_folds <- function(df_train) {
    subjects <- unique(as.character(df_train$global_subjid))

    map(subjects, function(s) {
        list(
            subject_id = s,
            train = filter(df_train, as.character(global_subjid) != s),
            test = filter(df_train, as.character(global_subjid) == s)
        )
    }) %>%
    setNames(subjects)
}

# =============================================================================
# load_heldout_data()
# -----------------------------------------------------------------------------
# Loads held-out data for one level/group/metric combination.
# A thin wrapper that applies QC filtering and z-scoring consistently,
# matching exactly what load_cv_data() does for training data - so
# evaluation metrics are computed on comparably-prepared data.
#
# Returns NULL if the metric is not present in the pivoted holdout data
# (e.g. delta_erd is absent for some subjects - handled gracefully).
# =============================================================================
load_holdout_data <- function(data_con, level, group_value, this_metric) {
    load_cv_data(data_con, level, group_value, this_metric, split_group = "holdout")
}

# =============================================================================
# summarize_fold_results()
# -----------------------------------------------------------------------------
# Aggregates per-fold eval_metrics results across all LOSO folds.
# Returns mean +- SD for each metric across folds, plus n_folds and
# the number of folds where prediction was successful (non-NA).
# Called at the end of run_cv_internal.R to produce a fold-level summary
# row alongside the per-fold rows already written to cv_internal.
# =============================================================================
summarize_fold_results <- function(fold_rows) {
    if (is.null(fold_rows) || nrow(fold_rows) == 0) return(NULL)

    fold_rows %>%
        summarise(
            n_folds = n(),
            n_folds_complete = sum(!is.na(pearson_r)),
            pearson_r_mean = mean(pearson_r, na.rm = TRUE),
            pearson_r_sd = sd(pearson_r, na.rm = TRUE),
            rmse_mean = mean(rmse, na.rm = TRUE),
            rmse_sd = sd(rmse, na.rm = TRUE),
            r_sq_mean = mean(r_sq, na.rm = TRUE),
            r_sq_sd = sd(r_sq, na.rm = TRUE),
            mae_mean = mean(mae, na.rm = TRUE),
            mae_sd = sd(mae, na.rm = TRUE),
            spearman_rho_mean = mean(spearman_rho, na.rm = TRUE),
            spearman_rho_sd = sd(spearman_rho, na.rm = TRUE)
        )
}