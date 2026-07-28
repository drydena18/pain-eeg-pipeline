# =============================================================================
# run_cv_external.R
# -----------------------------------------------------------------------------
# External validation on the TRUE held-out set (subject_split = 'holdout').
# Takes the final fitted bam() models from source_core.R and channel_core.R
# (trained on all training subjects, never modified here) and evaluates
# their predictions on held-out subjects.
#
# This is the PRIMARY generalization claim for the thesis. It answers:
# "do these alpha metrics predict pain ratings in completely unseen
# subjects, drawn from the same population but never seen during fitting?"
#
# The model is NOT refit here. This is evaluation-only. The training
# boundary is enforced by load_holdout_source_metrics_wide() /
# load_holdout_channel_metrics_wide() in db_helpers.R, which filter to
# split_group = 'holdout' at the SQL level.
#
# Writes to results.duckdb: cv_external table (one row per subject +
# one aggregate row per metrics x ROI).
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
cv_ext_cfg <- list(
    p_sig = 0.05
)

# =============================================================================
# CONNECT
# =============================================================================
data_con <- connect_data_db(gamm_cfg$data_db_path, read_only = TRUE)
results_con <- connect_results_db(gamm_cfg$results_db_path)
init_cv_table(results_con, "cv_external")

# Check holdout subjects exist
holdout_n <- dbGetQuery(data_con, "
    SELECT COUNT(*) AS n FROM subject_split WHERE split_group = 'holdout'"
)$n
if (holdout_n == 0) stop("No holdout subjects found in subject_split.")
message("Holdout subjects available: ", holdout_n)

# =============================================================================
# SELECT MODELS
# =============================================================================
models <- select_cv_models(results_con, p_sig = cv_ext_cfg$p_sig)
if (nrow(models) == 0) stop("No significant models - run source_core.R first.")
message("External CV: ", nrow(models), " model(s).\n")

# =============================================================================
# PER-MODEL LOOP
# =============================================================================
for (i in seq_len(nrow(models))) {
    mod_row <- models[i, ]
    this_metric <- mod_row$metrics
    group_value <- mod_row$group_value
    level <- mod_row$level
    mod_bam <- readRDS(mod_row$rds_path)

    message(sprintf("[%d/%d] External: %s @ %s (%s) r_sq = %.3f",
                    i, nrow(models), this_metric, group+value, level, mod_row$r_sq))

    # Load held-out data
    df_holdout <- load_holdout_data(data_con, level, group_value, this_metric)
    if (is.null(df_holdout) || nrow(df_holdout) == 0) {
        message("Skipping - no holdout data available.")
        next
    }

    n_ho_subj <- n_distinct(df_holdout$global_subjid)
    cap_ho <- df_holdout %>%
        distinct(global_subjid, cap_size) %>% count(cap_size)
    message("Holdout subjects: ", n_ho_subj, " Trials: ", nrow(df_holdout))
    for (j in seq_len(nrow(cap_ho))) {
        message(" cap_size ", cap_ho$cap_size[j], ": ", cap_ho$n[j])
    }

    # Predict - model is NOT refit, just evaluated
    df_holdout <- predict_holdout(mod_bam, df_holdout)

    # Per-subject rows
    subj_rows <- df_holdout %>%
        group_by(global_subjid, cap_size) %>%
        group_map(function(df_s, key) {
            m <- eval_metrics(df_s$pain_rating, df_s$predicted)
            make_cv_row("expernal_per_subject", level, group_value, this_metric, key$global_subjid, 1L, m) %>%
                mutate(cap_size = key$cap_size)
        }) %>%
        bind_rows()

    write_cv_rows(resilts_con, "cv_external", subj_rows %>% select(-cap_size))

    # Per cap_size aggregate rows (key for answering the cap_size adequacy question)
    for (this_cap in unique(df_holdout$cap_size)) {
        df_cap <- filter(df_holdout, cap_size == this_cap)
        n_cap <- n_distinct(df_cap$global_subjid)
        m_cap <- eval_metrics(df_cap$pain_rating, df_cap$predicted)
        cap_row <- make_cv_row(
            paste0("external_cap", this_cap), level, group_value,
            this_metrics, paste0("cap_size_", this_cap), n_cap, m_cap
        )
        write_cv_rows(results_con, "cv_external", cap_row)
    }

    message(sprintf(" Holdout: r = %.3f RMSE = %.3f R^2 = %.3f MAE = %.3f p = %.3f",
                    m_all$pearson_r, m_all$rmse, m_all$r_sq,
                    m_all$mae, m_all$spearman_rho))
}

disconnect_db(data_con)
close_results_db(results_con)
message("\nExternal CV compelte. Results: cv_external table in ", gamm_cfg$results_db_path)