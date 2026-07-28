# =============================================================================
# source_core.R
# =============================================================================
# Source-space GAMM combo screen. Screens every available alpha metrics
# against every ROI present in source_metrics (by the time this runs, the
# upstream pipeline has already restricted that table to the sex-region
# pain neuromatrix - this script does not hardcode or assume which ROIs
# exist, it just reads whatever's there).
#
# Reads: data.duckdb (subjects, trials, source_metrics - read-only)
# Writes: results.duckdb (gamm_fitted table, incremental per-model writes)
#         <out_dir>/<roi>/<model_name>.rds (every fitted model object)
#
# combo_sizes = 1L (single_metric models). Flip gamm_defaults.R's
# combo_sizes to c(1L, 2L) for pairwise interaction screening - no changes
# needed here or in gamm_combo_core.R.
#
# Multiple-comparisons stance: screen uncorrected on training data; held-out
# replication is the correction (see gamm_defaults.R: apply_fdr_correction).
#
# Figure generation deliberately NOT done here - a separate script
# (export_significant_rois.R) queries gamm_fitted for hits
# and triggers rendering, per the existing significance-gating convention.
# =============================================================================

library(dplyr)
library(stringr)
library(purrr)

source("gamm_defaults.R")
source("helpers/gamm_combo_core.R")
source("helpers/db_helpers.R")

zscore <- function(x) as.numeric(scale(x))
safe_factor <- function(x) as.factor(as.character(x))

# =============================================================================
# CONNECT
# =============================================================================
data_con <- connect_data_db(gamm_cfg$data_db_path, read_only = TRUE)
results_con <- connect_results_db(gamm_cfg$results_db_path)

dir.create(gamm_cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# -- Confirm split is populated and show the train/holdout breakdown ----
split_summary <- dbGetQuery(data_con, "
    SELECT cap_size, split_group, COUNT(*) AS n
    FROM subject_split
    GROUP BY cap_size, split_group
    ORDER BY cap_size, split_group
    ")
if (nrow(split_summary) == 0) {
    stop("subject_split table is empty - run merge_core.R first.")
}
message("Subject split (GAMM fits on TRAIN subjects only):")
for (i in seq_len(nrow(split_summary))) {
    message(sprintf(" cap_size %-4s | %-8s | n = %d",
                    split_summary$cap_size[i],
                    split_summary$split_group[i],
                    split_summary$n[i]))
}
n_train <- sum(split_summary$n[split_summary$split_group == "train"])
n_holdout <- sum(split_summary$n[split_summary$split_group == "holdout"])
message(sprintf("TOTAL train = %d holdout = %d (%.1f%% / %.1f%%)",
                n_train, n_holdout,
                100 * n_train / (n_train + n_holdout),
                100 * n_holdout / (n_train + n_holdout)))

rois <- list_available_rois(data_con)
if (length(rois) == 0) stop("No ROIs found in source_metrics for training subjects - check data.duckdb.")
message("\nScreening ", length(rois), " ROI(s): ", paste(rois, collapse = ", "))

# =============================================================================
# PER-ROI LOOP
# =============================================================================
for (this_roi in rois) {
    message("\n", strrep("=", 60))
    message(" ROI: ", this_roi)
    message(strrep("=", 60))

    roi_dir <- file.path(gamm_cfg$out_dir, this_roi)
    dir.create(roi_dir, recursive = TRUE, showWarnings = FALSE)

    # -- Load + Pivot (long source_metrics -> wide, joined to subjects/trials) ----
    df_roi <- load_source_metrics_wide(data_con, this_roi)

    df_roi <- df_roi %>%
        mutate(
            experiment_id = safe_factor(experiment_id),
            subjid_uid = safe_factor(subjid_uid),
            global_subjid = safe_factor(global_subjid),
            sex = safe_factor(sex),
            cap_size = safe_factor(cap_size)
        ) %>%
        filter(
            !is.na(pain_rating), !is.na(laser_power),
            !is.na(trial_index), !is.na(global_subjid),
            !is.na(experiment_id), !is.na(age), !is.na(sex), !is.na(cap_size)
        )

    if ("fooof_r2" %in% names(df_roi)) {
        df_roi <- filter(df_roi, is.na(fooof_r2) | fooof_r2 >= gamm_cfg$min_fooof_r2)
    }

    if (nrow(df_roi) < gamm_cfg$min_trials_grp) {
        message(" Skipping: on;y ", nrow(df_roi), " trials after QC (min = ", gamm_cfg$min_trials_grp, ").")
        next
    }

    cap_counts <- df_roi %>%
        distinct(global_subjid, cap_size) %>%
        count(cap_size, name = "n_subjects")
    message(" Training subjects after QC: ", n_distinct(df_roi$global_subjid),
            " (", nrow(df_roi), " trials)")
    for (i in seq_len(nrow(cap_counts))) {
        message(" cap_size ", cap_counts$cap_size[i], ": ", cap_counts$n_subjects[i])
    }

    # -- Z-score base covariates + whichever metrics came back from the pivot ----
    df_roi <- df_roi %>%
        mutate(
            age_z = zscore(age),
            laser_power_z = zscore(laser_power),
            trial_index_z = zscore(trial_index)
        )

    # Metric columns are whatever the PIVOT produces beyong the knowns join/
    # covariate columns - this is what makes new metrics from a future
    # spectral/source pipeline revision show up here automatically, with no
    # hardcoded metric list to maintain.
    known_cols <- c("global_subjid", "subjid", "subjid_uid", "experiment_id",
                    "experiment_name", "age", "sex", "cap_size", "trial",
                    "trial_index", "laser_power", "pain_rating", "roi",
                    "age_z", "laser_power_z", "trial_index_z")
    raw_metric_cols <- setdiff(names(df_roi), known_cols)

    metrics_present <- character(0)
    for (m in raw_metric_cols) {
        if (sum(!is.na(df_roi[[m]])) > 1) {
            zcol <- paste0(m, "_z")
            df_roi[[zcol]] <- zscore(df_roi[[m]])
            metrics_present <- c(metrics_present, zcol)
        } 
    }

    if (length(metrics_present) == 0) {
        message(" No usable metrics for this ROI - skipping.")
        next
    }
    message(" Metrics available: ", paste(metrics_present, collapse = ", "))

    # -- Build base formula (covariates + REs), with level-guards so a
    #    single-level factor (e.g. only one cap_size at this ROI after QC)
    #    doesn't blow up bam() with "contrasts not defined" ----
    k_laser <- safe_k(df_roi$laser_power_z, k_max = gamm_cfg$k_max)
    k_trial <- safe_k(df_roi$trial_index_z, k_max = gamm_cfg$k_max)

    fixed_terms <- "age_z"
    if (nlevels(df_roi$sex) > 1L) fixed_terms <- c(fixed_terms, "sex")
    if (nlevels(df_roi$cap_size) > 1L) fixed_terms <- c(fixed_terms, "cap_size")

    re_terms <- character(0)
    if (n_distinct(df_roi$global_subjid) > 1L)
        re_terms <- c(re_terms, "s(global_subjid, bs = 're')")
    if (n_distinct(df_roi$experiment_id) > 1L)
        re_terms <- c(re_terms, "s(experiment_id, bs = 're')")

    base_formula <- as.formula(paste(
        "pain_rating ~",
        sprintf("s(laser_power_z, k = %d)", k_laser), "+",
        sprintf("s(trial_index_z, k = %d", k_trial), "+",
        paste(c(fixed_terms, re_terms), collapse = " + ")
    ))

    # -- Generate Combos and Fit ----
    combos <- generate_metric_combos(metrics_present, sizes = gamm_cfg$combo_sizes, include_baseline = TRUE)
    message(" Models to fit: ", nrow(combos))

    for (i in seq_len(nrow(combos))) {
        model_name <- combos$model_name[i]
        metrics <- combos$metrics[[i]]

        message(" Fitting ", model_name, " (", i, "/", nrow(combos), ")...")

        fir <- fit_gamm_combo(
            data = df_roi,
            level = "source",
            group_value = this_roi,
            metrics = metrics,
            base_formula = base_formula,
            model_name = model_name,
            out_dir = roi_dir,
            tensor_pairs = gamm_cfg$tensor_pairs,
            k_max = gamm_cfg$k_max
        )

        if (is.null(fit)) next

        # Incremental write - immediately after each successful fit.
        write_result_row(results_con, fit$result_row)
    }

    message(" Done: ", this_roi)
}

disconnect_db(data_con)
close_results_db(results_con)

message("\nSource GAMM combo screen complete.")
message("Results table: ", gamm_cfg$results_db_path, " (table: gamm_fitted)")
message("Model objects under: ", gamm_cfg$out_dir, "/<roi>/<model_name>.rds")