# =============================================================================
# channel_core.R
# -----------------------------------------------------------------------------
# Channel-space GAMM combo screen. Mirror of source_core.R — same logic,
# reads spectral_metrics (keyed by channel) instead of source_metrics
# (keyed by roi). Kept as a near-duplicate rather than a parameterized
# shared function deliberately: the two differ in their grouping column
# name (channel vs roi) and their loader (load_channel_metrics_wide vs
# load_source_metrics_wide), and forcing them through one generic function
# would trade a small amount of duplication for a worse abstraction. If
# a third grouping level is ever added, revisit this call.
#
# Reads:  data.duckdb     (subjects, trials, spectral_metrics — read-only)
# Writes: results.duckdb  (gamm_fitted table, incremental per-model writes)
#         <out_dir>/<channel>/<model_name>.rds
# =============================================================================

library(dplyr)
library(stringr)
library(purrr)

source("gamm_defaults.R")
source("helpers/gamm_combo_core.R")
source("helpers/db_helpers.R")

zscore      <- function(x) as.numeric(scale(x))
safe_factor <- function(x) as.factor(as.character(x))

# =============================================================================
# CONNECT
# =============================================================================
data_con    <- connect_data_db(gamm_cfg$data_db_path, read_only = TRUE)
results_con <- init_results_db(gamm_cfg$results_db_path)

dir.create(gamm_cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Confirm split is populated and show the train/holdout breakdown ─────────
split_summary <- dbGetQuery(data_con, "
  SELECT cap_size, split_group, COUNT(*) AS n
  FROM subject_split
  GROUP BY cap_size, split_group
  ORDER BY cap_size, split_group
")
if (nrow(split_summary) == 0) {
  stop("subject_split table is empty — run merge_core.R first.")
}
message("Subject split (GAMM fits on TRAIN subjects only):")
for (i in seq_len(nrow(split_summary))) {
  message(sprintf("  cap_size %-4s | %-8s | n = %d",
                  split_summary$cap_size[i],
                  split_summary$split_group[i],
                  split_summary$n[i]))
}
n_train   <- sum(split_summary$n[split_summary$split_group == "train"])
n_holdout <- sum(split_summary$n[split_summary$split_group == "holdout"])
message(sprintf("  TOTAL  train = %d  holdout = %d  (%.1f%% / %.1f%%)",
                n_train, n_holdout,
                100 * n_train   / (n_train + n_holdout),
                100 * n_holdout / (n_train + n_holdout)))

channels <- list_available_channels(data_con)
if (length(channels) == 0) stop("No channels found in spectral_metrics for training subjects — check data.duckdb.")
message("\nScreening ", length(channels), " channel(s).")

# =============================================================================
# PER-CHANNEL LOOP
# =============================================================================
for (this_channel in channels) {
  message("\n", strrep("=", 60))
  message("  Channel: ", this_channel)
  message(strrep("=", 60))

  chan_dir <- file.path(gamm_cfg$out_dir, this_channel)
  dir.create(chan_dir, recursive = TRUE, showWarnings = FALSE)

  # ── Load + pivot (long spectral_metrics -> wide, joined to subjects/trials) ──
  df_chan <- load_channel_metrics_wide(data_con, this_channel)

  df_chan <- df_chan %>%
    mutate(
      experiment_id = safe_factor(experiment_id),
      subjid_uid    = safe_factor(subjid_uid),
      global_subjid = safe_factor(global_subjid),
      sex           = safe_factor(sex),
      cap_size      = safe_factor(cap_size)
    ) %>%
    filter(
      !is.na(pain_rating), !is.na(laser_power),
      !is.na(trial_index), !is.na(global_subjid),
      !is.na(experiment_id), !is.na(age), !is.na(sex), !is.na(cap_size)
    )

  if ("fooof_r2" %in% names(df_chan)) {
    df_chan <- filter(df_chan, is.na(fooof_r2) | fooof_r2 >= gamm_cfg$min_fooof_r2)
  }

  if (nrow(df_chan) < gamm_cfg$min_trials_grp) {
    message("  Skipping: only ", nrow(df_chan), " trials after QC (min = ",
            gamm_cfg$min_trials_grp, ").")
    next
  }

  cap_counts <- df_chan %>%
    distinct(global_subjid, cap_size) %>%
    count(cap_size, name = "n_subjects")
  message("  Training subjects after QC: ", n_distinct(df_chan$global_subjid),
          " (", nrow(df_chan), " trials)")
  for (i in seq_len(nrow(cap_counts))) {
    message("    cap_size ", cap_counts$cap_size[i], ": ", cap_counts$n_subjects[i])
  }

  # ── Z-score base covariates + whichever metrics came back from the pivot ──
  df_chan <- df_chan %>%
    mutate(
      age_z         = zscore(age),
      laser_power_z = zscore(laser_power),
      trial_index_z = zscore(trial_index)
    )

  known_cols <- c("global_subjid", "subjid", "subjid_uid", "experiment_id",
                   "experiment_name", "age", "sex", "cap_size", "trial",
                   "trial_index", "laser_power", "pain_rating", "channel",
                   "age_z", "laser_power_z", "trial_index_z")
  raw_metric_cols <- setdiff(names(df_chan), known_cols)

  metrics_present <- character(0)
  for (m in raw_metric_cols) {
    if (sum(!is.na(df_chan[[m]])) > 1) {
      zcol <- paste0(m, "_z")
      df_chan[[zcol]] <- zscore(df_chan[[m]])
      metrics_present <- c(metrics_present, zcol)
    }
  }

  if (length(metrics_present) == 0) {
    message("  No usable metrics for this channel — skipping.")
    next
  }
  message("  Metrics available: ", paste(metrics_present, collapse = ", "))

  # ── Build base formula (covariates + REs), with level-guards ───────────────
  k_laser <- safe_k(df_chan$laser_power_z, k_max = gamm_cfg$k_max)
  k_trial <- safe_k(df_chan$trial_index_z, k_max = gamm_cfg$k_max)

  fixed_terms <- "age_z"
  if (nlevels(df_chan$sex)      > 1L) fixed_terms <- c(fixed_terms, "sex")
  if (nlevels(df_chan$cap_size) > 1L) fixed_terms <- c(fixed_terms, "cap_size")

  re_terms <- character(0)
  if (n_distinct(df_chan$global_subjid) > 1L)
    re_terms <- c(re_terms, "s(global_subjid, bs = 're')")
  if (n_distinct(df_chan$experiment_id) > 1L)
    re_terms <- c(re_terms, "s(experiment_id,  bs = 're')")

  base_formula <- as.formula(paste(
    "pain_rating ~",
    sprintf("s(laser_power_z, k = %d)", k_laser), "+",
    sprintf("s(trial_index_z, k = %d)", k_trial), "+",
    paste(c(fixed_terms, re_terms), collapse = " + ")
  ))

  # ── Generate combos and fit ─────────────────────────────────────────────
  combos <- generate_metric_combos(metrics_present, sizes = gamm_cfg$combo_sizes,
                                    include_baseline = TRUE)
  message("  Models to fit: ", nrow(combos))

  for (i in seq_len(nrow(combos))) {
    model_name <- combos$model_name[i]
    metrics    <- combos$metrics[[i]]

    message("  Fitting ", model_name, " (", i, "/", nrow(combos), ")...")

    fit <- fit_gamm_combo(
      data         = df_chan,
      level        = "channel",
      group_value  = this_channel,
      metrics      = metrics,
      base_formula = base_formula,
      model_name   = model_name,
      out_dir      = chan_dir,
      tensor_pairs = gamm_cfg$tensor_pairs,
      k_max        = gamm_cfg$k_max
    )

    if (is.null(fit)) next

    write_result_row(results_con, fit$result_row)
  }

  message("  Done: ", this_channel)
}

disconnect_db(data_con)
close_results_db(results_con)

message("\nChannel GAMM combo screen complete.")
message("Results table: ", gamm_cfg$results_db_path, "  (table: gamm_fitted)")
message("Model objects under: ", gamm_cfg$out_dir, "/<channel>/<model_name>.rds")