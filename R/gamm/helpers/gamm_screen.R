# =============================================================================
# gamm_screen.R
# -----------------------------------------------------------------------------
# run_gamm_screen(level): the per-unit GAMM combo screen shared by
# source_core.R (level = "source", unit = ROI) and channel_core.R
# (level = "channel", unit = channel label).
#
# V2 (Oct 2026): the screen used to be written out twice, once per level.
# The two copies had drifted (source_core.R had `fir <-` vs `fit`, so no
# result was ever written; a missing ")" in s(trial_index_z, ...; and a
# call to connect_results_db() that never created gamm_fitted), while
# channel_core.R was correct. One implementation removes that failure mode;
# the level-specific parts (table, unit column, output folder) come from
# level_spec() in db_helpers.R.
#
# Per unit:
#   1. Load the screened metrics (gamm_cfg$metric_include) for TRAIN
#      subjects, pivoted wide and joined to subjects + trials.
#      Only subjects that passed the merge-stage alignment check have rows.
#   2. Drop rows with missing outcome/covariates; optional p5_flag and
#      subject-level FOOOF filters.
#   3. z-score covariates and metrics (constant or all-NA metrics dropped).
#      With metric_decomposition = "between_within" (default) each metric is
#      split into the subject's mean (_between_z) and the trial deviation
#      from it (_within_z); both enter the model, and the p-value written to
#      gamm_fitted is for primary_part (default "between": does a person's
#      typical value predict how much pain they report?).
#   4. Base model: pain_rating ~ s(laser_power_z) + s(trial_index_z) + age_z
#      (a covariate with < 3 unique values enters linearly)
#      [+ sex] [+ cap_size] [+ s(global_subjid, bs="re")]
#      [+ s(experiment_id, bs="re")]
#   5. Fit the baseline and every metric combo; write each row to
#      gamm_fitted immediately.
# =============================================================================

library(dplyr)
library(stringr)
library(purrr)

zscore      <- function(x) as.numeric(scale(x))
`%||%`      <- function(a, b) if (is.null(a)) b else a
safe_factor <- function(x) as.factor(as.character(x))

print_split_summary <- function(data_con) {
  split_summary <- dbGetQuery(data_con, "
    SELECT cap_size, split_group, COUNT(*) AS n
    FROM subject_split
    GROUP BY cap_size, split_group
    ORDER BY cap_size, split_group")
  if (nrow(split_summary) == 0) stop("subject_split table is empty - run merge_core.R first.")

  message("Subject split (GAMM fits on TRAIN subjects only):")
  for (i in seq_len(nrow(split_summary))) {
    message(sprintf("  cap_size %-4s | %-8s | n = %d",
                    split_summary$cap_size[i], split_summary$split_group[i],
                    split_summary$n[i]))
  }
  n_train   <- sum(split_summary$n[split_summary$split_group == "train"])
  n_holdout <- sum(split_summary$n[split_summary$split_group == "holdout"])
  message(sprintf("  TOTAL  train = %d  holdout = %d  (%.1f%% / %.1f%%)",
                  n_train, n_holdout,
                  100 * n_train / (n_train + n_holdout),
                  100 * n_holdout / (n_train + n_holdout)))

  # Subjects in the split but with no EEG rows (alignment exclusions,
  # not yet processed) are listed so the effective N is never a surprise.
  if (dbExistsTable(data_con, "alignment_qc")) {
    qc <- dbGetQuery(data_con, "
      SELECT q.status, COUNT(*) AS n
      FROM alignment_qc q JOIN subject_split s USING (global_subjid)
      WHERE s.split_group = 'train'
      GROUP BY q.status ORDER BY q.status")
    if (nrow(qc) > 0) {
      message("  Train subjects by alignment status: ",
              paste(sprintf("%s = %d", qc$status, qc$n), collapse = ", "))
    }
  }
}

run_gamm_screen <- function(level = c("source", "channel"), cfg = gamm_cfg) {
  level    <- match.arg(level)
  sp       <- level_spec(level)
  decomposition <- match.arg(cfg$metric_decomposition %||% "between_within",
                             c("between_within", "pooled"))
  primary_part  <- match.arg(cfg$primary_part %||% "between", c("between", "within"))
  unit_col <- sp$unit_col
  unit_lbl <- if (level == "source") "ROI" else "Channel"

  data_con    <- connect_data_db(cfg$data_db_path, read_only = TRUE)
  on.exit(try(disconnect_db(data_con), silent = TRUE), add = TRUE)
  results_con <- init_results_db(cfg$results_db_path)
  on.exit(try(close_results_db(results_con), silent = TRUE), add = TRUE)

  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  print_split_summary(data_con)

  # -- Units ----
  units <- select_units(list_available_units(data_con, level),
                        if (level == "source") cfg$roi_include else cfg$channel_include)
  if (length(units) == 0) stop("No ", unit_lbl, "s found in ", sp$table,
                               " for training subjects - check data.duckdb / *_include.")
  message("\nScreening ", length(units), " ", unit_lbl, "(s): ", paste(units, collapse = ", "))

  # -- Metrics ----
  available <- list_available_metrics(data_con, level)
  metrics_sel <- select_metrics(available, cfg$metric_include, cfg$metric_exclude)
  if (length(metrics_sel) == 0) stop("No metric in ", sp$table, " matches gamm_cfg$metric_include.")
  message("Metrics selected: ", length(metrics_sel), " of ", length(available), " in ", sp$table)
  message("Metric decomposition: ", decomposition,
          if (decomposition == "between_within")
            paste0(" (p-value reported for the ", primary_part, "-subject part)") else "")
  skipped <- setdiff(available, metrics_sel)
  if (length(skipped) > 0) message("  Not screened: ", paste(skipped, collapse = ", "))

  load_cols <- metrics_sel
  if (isTRUE(cfg$exclude_p5_flag)) {
    if (!"p5_flag" %in% available) stop("exclude_p5_flag = TRUE but p5_flag is not in ", sp$table)
    load_cols <- c(load_cols, "p5_flag")
  }

  # -- Subject-level FOOOF QC (optional) ----
  fooof_bad <- NULL
  if (is.numeric(cfg$min_fooof_r2)) {
    fq <- load_subject_metrics_wide(data_con, level, unit = units, metrics = "fooof_r2")
    if (!"fooof_r2" %in% names(fq)) {
      message("  [WARN] min_fooof_r2 set but no fooof_r2 in the subject-level table; filter skipped.")
    } else {
      fooof_bad <- fq %>% filter(is.na(fooof_r2) | fooof_r2 < cfg$min_fooof_r2) %>%
        select(global_subjid, all_of(unit_col))
      message("  FOOOF QC: ", nrow(fooof_bad), " subject x ", unit_lbl,
              " pairs below r2 = ", cfg$min_fooof_r2, " will be dropped.")
    }
  }

  # ===========================================================================
  # PER-UNIT LOOP
  # ===========================================================================
  for (this_unit in units) {
    message("\n", strrep("=", 60))
    message("  ", unit_lbl, ": ", this_unit)
    message(strrep("=", 60))

    unit_dir <- file.path(cfg$out_dir, level, this_unit)
    dir.create(unit_dir, recursive = TRUE, showWarnings = FALSE)

    df <- load_metrics_wide(data_con, level, this_unit, metrics = load_cols, split = "train")

    df <- df %>%
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

    if (isTRUE(cfg$exclude_p5_flag) && "p5_flag" %in% names(df)) {
      n0 <- nrow(df)
      df <- df %>% filter(is.na(p5_flag) | p5_flag < 0.5)
      message("  p5_flag: dropped ", n0 - nrow(df), " of ", n0, " trials.")
    }
    if (!is.null(fooof_bad)) {
      bad_here <- as.character(fooof_bad$global_subjid[fooof_bad[[unit_col]] == this_unit])
      if (length(bad_here) > 0) {
        df <- df %>% filter(!as.character(global_subjid) %in% bad_here)
        message("  FOOOF QC: dropped ", length(bad_here), " subject(s).")
      }
    }
    df <- droplevels(df)

    if (nrow(df) < cfg$min_trials_grp) {
      message("  Skipping: only ", nrow(df), " trials after QC (min = ", cfg$min_trials_grp, ").")
      next
    }

    cap_counts <- df %>% distinct(global_subjid, cap_size) %>% count(cap_size, name = "n_subjects")
    message("  Training subjects after QC: ", n_distinct(df$global_subjid),
            " (", nrow(df), " trials)")
    for (i in seq_len(nrow(cap_counts))) {
      message("    cap_size ", cap_counts$cap_size[i], ": ", cap_counts$n_subjects[i])
    }

    # -- z-score covariates + metrics ----
    df <- df %>%
      mutate(
        age_z         = zscore(age),
        laser_power_z = zscore(laser_power),
        trial_index_z = zscore(trial_index)
      )

    # Mean / SD behind every z-scored column, saved with each model so
    # plots and holdout predictions can use the original units.
    zscale <- lapply(c(age = "age", laser_power = "laser_power", trial_index = "trial_index"),
                     function(cn) c(mean = mean(df[[cn]], na.rm = TRUE),
                                    sd = stats::sd(df[[cn]], na.rm = TRUE)))

    split_bw <- identical(decomposition, "between_within")
    metrics_present <- character(0)
    dropped <- character(0)
    for (m in intersect(metrics_sel, names(df))) {
      x <- df[[m]]
      if (sum(!is.na(x)) > 1 && isTRUE(stats::sd(x, na.rm = TRUE) > 0)) {
        zcol <- paste0(m, "_z")
        df[[zcol]] <- zscore(x)
        zscale[[m]] <- c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE))
        if (split_bw) {
          df <- add_between_within(df, zcol)
          pv <- metric_parts(zcol, decomposition)
          b_ok <- isTRUE(stats::sd(unique(df[!is.na(df[[pv[["between"]]]]), c("global_subjid", pv[["between"]])])[[2]]) > 0)
          w_ok <- isTRUE(stats::sd(df[[pv[["within"]]]], na.rm = TRUE) > 0)
          if (!(b_ok && w_ok)) { dropped <- c(dropped, m); next }
        }
        metrics_present <- c(metrics_present, zcol)
      }
    }
    if (length(dropped) > 0)
      message("  No between- or within-subject variance, skipped: ", paste(dropped, collapse = ", "))
    if (length(metrics_present) == 0) {
      message("  No usable metrics for this ", unit_lbl, " - skipping.")
      next
    }
    message("  Metrics usable: ", length(metrics_present))

    # -- Base formula (covariates + REs), with level guards ----

    fixed_terms <- "age_z"
    if (nlevels(df$sex)      > 1L) fixed_terms <- c(fixed_terms, "sex")
    if (nlevels(df$cap_size) > 1L) fixed_terms <- c(fixed_terms, "cap_size")

    re_terms <- character(0)
    if (n_distinct(df$global_subjid) > 1L) re_terms <- c(re_terms, "s(global_subjid, bs = 're')")
    if (n_distinct(df$experiment_id) > 1L) re_terms <- c(re_terms, "s(experiment_id, bs = 're')")

    base_formula <- as.formula(paste(
      "pain_rating ~",
      smooth_or_linear("laser_power_z", df$laser_power_z, k_max = cfg$k_max), "+",
      smooth_or_linear("trial_index_z", df$trial_index_z, k_max = cfg$k_max), "+",
      paste(c(fixed_terms, re_terms), collapse = " + ")
    ))

    # -- Combos ----
    combos <- generate_metric_combos(metrics_present, sizes = cfg$combo_sizes,
                                     include_baseline = TRUE)
    message("  Models to fit: ", nrow(combos))

    n_ok <- 0L
    for (i in seq_len(nrow(combos))) {
      model_name <- combos$model_name[i]
      metrics    <- combos$metrics[[i]]

      fit <- fit_gamm_combo(
        data         = df,
        level        = level,
        group_value  = this_unit,
        metrics      = metrics,
        base_formula = base_formula,
        model_name   = model_name,
        out_dir      = unit_dir,
        tensor_pairs = cfg$tensor_pairs,
        k_max        = cfg$k_max,
        feature_set  = cfg$feature_set,
        nthreads     = cfg$nthreads,
        discrete     = !isFALSE(cfg$bam_discrete),
        zscale       = zscale,
        decomposition = decomposition,
        primary_part = primary_part
      )
      if (is.null(fit)) next

      write_result_row(results_con, fit$result_row)
      n_ok <- n_ok + 1L
    }

    message("  Done: ", this_unit, " (", n_ok, "/", nrow(combos), " models written)")
  }

  message("\n", unit_lbl, " GAMM combo screen complete.")
  message("Results table: ", cfg$results_db_path, "  (table: gamm_fitted, feature_set = ",
          cfg$feature_set, ")")
  message("Model objects under: ", file.path(cfg$out_dir, level), "/<unit>/<model_name>.rds")
  invisible(NULL)
}