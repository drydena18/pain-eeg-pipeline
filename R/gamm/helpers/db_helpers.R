# =============================================================================
# db_helpers.R
# -----------------------------------------------------------------------------
# Connection + read helpers for the two-database architecture:
#
#   data.duckdb     subjects, trials, epochs, alignment_qc,
#                   spectral_metrics / source_metrics           (trial level)
#                   spectral_subject_metrics / source_subject_metrics
#                                                               (subject level)
#                   subject_split
#                   All metric tables are LONG/tidy: one row per
#                   global_subjid x [trial x] channel-or-roi x metric_name.
#                   New metrics from pipeline revisions appear as new
#                   metric_name values - no schema change.
#
#   results.duckdb  gamm_fitted (gamm_combo_core.R: init_results_db /
#                   write_result_row)
#
# This file owns the long-to-wide PIVOT so the GAMM scripts never see
# long-format data or write SQL themselves.
#
# V2 (Oct 2026):
#   - One generic loader, load_metrics_wide(con, level, unit, metrics, split),
#     behind the old load_source_metrics_wide / load_channel_metrics_wide /
#     load_holdout_source_metrics_wide names (kept as thin wrappers).
#   - `metrics` is pushed into SQL (metric_name IN (...)), so only the
#     screened metrics are pivoted, not all 58 v3.1 columns.
#   - String values are quoted with dbQuoteString (an ROI or channel label
#     containing a quote no longer breaks the query).
#   - load_subject_metrics_wide(): subject-level tables (TVI, ITC, FOOOF,
#     trial-mean summaries), one row per subject x unit.
#   - connect_results_db() is kept for older callers, but the GAMM and
#     validate stages now open results.duckdb through init_results_db()
#     (gamm_combo_core.R), which also creates gamm_fitted and refuses a
#     pre-v3.1 table.
# =============================================================================

library(DBI)
library(duckdb)
library(dplyr)

# =============================================================================
# connect_data_db()
# -----------------------------------------------------------------------------
# read_only = TRUE on the data connection by design — the gamm stage never
# writes to data.duckdb, only reads. DuckDB enforces it.
# =============================================================================
connect_data_db <- function(db_path, read_only = TRUE) {
  if (!file.exists(db_path)) {
    stop("data.duckdb not found at: ", db_path,
         "\nRun the merge stage (R/merge/merge_core.R) first.")
  }
  dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = read_only)
}

connect_results_db <- function(db_path) {
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
  dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = FALSE)
}

disconnect_db <- function(con) {
  dbDisconnect(con, shutdown = TRUE)
}

# -- Level lookup ----
level_spec <- function(level = c("source", "channel"), subject_level = FALSE) {
  level <- match.arg(level)
  if (level == "source") {
    list(table = if (subject_level) "source_subject_metrics" else "source_metrics",
         unit_col = "roi")
  } else {
    list(table = if (subject_level) "spectral_subject_metrics" else "spectral_metrics",
         unit_col = "channel")
  }
}

split_filter_sql <- function(split) {
  split <- match.arg(split, c("train", "holdout", "all"))
  if (split == "all") return("SELECT global_subjid FROM subject_split")
  sprintf("SELECT global_subjid FROM subject_split WHERE split_group = '%s'", split)
}

sql_in <- function(con, x) paste(vapply(x, function(v) as.character(dbQuoteString(con, v)),
                                        character(1)), collapse = ", ")

# =============================================================================
# list_available_units() / list_available_metrics()
# -----------------------------------------------------------------------------
# Only units/metrics with at least one subject in the requested split, so
# the fitting loop never attempts a model on a group with no training data.
# =============================================================================
list_available_units <- function(con, level = c("source", "channel"), split = "train") {
  sp <- level_spec(level)
  dbGetQuery(con, sprintf("
    SELECT DISTINCT m.%s AS unit
    FROM %s m
    WHERE m.global_subjid IN (%s)
    ORDER BY unit", sp$unit_col, sp$table, split_filter_sql(split)))$unit
}

list_available_rois     <- function(con) list_available_units(con, "source")
list_available_channels <- function(con) list_available_units(con, "channel")

list_available_metrics <- function(con, level = c("source", "channel"), split = "train",
                                   subject_level = FALSE) {
  sp <- level_spec(level, subject_level)
  dbGetQuery(con, sprintf("
    SELECT DISTINCT m.metric_name
    FROM %s m
    WHERE m.global_subjid IN (%s)
    ORDER BY m.metric_name", sp$table, split_filter_sql(split)))$metric_name
}

# Apply include/exclude regexes to a vector of metric names.
select_metrics <- function(available, include, exclude = character(0)) {
  keep <- Reduce(`|`, lapply(include, function(p) grepl(p, available)), FALSE)
  if (length(exclude) > 0) {
    drop <- Reduce(`|`, lapply(exclude, function(p) grepl(p, available)), FALSE)
    keep <- keep & !drop
  }
  available[keep]
}

# Apply an include spec (NULL, exact names, or a single regex) to units.
select_units <- function(available, include) {
  if (is.null(include) || length(include) == 0) return(available)
  if (length(include) == 1 && grepl("[\\^\\$\\|\\(\\[\\*\\+\\?]", include)) {
    return(available[grepl(include, available)])
  }
  missing <- setdiff(include, available)
  if (length(missing) > 0) message("  [WARN] Requested unit(s) not in database: ",
                                   paste(missing, collapse = ", "))
  intersect(available, include)
}

# =============================================================================
# load_metrics_wide()
# -----------------------------------------------------------------------------
# Trial-level metrics for one unit, pivoted wide and joined to subjects and
# trials. The split filter is applied BEFORE pivoting — held-out subjects'
# data is physically excluded from every data frame handed to the fitting
# stage, so the GAMM scripts never see held-out data.
#
# metrics : character vector of metric_name values to pivot (NULL = all).
#
# `w.* EXCLUDE (...)` avoids duplicate join-key columns (SELECT s.*, w.*
# produces global_subjid_1 duplicates that R silently mangles).
# =============================================================================
load_metrics_wide <- function(con, level = c("source", "channel"), unit,
                              metrics = NULL, split = "train") {
  sp <- level_spec(level)
  metric_sql <- if (is.null(metrics)) "" else
    sprintf("AND m.metric_name IN (%s)", sql_in(con, metrics))

  query <- sprintf("
    WITH wide_metrics AS (
      PIVOT (
        SELECT m.global_subjid, m.trial, m.%1$s, m.metric_name, m.metric_value
        FROM %2$s m
        WHERE m.%1$s = %3$s
          AND m.global_subjid IN (%4$s)
          %5$s
      )
      ON metric_name
      USING FIRST(metric_value)
      GROUP BY global_subjid, trial, %1$s
    )
    SELECT
      s.global_subjid, s.subjid, s.subjid_uid, s.experiment_id, s.experiment_name,
      s.age, s.sex, s.cap_size,
      t.trial, t.trial_index, t.laser_power, t.pain_rating,
      w.%1$s,
      w.* EXCLUDE (global_subjid, trial, %1$s)
    FROM wide_metrics w
    JOIN subjects s USING (global_subjid)
    JOIN trials   t USING (global_subjid, trial)
    ORDER BY s.global_subjid, t.trial",
    sp$unit_col, sp$table, as.character(dbQuoteString(con, unit)),
    split_filter_sql(split), metric_sql)

  as_tibble(dbGetQuery(con, query))
}

# -- Backwards-compatible wrappers ----
load_source_metrics_wide  <- function(con, roi, metrics = NULL)
  load_metrics_wide(con, "source", roi, metrics, "train")
load_channel_metrics_wide <- function(con, channel, metrics = NULL)
  load_metrics_wide(con, "channel", channel, metrics, "train")
load_holdout_source_metrics_wide <- function(con, roi, metrics = NULL)
  load_metrics_wide(con, "source", roi, metrics, "holdout")

# =============================================================================
# load_subject_metrics_wide()
# -----------------------------------------------------------------------------
# Subject-level tables (one row per subject x unit): TVI, ITC, FOOOF and the
# trial-mean window summaries. unit = NULL returns every unit.
# Joined to subjects only (there is no trial at this level).
# =============================================================================
load_subject_metrics_wide <- function(con, level = c("source", "channel"), unit = NULL,
                                      metrics = NULL, split = "train") {
  sp <- level_spec(level, subject_level = TRUE)
  unit_sql <- if (is.null(unit)) "" else
    sprintf("AND m.%s IN (%s)", sp$unit_col, sql_in(con, unit))
  metric_sql <- if (is.null(metrics)) "" else
    sprintf("AND m.metric_name IN (%s)", sql_in(con, metrics))

  query <- sprintf("
    WITH wide_metrics AS (
      PIVOT (
        SELECT m.global_subjid, m.%1$s, m.metric_name, m.metric_value
        FROM %2$s m
        WHERE m.global_subjid IN (%3$s) %4$s %5$s
      )
      ON metric_name
      USING FIRST(metric_value)
      GROUP BY global_subjid, %1$s
    )
    SELECT
      s.global_subjid, s.subjid, s.subjid_uid, s.experiment_id, s.experiment_name,
      s.age, s.sex, s.cap_size,
      w.%1$s,
      w.* EXCLUDE (global_subjid, %1$s)
    FROM wide_metrics w
    JOIN subjects s USING (global_subjid)
    ORDER BY s.global_subjid, w.%1$s",
    sp$unit_col, sp$table, split_filter_sql(split), unit_sql, metric_sql)

  as_tibble(dbGetQuery(con, query))
}