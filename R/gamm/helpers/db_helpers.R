# =============================================================================
# db_helpers.R
# -----------------------------------------------------------------------------
# Connection + read helpers for the two-database architecture:
#
#   data.duckdb     subjects, trials, spectral_metrics, source_metrics
#                   (spectral_metrics/source_metrics are LONG/tidy format:
#                    one row per global_subjid x trial x channel-or-roi x
#                    metric_name. New metrics from pipeline revisions just
#                    appear as new metric_name values — no schema change,
#                    no ALTER TABLE, source/channel_core.R need no edits.)
#
#   results.duckdb  gamm_fitted (written by gamm_combo_core.R's
#                   init_results_db()/write_result_row(), unchanged from
#                   before)
#
# This file owns the long-to-wide PIVOT so source_core.R/channel_core.R
# never see long-format data or write SQL themselves — they call
# load_source_metrics_wide(con, roi) or load_channel_metrics_wide(con, channel)
# and get back a normal data frame, same shape as the old
# source_pain_master.csv read used to produce.
# =============================================================================

library(DBI)
library(duckdb)
library(dplyr)

# =============================================================================
# connect_data_db() / connect_results_db()
# -----------------------------------------------------------------------------
# read_only = TRUE on the data connection by design — the gamm stage never
# writes to data.duckdb, only reads. This is a real guard, not just
# documentation: DuckDB enforces it, so an accidental write attempt from
# this stage fails loudly instead of silently corrupting upstream data
# that the merge stage owns.
# =============================================================================
connect_data_db <- function(db_path, read_only = TRUE) {
  if (!file.exists(db_path)) {
    stop("data.duckdb not found at: ", db_path,
         "\nRun the merge stage (build_data_db.R) first.")
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

# =============================================================================
# list_available_rois() / list_available_channels()
# -----------------------------------------------------------------------------
# Returns only ROIs/channels that have at least one training-set subject's
# data — prevents the fitting loop from attempting a model on a group where
# every subject happened to land in holdout (possible at small N, though
# quota-based stratification makes it unlikely).
# =============================================================================
list_available_rois <- function(con) {
  dbGetQuery(con, "
    SELECT DISTINCT sm.roi
    FROM source_metrics sm
    JOIN subject_split ss USING (global_subjid)
    WHERE ss.split_group = 'train'
    ORDER BY sm.roi
  ")$roi
}

list_available_channels <- function(con) {
  dbGetQuery(con, "
    SELECT DISTINCT sm.channel
    FROM spectral_metrics sm
    JOIN subject_split ss USING (global_subjid)
    WHERE ss.split_group = 'train'
    ORDER BY sm.channel
  ")$channel
}

# =============================================================================
# load_source_metrics_wide() / load_channel_metrics_wide()
# -----------------------------------------------------------------------------
# Filters to training subjects only (split_group = 'train') via JOIN against
# subject_split before pivoting — held-out subjects' data is physically
# excluded from every data frame handed to the GAMM fitting stage.
#
# This is the correct place to enforce the train/holdout boundary: the
# GAMM scripts themselves never see held-out data, so there's no risk of
# accidentally fitting on it if someone modifies source_core.R later.
#
# Uses `w.* EXCLUDE (...)` to avoid duplicate join-key columns (tested and
# confirmed necessary — SELECT s.*, w.* produces global_subjid_1 duplicates
# that R silently mangles via make.unique()).
# =============================================================================
load_source_metrics_wide <- function(con, roi) {
  query <- sprintf("
    WITH train_subjects AS (
      SELECT global_subjid FROM subject_split WHERE split_group = 'train'
    ),
    wide_metrics AS (
      PIVOT (
        SELECT sm.* FROM source_metrics sm
        JOIN train_subjects ts USING (global_subjid)
        WHERE sm.roi = '%s'
      )
      ON metric_name
      USING FIRST(metric_value)
      GROUP BY global_subjid, trial, roi
    )
    SELECT
      s.global_subjid, s.subjid, s.subjid_uid, s.experiment_id, s.experiment_name,
      s.age, s.sex, s.cap_size,
      t.trial, t.trial_index, t.laser_power, t.pain_rating,
      w.roi,
      w.* EXCLUDE (global_subjid, trial, roi)
    FROM wide_metrics w
    JOIN subjects s USING (global_subjid)
    JOIN trials   t USING (global_subjid, trial)
  ", roi)

  as_tibble(dbGetQuery(con, query))
}

load_channel_metrics_wide <- function(con, channel) {
  query <- sprintf("
    WITH train_subjects AS (
      SELECT global_subjid FROM subject_split WHERE split_group = 'train'
    ),
    wide_metrics AS (
      PIVOT (
        SELECT sm.* FROM spectral_metrics sm
        JOIN train_subjects ts USING (global_subjid)
        WHERE sm.channel = '%s'
      )
      ON metric_name
      USING FIRST(metric_value)
      GROUP BY global_subjid, trial, channel
    )
    SELECT
      s.global_subjid, s.subjid, s.subjid_uid, s.experiment_id, s.experiment_name,
      s.age, s.sex, s.cap_size,
      t.trial, t.trial_index, t.laser_power, t.pain_rating,
      w.channel,
      w.* EXCLUDE (global_subjid, trial, channel)
    FROM wide_metrics w
    JOIN subjects s USING (global_subjid)
    JOIN trials   t USING (global_subjid, trial)
  ", channel)

  as_tibble(dbGetQuery(con, query))
}

# =============================================================================
# load_holdout_metrics_wide()
# -----------------------------------------------------------------------------
# Counterpart to load_source_metrics_wide() for the CV evaluation stage —
# returns held-out subjects only (split_group = 'holdout'). Not used by
# source_core.R/channel_core.R; reserved for the CV scripts (not yet
# written) so the query pattern is in one place rather than duplicated
# across CV and fitting scripts.
# =============================================================================
load_holdout_source_metrics_wide <- function(con, roi) {
  query <- sprintf("
    WITH holdout_subjects AS (
      SELECT global_subjid FROM subject_split WHERE split_group = 'holdout'
    ),
    wide_metrics AS (
      PIVOT (
        SELECT sm.* FROM source_metrics sm
        JOIN holdout_subjects hs USING (global_subjid)
        WHERE sm.roi = '%s'
      )
      ON metric_name
      USING FIRST(metric_value)
      GROUP BY global_subjid, trial, roi
    )
    SELECT
      s.global_subjid, s.subjid, s.subjid_uid, s.experiment_id, s.experiment_name,
      s.age, s.sex, s.cap_size,
      t.trial, t.trial_index, t.laser_power, t.pain_rating,
      w.roi,
      w.* EXCLUDE (global_subjid, trial, roi)
    FROM wide_metrics w
    JOIN subjects s USING (global_subjid)
    JOIN trials   t USING (global_subjid, trial)
  ", roi)

  as_tibble(dbGetQuery(con, query))
}

# =============================================================================
# list_available_metrics()
# -----------------------------------------------------------------------------
# Replaces the old candidate_metrics_raw hardcoded vector. Queries
# metric_name directly from the long-format table (training subjects only)
# so newly-added metrics appear automatically without any hardcoded list.
# =============================================================================
list_available_metrics <- function(con, level = c("source", "channel")) {
  level <- match.arg(level)
  tbl <- if (level == "source") "source_metrics" else "spectral_metrics"
  dbGetQuery(con, sprintf("
    SELECT DISTINCT sm.metric_name
    FROM %s sm
    JOIN subject_split ss USING (global_subjid)
    WHERE ss.split_group = 'train'
    ORDER BY sm.metric_name
  ", tbl))$metric_name
}