# =============================================================================
# db_write_helpers.R
# -----------------------------------------------------------------------------
# Schema creation + upsert helpers for data.duckdb. Used by merge_core.R.
# =============================================================================

library(DBI)
library(duckdb)
library(dplyr)

connect_or_create_data_db <- function(db_path) {
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
  dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = FALSE)
}

init_data_schema <- function(con) {
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS subjects (
      global_subjid   VARCHAR PRIMARY KEY,
      subjid          INTEGER,
      subjid_uid      VARCHAR,
      experiment_id   INTEGER,
      experiment_name VARCHAR,
      age             DOUBLE,
      sex             VARCHAR,
      cap_size        VARCHAR
    )
  ")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS trials (
      global_subjid VARCHAR,
      trial         INTEGER,
      trial_index   INTEGER,
      laser_power   DOUBLE,
      pain_rating   DOUBLE,
      PRIMARY KEY (global_subjid, trial)
    )
  ")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS source_metrics (
      global_subjid VARCHAR,
      trial         INTEGER,
      roi           VARCHAR,
      metric_name   VARCHAR,
      metric_value  DOUBLE,
      PRIMARY KEY (global_subjid, trial, roi, metric_name)
    )
  ")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS spectral_metrics (
      global_subjid VARCHAR,
      trial         INTEGER,
      channel       VARCHAR,
      metric_name   VARCHAR,
      metric_value  DOUBLE,
      PRIMARY KEY (global_subjid, trial, channel, metric_name)
    )
  ")

  invisible(NULL)
}

# =============================================================================
# upsert_subjects() / upsert_trials()
# -----------------------------------------------------------------------------
# Delete-then-insert per global_subjid, same pattern as gamm_fitted's upsert.
# Subjects/trials are small enough that re-ingesting a whole subject's worth
# of rows on re-run is cheap and avoids partial-update bugs.
# =============================================================================
upsert_subjects <- function(con, df) {
  if (nrow(df) == 0) return(invisible(NULL))
  ids <- unique(df$global_subjid)
  placeholders <- paste(rep("?", length(ids)), collapse = ",")
  dbExecute(con, sprintf("DELETE FROM subjects WHERE global_subjid IN (%s)", placeholders),
            params = as.list(ids))
  dbWriteTable(con, "subjects", df, append = TRUE)
  invisible(NULL)
}

upsert_trials <- function(con, df) {
  if (nrow(df) == 0) return(invisible(NULL))
  ids <- unique(df$global_subjid)
  placeholders <- paste(rep("?", length(ids)), collapse = ",")
  dbExecute(con, sprintf("DELETE FROM trials WHERE global_subjid IN (%s)", placeholders),
            params = as.list(ids))
  dbWriteTable(con, "trials", df, append = TRUE)
  invisible(NULL)
}

# =============================================================================
# upsert_metrics_long()
# -----------------------------------------------------------------------------
# Shared by source_metrics and spectral_metrics — same delete-then-insert
# pattern, parameterized by table name and grouping column name (roi vs
# channel), since the two tables are structurally identical apart from that
# one column.
# =============================================================================
upsert_metrics_long <- function(con, table_name, df) {
  if (nrow(df) == 0) return(invisible(NULL))
  ids <- unique(df$global_subjid)
  placeholders <- paste(rep("?", length(ids)), collapse = ",")
  dbExecute(con, sprintf("DELETE FROM %s WHERE global_subjid IN (%s)", table_name, placeholders),
            params = as.list(ids))
  dbWriteTable(con, table_name, df, append = TRUE)
  invisible(NULL)
}