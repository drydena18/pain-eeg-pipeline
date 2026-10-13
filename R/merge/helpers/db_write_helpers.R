# =============================================================================
# db_write_helpers.R
# -----------------------------------------------------------------------------
# Schema creation, rebuild and write helpers for data.duckdb (merge stage).
#
# V2 (Oct 2026):
#   - New tables: epochs, alignment_qc, spectral_subject_metrics,
#     source_subject_metrics, ingest_log.
#   - drop_derived_tables(): drops everything except subject_split.
#   - ingest_feature_csv(): DuckDB reads one per-subject CSV, unpivots it to
#     long format and inserts it in a single SQL statement (all-varchar read
#     + TRY_CAST, so a column that is entirely NaN/empty in one file cannot
#     break type inference). NaN / Inf / empty values are not stored.
# =============================================================================

library(DBI)
library(duckdb)
library(dplyr)

DERIVED_TABLES <- c("subjects", "trials", "epochs", "alignment_qc",
                    "spectral_metrics", "source_metrics",
                    "spectral_subject_metrics", "source_subject_metrics",
                    "ingest_log")

EEG_TABLES <- c("epochs", "spectral_metrics", "source_metrics",
                "spectral_subject_metrics", "source_subject_metrics")

connect_or_create_data_db <- function(db_path) {
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
  dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = FALSE)
}

# Never touches subject_split.
drop_derived_tables <- function(con) {
  for (tb in DERIVED_TABLES) dbExecute(con, sprintf("DROP TABLE IF EXISTS %s", tb))
  invisible(NULL)
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
    )")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS trials (
      global_subjid VARCHAR,
      trial         INTEGER,
      trial_index   INTEGER,
      laser_power   DOUBLE,
      pain_rating   DOUBLE,
      PRIMARY KEY (global_subjid, trial)
    )")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS epochs (
      global_subjid VARCHAR,
      trial         INTEGER,
      event_code    VARCHAR,
      condition     VARCHAR,
      intensity     VARCHAR,
      PRIMARY KEY (global_subjid, trial)
    )")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS alignment_qc (
      global_subjid       VARCHAR PRIMARY KEY,
      experiment_name     VARCHAR,
      subjid              INTEGER,
      n_behav             INTEGER,
      n_epochs            INTEGER,
      n_spectral_trials   INTEGER,
      n_source_trials     INTEGER,
      intensity_check     VARCHAR,
      intensity_agreement DOUBLE,
      status              VARCHAR,
      reason              VARCHAR,
      ingested            BOOLEAN,
      checked_at          VARCHAR
    )")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS source_metrics (
      global_subjid VARCHAR,
      trial         INTEGER,
      roi           VARCHAR,
      metric_name   VARCHAR,
      metric_value  DOUBLE,
      PRIMARY KEY (global_subjid, trial, roi, metric_name)
    )")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS spectral_metrics (
      global_subjid VARCHAR,
      trial         INTEGER,
      channel       VARCHAR,
      metric_name   VARCHAR,
      metric_value  DOUBLE,
      PRIMARY KEY (global_subjid, trial, channel, metric_name)
    )")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS source_subject_metrics (
      global_subjid VARCHAR,
      roi           VARCHAR,
      metric_name   VARCHAR,
      metric_value  DOUBLE,
      PRIMARY KEY (global_subjid, roi, metric_name)
    )")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS spectral_subject_metrics (
      global_subjid VARCHAR,
      channel       VARCHAR,
      metric_name   VARCHAR,
      metric_value  DOUBLE,
      PRIMARY KEY (global_subjid, channel, metric_name)
    )")

  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS ingest_log (
      global_subjid VARCHAR,
      kind          VARCHAR,
      path          VARCHAR,
      n_rows        BIGINT,
      feature_set   VARCHAR,
      ingested_at   VARCHAR
    )")

  invisible(NULL)
}

# =============================================================================
# Delete-then-insert helpers
# =============================================================================
delete_by_subject <- function(con, table_name, ids) {
  if (length(ids) == 0) return(invisible(NULL))
  placeholders <- paste(rep("?", length(ids)), collapse = ",")
  dbExecute(con, sprintf("DELETE FROM %s WHERE global_subjid IN (%s)", table_name, placeholders),
            params = as.list(ids))
  invisible(NULL)
}

upsert_by_subject <- function(con, table_name, df) {
  if (nrow(df) == 0) return(invisible(NULL))
  delete_by_subject(con, table_name, unique(df$global_subjid))
  dbAppendTable(con, table_name, as.data.frame(df))
  invisible(NULL)
}

upsert_subjects <- function(con, df) upsert_by_subject(con, "subjects", df)
upsert_trials   <- function(con, df) upsert_by_subject(con, "trials", df)

# Remove every EEG-derived row for these subjects (used before re-ingesting
# a subject, and to purge subjects that fail alignment when rebuild is off).
purge_subject_eeg <- function(con, ids) {
  for (tb in EEG_TABLES) delete_by_subject(con, tb, ids)
  invisible(NULL)
}

# =============================================================================
# SQL helpers for direct CSV reads
# =============================================================================
sql_str <- function(con, x) as.character(dbQuoteString(con, x))
sql_id  <- function(con, x) as.character(dbQuoteIdentifier(con, x))

csv_columns <- function(con, path) {
  names(dbGetQuery(con, sprintf(
    "SELECT * FROM read_csv(%s, all_varchar = true, header = true) LIMIT 0", sql_str(con, path))))
}

# Number of distinct trials (and the max trial index) in a trial-level CSV.
csv_trial_count <- function(con, path) {
  r <- dbGetQuery(con, sprintf(
    "SELECT COUNT(DISTINCT TRY_CAST(TRY_CAST(trial AS DOUBLE) AS INTEGER)) AS n,
            MAX(TRY_CAST(TRY_CAST(trial AS DOUBLE) AS INTEGER)) AS max_trial
     FROM read_csv(%s, all_varchar = true, header = true)", sql_str(con, path)))
  list(n = as.integer(r$n), max_trial = as.integer(r$max_trial))
}

# =============================================================================
# ingest_feature_csv()
# -----------------------------------------------------------------------------
# table     : destination long table
# path      : per-subject CSV
# gsid      : global_subjid to stamp on every row
# key_map   : named character vector, destination column = CSV column, e.g.
#             c(trial = "trial", channel = "chan_label") or c(roi = "roi")
# drop_cols : CSV columns that are neither keys nor features
#             (subject / subjid / chan_idx)
#
# Every remaining column is a feature: metric_name = lower(column name),
# metric_value = TRY_CAST(value AS DOUBLE). Non-finite and empty values are
# skipped. Returns the number of rows inserted.
# =============================================================================
ingest_feature_csv <- function(con, table, path, gsid, key_map,
                               drop_cols = c("subject", "subjid", "chan_idx")) {
  cols <- csv_columns(con, path)
  miss <- setdiff(unname(key_map), cols)
  if (length(miss) > 0) stop("Missing key column(s) ", paste(miss, collapse = ", "), " in ", path)

  feat_cols <- setdiff(cols, c(unname(key_map), drop_cols))
  if (length(feat_cols) == 0) stop("No feature columns in ", path)

  key_src <- vapply(unname(key_map), function(k) sql_id(con, k), character(1))
  key_sel <- vapply(seq_along(key_map), function(i) {
    dst <- names(key_map)[i]
    src <- sql_id(con, key_map[[i]])
    if (dst == "trial") sprintf("TRY_CAST(TRY_CAST(%s AS DOUBLE) AS INTEGER)", src) else src
  }, character(1))
  feat_q <- paste(vapply(feat_cols, function(k) sql_id(con, k), character(1)), collapse = ", ")

  sql <- sprintf("
    INSERT INTO %s (global_subjid, %s, metric_name, metric_value)
    SELECT %s, %s, lower(metric_name), v
    FROM (
      SELECT *, TRY_CAST(metric_value AS DOUBLE) AS v
      FROM (
        UNPIVOT (SELECT %s, %s FROM read_csv(%s, all_varchar = true, header = true))
        ON %s
        INTO NAME metric_name VALUE metric_value
      )
    )
    WHERE v IS NOT NULL AND isfinite(v)",
    table,
    paste(names(key_map), collapse = ", "),
    sql_str(con, gsid), paste(key_sel, collapse = ", "),
    paste(key_src, collapse = ", "), feat_q, sql_str(con, path),
    feat_q)

  dbExecute(con, sql)
}

log_ingest <- function(con, gsid, kind, path, n_rows, feature_set) {
  dbAppendTable(con, "ingest_log", data.frame(
    global_subjid = gsid, kind = kind, path = path, n_rows = as.numeric(n_rows),
    feature_set = feature_set, ingested_at = as.character(Sys.time()),
    stringsAsFactors = FALSE))
  invisible(NULL)
}