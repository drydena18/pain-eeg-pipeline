# =============================================================================
# lateral_helpers.R                                                   V1.0.0
# -----------------------------------------------------------------------------
# Bilateral-stimulation experiments (both hands stimulated, hand per trial
# unknown: 95ByBP, 29ByBP). Spectral labels lateral clusters relative to a
# right-hand stimulus (S1_contra = left hemisphere), which is wrong for about
# half of these trials. merge_core.R calls drop_lateral_units() inside each
# such subject's ingest transaction, after the feature CSVs are loaded.
#
# Only channel-level rows whose label matches `pattern` (default: names
# ending in _contra / _ipsi) are removed. Single electrodes (C3, C4, ...) and
# non-lateral clusters (SM_furman) are kept. Source ROIs are untouched: their
# _lh / _rh labels are anatomical, not relative to the stimulated side.
# =============================================================================

drop_lateral_units <- function(con, gsid,
                               tables = c("spectral_metrics", "spectral_subject_metrics"),
                               pattern = "_(contra|ipsi)$") {
  n <- 0
  for (tb in tables) {
    if (!DBI::dbExistsTable(con, tb)) next
    ch <- DBI::dbGetQuery(con, sprintf("SELECT DISTINCT channel FROM %s WHERE global_subjid = ?", tb),
                          params = list(gsid))$channel
    lat <- ch[grepl(pattern, ch)]
    if (length(lat) == 0) next
    ph <- paste(rep("?", length(lat)), collapse = ", ")
    n <- n + DBI::dbExecute(con,
      sprintf("DELETE FROM %s WHERE global_subjid = ? AND channel IN (%s)", tb, ph),
      params = c(list(gsid), as.list(lat)))
  }
  n
}