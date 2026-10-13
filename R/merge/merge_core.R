# =============================================================================
# merge_core.R
# =============================================================================
# Builds data.duckdb from the per-subject pipeline outputs. Run from R/merge/.
#
# Tables (see helpers/db_write_helpers.R):
#   subjects, trials                      behaviour + participants + cap size
#   epochs                                07_epoch metadata (event code, low/high)
#   alignment_qc                          one row per subject: trial-alignment status
#   spectral_metrics, source_metrics      trial x channel / ROI features (long)
#   spectral_subject_metrics,
#   source_subject_metrics                subject x channel / ROI features (long:
#                                         TVI, ITC, FOOOF, trial-mean summaries)
#   ingest_log                            file -> rows provenance
#   subject_split                         train / holdout (persistent, never dropped)
#
# V2 (Oct 2026), feature set v3.1:
#   - Fixes: subejcts_df typo (crash), cap_size never joined (crash),
#     `x.` / `messsage` typos in the error paths, per-file errors silently
#     downgraded to warnings.
#   - Trial alignment check before any EEG row is written (alignment_helpers.R).
#   - Feature CSVs are unpivoted inside DuckDB, one subject at a time, in a
#     transaction (a failure leaves no partial subject behind).
#   - Subject-level summaries are ingested.
#   - Subjects with a rating outside merge_cfg$rating_range (default 0-10)
#     get no EEG rows (status "rating_out_of_range" in alignment_qc).
#   - rebuild_derived_tables: database mirrors the files on disk each run;
#     subject_split is preserved.
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

source("merge_defaults.R")
source("helpers/read_helpers.R")
source("helpers/db_write_helpers.R")
source("helpers/alignment_helpers.R")
source("helpers/split_helpers.R")

# =============================================================================
# CONNECT + SCHEMA
# =============================================================================
exp_lookup <- merge_cfg$experiment_lookup
if (is.null(exp_lookup) || nrow(exp_lookup) == 0) stop("merge_cfg$experiment_lookup is empty.")

con <- connect_or_create_data_db(merge_cfg$data_db_path)
# Closed explicitly at the end (on.exit() has no effect at the top level of a
# sourced script). If the run stops early, close it before re-running in the
# same R session: dbDisconnect(con, shutdown = TRUE)

if (isTRUE(merge_cfg$rebuild_derived_tables)) {
  message("Rebuild: dropping derived tables (subject_split is kept).")
  drop_derived_tables(con)
}
init_data_schema(con)

# =============================================================================
# STEP 1: BEHAVIOUR
# =============================================================================
message("== Step 1: Behaviour ==")
behav_files <- list.files(merge_cfg$behav_dir, pattern = merge_cfg$behav_pattern, full.names = TRUE)
if (length(behav_files) == 0) stop("No behavioural CSVs found in: ", merge_cfg$behav_dir)

behav_errors <- character(0)
behaviour_all <- map(behav_files, function(f) {
  tryCatch(read_behavioural_file(f, exp_lookup), error = function(e) {
    behav_errors <<- c(behav_errors, paste0(basename(f), ": ", conditionMessage(e)))
    NULL
  })
}) %>% compact() %>% bind_rows()

if (length(behav_errors) > 0) message("  [ERROR] ", paste(behav_errors, collapse = "\n  [ERROR] "))
if (nrow(behaviour_all) == 0) stop("No behavioural rows could be read.")

behaviour_all <- behaviour_all %>%
  arrange(experiment_id, subjid, trial) %>%
  group_by(experiment_id, subjid) %>%
  mutate(trial_index = row_number()) %>%
  ungroup() %>%
  mutate(global_subjid = sprintf("E%02d_S%03d", experiment_id, subjid))

message("  Rows: ", nrow(behaviour_all), "  Subjects: ", n_distinct(behaviour_all$global_subjid))

# Subjects with any rating outside merge_cfg$rating_range are excluded from
# every EEG table (behaviour and split are kept, as for alignment failures).
rr <- merge_cfg$rating_range
bad_rating <- behaviour_all %>%
  filter(!is.na(pain_rating), pain_rating < rr[1] | pain_rating > rr[2]) %>%
  group_by(global_subjid) %>%
  summarise(reason = sprintf("%d rating(s) outside %g-%g (e.g. %s)", n(), rr[1], rr[2],
                             paste(head(signif(unique(pain_rating), 3), 3), collapse = ", ")),
            .groups = "drop")
if (nrow(bad_rating) > 0) {
  message("  Excluding ", nrow(bad_rating), " subject(s) with ratings outside ", rr[1], "-", rr[2],
          if (nrow(bad_rating) > 20) " (first 20 shown):" else ":")
  for (k in seq_len(min(20, nrow(bad_rating))))
    message("    ", bad_rating$global_subjid[k], "  ", bad_rating$reason[k])
}

# =============================================================================
# STEP 2: SUBJECTS (participants + cap size) + TRIALS
# =============================================================================
message("\n== Step 2: Subjects + trials ==")

participants_all <- pmap(exp_lookup, function(experiment_name, experiment_id) {
  fp <- file.path(merge_cfg$data_root, experiment_name, "resource", merge_cfg$participants_filename)
  if (!file.exists(fp)) {
    if (experiment_name %in% behaviour_all$experiment_name)
      message("  [WARN] participants.tsv not found for: ", experiment_name)
    return(NULL)
  }
  tryCatch(read_participants_file(fp, experiment_name, experiment_id),
           error = function(e) { message("  [ERROR] ", experiment_name, ": ", conditionMessage(e)); NULL })
}) %>% compact() %>% bind_rows()
if (nrow(participants_all) == 0) {
  participants_all <- tibble(subjid = integer(), age = double(), sex = character(),
                             experiment_name = character(), experiment_id = integer())
}

if (!file.exists(merge_cfg$cap_size_file)) stop("cap_size file not found: ", merge_cfg$cap_size_file)
cap_df <- read_cap_size_file(merge_cfg$cap_size_file)

subject_key <- behaviour_all %>%
  distinct(global_subjid, experiment_id, experiment_name, subjid)

subjects_df <- subject_key %>%
  left_join(participants_all, by = c("experiment_name", "experiment_id", "subjid")) %>%
  left_join(cap_df, by = c("experiment_name", "experiment_id")) %>%
  transmute(global_subjid, subjid, subjid_uid = global_subjid, experiment_id,
            experiment_name, age, sex, cap_size)

no_cap <- unique(subjects_df$experiment_name[is.na(subjects_df$cap_size)])
if (length(no_cap) > 0) stop("No cap_size for experiment(s): ", paste(no_cap, collapse = ", "),
                             " — add them to ", merge_cfg$cap_size_file)

message("  Subjects: ", nrow(subjects_df),
        "  missing age: ", sum(is.na(subjects_df$age)),
        "  missing sex: ", sum(is.na(subjects_df$sex)))

trials_df <- behaviour_all %>% transmute(global_subjid, trial, trial_index, laser_power, pain_rating)

upsert_subjects(con, subjects_df)
upsert_trials(con, trials_df)

# =============================================================================
# STEP 3: DISCOVER PIPELINE OUTPUTS
# =============================================================================
message("\n== Step 3: Pipeline outputs ==")
files_df <- discover_pipeline_files(merge_cfg$data_root, exp_lookup, merge_cfg$files)
if (nrow(files_df) > 0) {
  tab <- files_df %>% count(experiment_name, kind) %>%
    tidyr::pivot_wider(names_from = kind, values_from = n, values_fill = 0)
  print(as.data.frame(tab), row.names = FALSE)
} else {
  message("  [WARN] No pipeline outputs found under ", merge_cfg$data_root)
}
files_df <- files_df %>% mutate(global_subjid = sprintf("E%02d_S%03d", experiment_id, subjid))

# =============================================================================
# STEP 4: PER-SUBJECT ALIGNMENT + INGEST
# =============================================================================
message("\n== Step 4: Alignment check + feature ingest ==")

all_subjects <- bind_rows(
  subject_key %>% select(global_subjid, experiment_name, subjid),
  files_df    %>% distinct(global_subjid, experiment_name, subjid)
) %>% distinct(global_subjid, .keep_all = TRUE) %>% arrange(global_subjid)

ingest_spec <- list(
  spectral_trial   = list(table = "spectral_metrics",         key = c(trial = "trial", channel = "chan_label")),
  source_trial     = list(table = "source_metrics",           key = c(trial = "trial", roi = "roi")),
  spectral_subject = list(table = "spectral_subject_metrics", key = c(channel = "chan_label")),
  source_subject   = list(table = "source_subject_metrics",   key = c(roi = "roi"))
)

qc_rows <- vector("list", nrow(all_subjects))
n_ingested_rows <- 0

for (i in seq_len(nrow(all_subjects))) {
  gsid <- all_subjects$global_subjid[i]
  sf   <- files_df %>% filter(global_subjid == gsid)
  path_of <- function(kind) { p <- sf$path[sf$kind == kind]; if (length(p)) p[1] else NA_character_ }

  behav  <- behaviour_all %>% filter(global_subjid == gsid) %>% select(trial, laser_power)
  epochs <- NULL
  n_spec <- NA_integer_
  n_src  <- NA_integer_
  read_err <- NULL

  tryCatch({
    if (!is.na(path_of("epoch_meta")))     epochs <- read_epoch_metadata(path_of("epoch_meta"))
    if (!is.na(path_of("spectral_trial"))) n_spec <- csv_trial_count(con, path_of("spectral_trial"))$n
    if (!is.na(path_of("source_trial")))   n_src  <- csv_trial_count(con, path_of("source_trial"))$n
  }, error = function(e) read_err <<- conditionMessage(e))

  a <- if (is.null(read_err)) assess_alignment(behav, epochs, n_spec, n_src) else
    list(n_behav = nrow(behav), n_epochs = NA_integer_, intensity_check = "not_run",
         intensity_agreement = NA_real_, status = "read_error", reason = read_err)
  if (gsid %in% bad_rating$global_subjid && a$status != "no_eeg") {
    a$status <- "rating_out_of_range"
    a$reason <- bad_rating$reason[bad_rating$global_subjid == gsid]
  }

  do_ingest <- a$status == "ok" ||
    (!isTRUE(merge_cfg$exclude_misaligned) && a$status %in% c("count_mismatch", "intensity_mismatch"))
  ingested <- FALSE

  if (!isTRUE(merge_cfg$rebuild_derived_tables)) purge_subject_eeg(con, gsid)

  if (do_ingest) {
    res <- tryCatch({
      dbBegin(con)
      n_sub <- 0
      if (!is.null(epochs)) {
        dbAppendTable(con, "epochs", as.data.frame(epochs %>% mutate(global_subjid = gsid, .before = 1)))
        log_ingest(con, gsid, "epoch_meta", path_of("epoch_meta"), nrow(epochs), merge_cfg$feature_set)
      }
      for (kind in names(ingest_spec)) {
        p <- path_of(kind)
        if (is.na(p)) next
        n <- ingest_feature_csv(con, ingest_spec[[kind]]$table, p, gsid, ingest_spec[[kind]]$key)
        if (n == 0) stop("0 feature rows ingested from ", p)
        log_ingest(con, gsid, kind, p, n, merge_cfg$feature_set)
        n_sub <- n_sub + n
      }
      dbCommit(con)
      n_sub
    }, error = function(e) {
      try(dbRollback(con), silent = TRUE)
      structure(conditionMessage(e), class = "ingest_error")
    })

    if (inherits(res, "ingest_error")) {
      a$status <- "ingest_error"
      a$reason <- as.character(res)
    } else {
      ingested <- TRUE
      n_ingested_rows <- n_ingested_rows + res
    }
  }

  qc_rows[[i]] <- tibble(
    global_subjid = gsid, experiment_name = all_subjects$experiment_name[i],
    subjid = all_subjects$subjid[i], n_behav = as.integer(a$n_behav),
    n_epochs = as.integer(a$n_epochs), n_spectral_trials = as.integer(n_spec),
    n_source_trials = as.integer(n_src), intensity_check = a$intensity_check,
    intensity_agreement = a$intensity_agreement, status = a$status, reason = a$reason,
    ingested = ingested, checked_at = as.character(Sys.time()))
}

qc_df <- bind_rows(qc_rows)
upsert_by_subject(con, "alignment_qc", qc_df)

# -- Summary ----
message("\n  Alignment status by experiment:")
st <- qc_df %>% count(experiment_name, status) %>%
  tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0)
print(as.data.frame(st), row.names = FALSE)

problems <- qc_df %>% filter(!status %in% c("ok", "no_eeg"))
if (nrow(problems) > 0) {
  message("\n  Subjects NOT ingested (first 20; full list: SELECT * FROM alignment_qc WHERE NOT ingested):")
  for (k in seq_len(min(20, nrow(problems)))) {
    message(sprintf("    %s  %-24s %s", problems$global_subjid[k], problems$status[k], problems$reason[k]))
  }
}
message("  Feature rows written: ", format(n_ingested_rows, big.mark = ","))

has_eeg_files <- nrow(files_df) > 0
if (has_eeg_files && sum(qc_df$ingested) == 0) {
  stop("EEG outputs were found but no subject was ingested — see alignment_qc.")
}

# =============================================================================
# STEP 5: TRAIN/HOLDOUT SPLIT — incremental, never reassigns existing subjects
# =============================================================================
message("\n== Step 5: Train/holdout split assignment ==")
assign_new_subjects_to_split(con, train_frac = merge_cfg$split_train_frac,
                             seed = merge_cfg$split_seed)

dbDisconnect(con, shutdown = TRUE)
message("\nMerge complete. Database: ", merge_cfg$data_db_path)