# =============================================================================
# merge_core.R
# =============================================================================
# Replaces all 5 old merge_*.R scripts. Reads directly from per-subject
# pipeline output CSVs (skipping the old master-CSV intermediate layer
# entirely) and writes to data.duckdb's normalized schema:
#   subjects, trials, source_metrics (long), spectral_metrics (long)
#
# Also runs incremental train/holdout split assignment for any new subjects
# (stratified by cap_size, quota-based - see helpers/split_helpers.R).
#
# Idempotent: re-running re-ingests all discovered files (upsert per
# subject), so it's safe to run repeatedly as new subjects clear the
# pipeline. Existing subject_split assignments are never touched.
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

source("merge_defaults.R")
source("helpers/read_helpers.R")
source("helpers/db_write_helpers.R")
source("helpers/split_helpers.R")

# =============================================================================
# CONNECT + INIT SCHEMA
# =============================================================================
con <- connect_or_create_data_db(merge_cfg$data_db_path)
init_data_schema(db)

exp_lookup <- merge_cfg$experiment_lookup

# =============================================================================
# STEP 1: BEHAVIOURAL (subjid, trial, laser_power, pain_rating per experiment)
# =============================================================================
message("== Step 1: Behavioural ==")
behav_files <- list.files(merge_cfg$behav_dir, pattern = merge_cfg$behav_pattern, recursive = TRUE, full.names = TRUE)
if (length(behav_files) == 0) stop("No behavioural CSVs found in: ", merge_cfg$behav_dir)
message("Found ", length(behav_files), " behavioural files.")

behaviour_all <- map(behav_files, ~tryCatch(
    read_behavioural_file(.x, exp_lookup),
    error = function(e) { warning(.x, ": ", conditionMessage(e)); NULL }
)) %>% compact() %>% bind_rows()

# trial_index: position within (experiment_id, subjid), same derivation as
# the old merge_behavioural.R
behaviour_all <- behaviour_all %>%
    arrange(experiment_id, subjid, trial) %>%
    group_by(experiment_id, subjid) %>%
    mutate(trial_index = row_number()) %>%
    ungroup() %>%
    mutate(subjid_uid = sprintf("E%02d_S%03d", experiment_id, subjid))

message(" Behavioural rows: ", nrow(behaviour_all), " Subjects: ", n_distinct(behaviour_all$subjid_uid))

# =============================================================================
# STEP 2: PARTICIPANTS (age/sex) + CAP SIZE -> build subjects table
# =============================================================================
message("\n== Step 2: Participants + cap size ==")

participants_all <- pmap(exp_lookup, function(experiment_name, experiment_id) {
    file_path <- file.path(merge_cfg$data_root, experiment_name, "resource", merge_cfg$participants_filename)
    if (!file.exists(file_path)) {
        warning("participants.tsv not found for: ", experiment_name)
        return(NULL)
    }
    tryCatch(
        read_participants_file(file_path, experiment_name, experiment_id),
        error = function(e) { warning(experiment_name, ": ", conditionMessage(e)); NULL }
    )
}) %>% compact() %>% bind_rows()

if (!file.exists(merge_cfg$cap_size_file)) {
    stop("cap_size file not found: ", merge_cfg$cap_size_file)
}
cap_df <- read_cap_size_file(merge_cfg$cap_size_file)

# =============================================================================
# STEP 3: ASSEMBLE subjects + trials, with stable global_subjid
# =============================================================================
# global_subjid must be STABLE across re-runs - a subject ingested today
# must keep the same global_subjid tomorow, or every downstream table
# (and subject_split) silently desyncs. Strategy: global_subjid is derived
# deterministically from (experiment_id, subjid), not from row order at
# ingestion time, so re-running produces identical IDs for existing
# subjects regardless of what new subjects get added in between.
# =============================================================================
message("\n== Step 3: Assemble subjects + trials ==")

subject_key <- behaviour_all %>%
    distinct(experiment_id, experiment_name, subjid, subjid_uid) %>%
    mutate(global_subjid = subjid_uid)

subjects_df <- subject_key %>%
    left_join(participants_all, by = c("experiment_name", "experiment_id", "subjid"))%>%
    transmute(
        global_subjid, subjid, subjid_uid, experiment_id, experiment_name, age, sex, cap_size
    )

message(" Subjects assembled: ", nrow(subjects_df))
message("  Missing age: ", sum(is.na(subjects_df$age)),
        "  Missing sex: ", sum(is.na(subjects_df$sex)),
        "  Missing cap_size: ", sum(is.na(subjects_df$cap_size)))

trials_df <- behaviour_all %>%
    left_join(subject_key %>% select(experiment_id, subjid, global_subjid), by = c("experiment_id", "subjid")) %>%
    transmute(global_subjid, trial, trial_index, laser_power, pain_rating)

upsert_subjects(con, subejcts_df)
upsert_trials(con, trials_df)
message(" Written to data.duckdb: subjects, trials.")

# =============================================================================
# STEP 4: SPECTRAL METRICS (long format)
# =============================================================================
message("\n== Step 4: Spectral Metrics ==")
spectral_files <- list.files(merge_cfg$data_root, pattern = merge_cfg$spectral_pattern, recursive = TRUE, full.names = TRUE)
message("Found ", length(spectral_files), " spectral files.")

if (length(spectral_files) > 0) {
    spectral_long_list <- map(spectral_files, ~tryCatch(
        read_spectral_file(.x, exp_lookup),
        error = function(e) { warning(x., ": ", conditionMessage(e)); NULL }
    )) %>% compact()

    spectral_long <- bind_rows(spectral_long_list) %>%
        left_join(subject_key %>% select(experiment_id, subjid, global_subjid), by = c("experiment_id", "subjid")) %>%
        filter(!is.na(global_subjid)) %>%
        transmute(global_subjid, trial, channel, metric_name, metric_value)

    upsert_metrics_long(con, "spectral_metrics", spectral_long)
    message(" Spectral metric rows written: ", nrow(spectral_long))
} else {
    messsage(" [WARN] No spectral files found.")
}

# =============================================================================
# STEP 5: SOURCE METRICS (long format)
# =============================================================================
message("\n== Step 5: Source metrics ==")
source_files <- list.files(merge_cfg$data_root, pattern = merge_cfg$source_pattern,
                            recursive = TRUE, full.names = TRUE)
message("Found ", length(source_files), " source files.")
 
if (length(source_files) > 0) {
  source_long_list <- map(source_files, ~tryCatch(
    read_source_file(.x, exp_lookup),
    error = function(e) { warning(.x, ": ", conditionMessage(e)); NULL }
  )) %>% compact()
 
  source_long <- bind_rows(source_long_list) %>%
    left_join(subject_key %>% select(experiment_id, subjid, global_subjid),
              by = c("experiment_id", "subjid")) %>%
    filter(!is.na(global_subjid)) %>%
    transmute(global_subjid, trial, roi, metric_name, metric_value)
 
  upsert_metrics_long(con, "source_metrics", source_long)
  message("  Source metric rows written: ", nrow(source_long))
} else {
  message("  [WARN] No source files found.")
}
 
# =============================================================================
# STEP 6: TRAIN/HOLDOUT SPLIT — incremental, never reassigns existing subjects
# =============================================================================
message("\n== Step 6: Train/holdout split assignment ==")
assign_new_subjects_to_split(con, train_frac = merge_cfg$split_train_frac,
                              seed = merge_cfg$split_seed)
 
# =============================================================================
# DONE
# =============================================================================
dbDisconnect(con, shutdown = TRUE)
message("\nMerge complete. Database: ", merge_cfg$data_db_path)