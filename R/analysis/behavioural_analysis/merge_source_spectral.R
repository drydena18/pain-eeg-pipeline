# =============================================================================
# merge_source_spectral.R
# -----------------------------------------------------------------------------
# Merges per-subject, per-trial source-space CSVs (produced by source_core.py)
# with the behavioural + demographic master and writes source_pain_master.csv.
#
# Mirrors merge_spectral_behaviour.R in structure and upsert logic.
#
# Input files expected (one per subject, per experiment):
#   <source_root>/<exp_out>/source/sub-XXX/csv/sub-XXX_source_trial.csv
#
# Key columns in source trial CSVs:
#   subject, trial, roi,
#   pow_slow, pow_fast, pow_alpha,
#   BI_pre, LR_pre, CoG_pre, psi_cog,
#   slow_phase, sin_phase, cos_phase,
#   pow_slow_post, pow_fast_post, pow_alpha_post,
#   ERD_slow, ERD_fast, delta_ERD,
#   slow_phase_post, sin_phase_post, cos_phase_post,
#   n2_amp, n2_lat_ms, p2_amp, p2_lat_ms, n2p2_amp, n2_mean, p2_mean
#
# Output:
#   source_pain_master.csv   — one row per (experiment, subject, trial, roi)
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# =============================================================================
# USER SETTINGS
# =============================================================================
source_root  <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/"
behav_file   <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/behavioural_demo_master.csv"
output_file  <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/source_pain_master.csv"
source_pattern <- "_source_trial\\.csv$"

# =============================================================================
# EXPERIMENT LOOKUP
# =============================================================================
experiment_lookup <- tibble::tribble(
  ~experiment_name,   ~experiment_id,
  "26ByBiosemi",      1L,
  "29ByANT",          2L,
  "39ByBP",           3L,
  "30ByANT",          4L,
  "65ByANT",          5L,
  "95ByBP",           6L,
  "142ByBiosemi",     7L,
  "223ByBP",          8L,
  "29ByBP",           9L
)

# =============================================================================
# HELPERS  (mirrored from merge_spectral_behaviour.R)
# =============================================================================
clean_names_local <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\^2",  "r2") %>%
    str_replace_all("\\.",   "_")
}

extract_experiment_name <- function(file_path) {
  # Path pattern: .../da-analysis/<exp_name>/source/sub-XXX/csv/...
  exp_name <- str_match(file_path, "da-analysis/([^/]+)/source/")[, 2]
  if (is.na(exp_name)) {
    stop("Could not extract experiment_name from path: ", file_path)
  }
  exp_name
}

# Source columns to pass through alongside the join keys.
# Any column not in this list that exists in the CSV is still kept —
# this list is used only for validation.
.required_source_cols <- c("subject", "trial", "roi")

read_one_source <- function(file_path, experiment_lookup) {
  message("Reading source: ", basename(file_path))
  
  df <- read_csv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))
  
  # Normalise the subject ID column name
  if ("subject" %in% names(df) && !"subjid" %in% names(df)) {
    df <- rename(df, subjid = subject)
  }
  
  missing <- setdiff(.required_source_cols, c(names(df), "subjid"))
  if ("subjid" %in% names(df)) missing <- setdiff(missing, "subject")
  if (length(missing) > 0) {
    stop("Source trial CSV missing required columns (", paste(missing, collapse = ", "),
         "): ", file_path)
  }
  
  experiment_name <- extract_experiment_name(file_path)
  experiment_id   <- experiment_lookup %>%
    filter(experiment_name == !!experiment_name) %>%
    pull(experiment_id)
  
  if (length(experiment_id) != 1L) {
    stop("Experiment name not found in lookup table: ", experiment_name)
  }
  
  df %>%
    mutate(
      experiment_name = experiment_name,
      experiment_id   = experiment_id,
      subjid          = as.integer(subjid),
      trial           = as.integer(trial),
      roi             = as.character(roi),
      subjid_uid      = sprintf("E%02d_S%03d", experiment_id, subjid)
    )
}

merge_one_source_subject <- function(src_df, behav_master) {
  this_subjid <- unique(src_df$subjid)
  this_exp_id <- unique(src_df$experiment_id)
  this_exp_nm <- unique(src_df$experiment_name)
  
  if (length(this_subjid) != 1L) {
    warning("Source file contains multiple subject IDs; skipping.")
    return(NULL)
  }
  if (length(this_exp_id) != 1L || length(this_exp_nm) != 1L) {
    warning("Source file has ambiguous experiment ID/name; skipping.")
    return(NULL)
  }
  
  behav_sub <- behav_master %>%
    filter(experiment_id == this_exp_id, subjid == this_subjid)
  
  if (nrow(behav_sub) == 0L) {
    warning("No behavioural rows for subjid ", this_subjid, " in ", this_exp_nm)
    return(NULL)
  }
  
  # Behavioural data has one row per trial. Source data has one row per
  # (trial, roi). Join on (experiment_name, experiment_id, subjid,
  # subjid_uid, trial) — roi is source-only.
  merged <- src_df %>%
    left_join(
      behav_sub,
      by = c("experiment_name", "experiment_id", "subjid", "subjid_uid", "trial")
    )
  
  id_cols    <- c("experiment_name", "experiment_id", "subjid", "subjid_uid", "roi")
  behav_cols <- c("global_subjid", "age", "sex", "cap_size",
                  "trial", "trial_index", "laser_power", "pain_rating")
  src_only   <- setdiff(names(src_df),
                        c("experiment_name", "experiment_id",
                          "subjid", "subjid_uid", "trial", "roi"))
  
  merged %>%
    select(all_of(id_cols), all_of(behav_cols), all_of(src_only))
}

# =============================================================================
# LOAD BEHAVIOURAL MASTER
# =============================================================================
if (!file.exists(behav_file)) {
  stop("Behavioural demo master not found: ", behav_file)
}

behav_master <- read_csv(behav_file, show_col_types = FALSE)
names(behav_master) <- clean_names_local(names(behav_master))

required_behav_cols <- c(
  "experiment_name", "experiment_id", "subjid", "subjid_uid",
  "global_subjid", "age", "sex", "cap_size", "trial", "trial_index",
  "laser_power", "pain_rating"
)
missing_behav <- setdiff(required_behav_cols, names(behav_master))
if (length(missing_behav) > 0) {
  stop("Behavioural demo master missing: ", paste(missing_behav, collapse = ", "))
}

behav_master <- behav_master %>%
  mutate(
    experiment_name = as.character(experiment_name),
    experiment_id   = as.integer(experiment_id),
    subjid          = as.integer(subjid),
    subjid_uid      = as.character(subjid_uid),
    global_subjid   = as.integer(global_subjid),
    age             = suppressWarnings(as.numeric(age)),
    sex             = as.character(sex),
    cap_size        = as.factor(cap_size),
    trial           = as.integer(trial),
    trial_index     = as.integer(trial_index)
  )

# =============================================================================
# FIND SOURCE TRIAL CSVs
# =============================================================================
source_files <- list.files(
  source_root,
  pattern   = source_pattern,
  recursive = TRUE,
  full.names = TRUE
)

if (length(source_files) == 0L) {
  stop("No source trial CSVs found under: ", source_root)
}
message("Found ", length(source_files), " source trial files.")

read_one_fooof_ga <- function(csv_path, experiment_id, subjid, subjid_uid) {
  if (!file.exists(csv_path)) return(NULL)
  df <- tryCatch(read_csv(csv_path, show_col_types = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) return(NULL)
  names(df) <- clean_names_local(names(df))
  
  # Rename roi_idx -> roi if needed (older output format)
  if ("roi_idx" %in% names(df) && !"roi" %in% names(df)) {
    df <- rename(df, roi = roi_idx)
  }
  if (!"roi" %in% names(df)) return(NULL)
  
  # Keep only FOOOF metric columns + roi
  fooof_cols <- intersect(
    c("roi", "fooof_offset", "fooof_exponent", "fooof_knee",
      "fooof_alpha_cf", "fooof_alpha_pw", "fooof_alpha_bw"),
    names(df)
  )
  if (length(fooof_cols) < 2) return(NULL)
  
  df %>%
    select(all_of(fooof_cols)) %>%
    mutate(
      experiment_id = experiment_id,
      subjid        = subjid,
      subjid_uid    = subjid_uid,
      roi           = as.character(roi)
    )
}

# =============================================================================
# READ + MERGE (now includes GA FOOOF broadcast)
# =============================================================================
merged_list <- source_files %>%
  map(~tryCatch(
    read_one_source(.x, experiment_lookup),
    error = function(e) { warning(.x, ": ", conditionMessage(e)); NULL }
  )) %>%
  compact() %>%
  map(~tryCatch(
    merge_one_source_subject(.x, behav_master),
    error = function(e) { warning(conditionMessage(e)); NULL }
  )) %>%
  compact()

if (length(merged_list) == 0L) {
  stop("No source subject files were successfully merged.")
}

new_data <- bind_rows(merged_list) %>%
  mutate(
    experiment_name = as.character(experiment_name),
    experiment_id   = as.integer(experiment_id),
    subjid          = as.integer(subjid),
    subjid_uid      = as.character(subjid_uid),
    global_subjid   = as.integer(global_subjid),
    age             = suppressWarnings(as.numeric(age)),
    sex             = as.character(sex),
    cap_size        = as.factor(cap_size),
    trial           = as.integer(trial),
    trial_index     = as.integer(trial_index),
    roi             = as.character(roi)
  )

# ── Broadcast GA FOOOF metrics to every trial row ────────────────────────────
# FOOOF is computed at the grand-average level (one value per subject × ROI).
# We join it here so models m12-m14 in run_gamm_source.R have access to
# fooof_offset / fooof_exponent on every trial row, enabling FOOOF-controlled
# GAMMs at the trial level.
message("Broadcasting GA FOOOF metrics to trial rows...")
fooof_ga_list <- list()

for (i in seq_len(nrow(merged_list |> bind_rows() |> distinct(experiment_id, subjid, subjid_uid)))) {
  # iterate via the actual unique subjects in new_data
}
# Cleaner: build FOOOF GA from the file system matching the source trial CSVs
fooof_sources <- source_files %>%
  map_dfr(function(fp) {
    fooof_path <- str_replace(fp, "_source_trial\\.csv$", "_source_ga_fooof.csv")
    exp_name   <- extract_experiment_name(fp)
    exp_id_val <- experiment_lookup %>%
      filter(experiment_name == exp_name) %>% pull(experiment_id)
    if (length(exp_id_val) != 1L) return(NULL)
    
    sub_str    <- str_extract(basename(fp), "sub-\\d+")
    sub_id     <- as.integer(str_remove(sub_str, "sub-"))
    sub_uid    <- sprintf("E%02d_S%03d", exp_id_val, sub_id)
    
    read_one_fooof_ga(fooof_path, exp_id_val, sub_id, sub_uid)
  }) %>%
  bind_rows()

if (nrow(fooof_sources) > 0) {
  new_data <- new_data %>%
    left_join(fooof_sources,
              by = c("experiment_id", "subjid", "subjid_uid", "roi"))
  message("  FOOOF columns joined: ",
          paste(intersect(names(fooof_sources),
                          c("fooof_offset","fooof_exponent","fooof_knee",
                            "fooof_alpha_cf","fooof_alpha_pw","fooof_alpha_bw")),
                collapse = ", "))
} else {
  message("  No GA FOOOF files found — fooof_offset/exponent will be absent.")
}

# =============================================================================
# UPSERT INTO EXISTING MASTER  (same anti-join pattern as merge_spectral)
# Key is (experiment_id, subjid, trial, roi) — unique row identity.
# =============================================================================
if (file.exists(output_file)) {
  old_data <- read_csv(output_file, show_col_types = FALSE)
  names(old_data) <- clean_names_local(names(old_data))
  
  old_data <- old_data %>%
    mutate(
      experiment_name = as.character(experiment_name),
      experiment_id   = as.integer(experiment_id),
      subjid          = as.integer(subjid),
      subjid_uid      = as.character(subjid_uid),
      global_subjid   = as.integer(global_subjid),
      age             = suppressWarnings(as.numeric(age)),
      sex             = as.character(sex),
      cap_size        = as.factor(cap_size),
      trial           = as.integer(trial),
      trial_index     = as.integer(trial_index),
      roi             = as.character(roi)
    )
  
  all_cols <- union(names(old_data), names(new_data))
  for (col in setdiff(all_cols, names(old_data))) old_data[[col]] <- NA
  for (col in setdiff(all_cols, names(new_data))) new_data[[col]] <- NA
  
  old_data <- select(old_data, all_of(all_cols))
  new_data <- select(new_data, all_of(all_cols))
  
  updated_master <- old_data %>%
    anti_join(
      new_data %>% select(experiment_id, subjid, trial, roi),
      by = c("experiment_id", "subjid", "trial", "roi")
    ) %>%
    bind_rows(new_data)
} else {
  updated_master <- new_data
}

updated_master <- updated_master %>%
  arrange(experiment_id, subjid, trial, roi)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(updated_master, output_file)

message("Saved source_pain_master: ", output_file)
message("  Rows     : ", nrow(updated_master))
message("  Subjects : ", n_distinct(updated_master$subjid_uid))
message("  ROIs     : ", n_distinct(updated_master$roi), " — ",
        paste(sort(unique(updated_master$roi)), collapse = ", "))
message("  Experiments: ", n_distinct(updated_master$experiment_id))