# =============================================================================
# read_helpers.R
# -----------------------------------------------------------------------------
# Per-file-type readers for the merge stage. Each function reads ONE file
# and returns a tidy data frame in a known shape — no cross-file joining
# here, that happens in merge_core.R once everything's been read.
#
# Column-name harmonization notes (confirmed against real sample files):
#   - Behavioural CSVs: ID/Trial_Num/Painrating -> subjid/trial/pain_rating
#     (verified against all 9 real behavioural CSVs)
#   - participants.tsv: Age/Gender already match expected names directly
#   - spectral_chan_by_trial.csv: already lowercase metric names
#     (bi_pre, cog_pre, lr_pre, delta_erd, ...) — canonical form
#   - source_trial.csv: mixed-case metric names (BI_pre, CoG_pre, LR_pre,
#     delta_ERD) — MUST be lowercased at ingestion to match spectral's
#     convention, or the same concept becomes two different metric_name
#     values in the long-format table (verified this was a real risk with
#     real sample files: BI_pre vs bi_pre, CoG_pre vs cog_pre, etc.)
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

clean_names_local <- function(x) {
  x %>% str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\^2",  "r2") %>%
    str_replace_all("\\.",   "_")
}

rename_cols <- function(df, mapping) {
  cur <- names(df)
  for (new_nm in names(mapping)) {
    candidates <- mapping[[new_nm]]
    hit <- intersect(candidates, cur)
    if (length(hit) > 0 && !(new_nm %in% cur)) {
      df <- rename(df, !!new_nm := !!hit[1])
      cur <- names(df)
    }
  }
  df
}

harmonize_sex <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  case_when(
    x %in% c("m", "male")   ~ "M",
    x %in% c("f", "female") ~ "F",
    .default = NA_character_
  )
}

extract_experiment_name_from_path <- function(file_path, valid_names) {
  fname <- basename(file_path)
  matches <- valid_names[str_detect(fname, fixed(valid_names))]
  if (length(matches) == 1L) return(matches)
  # fall back to parent-directory matching for files under <root>/<exp_name>/...
  matches2 <- valid_names[str_detect(file_path, fixed(valid_names))]
  if (length(matches2) == 1L) return(matches2)
  if (length(matches) == 0L && length(matches2) == 0L) {
    stop("Could not match file to an experiment name: ", file_path)
  }
  stop("Multiple experiment name matches found for: ", file_path)
}

# =============================================================================
# read_behavioural_file()
# -----------------------------------------------------------------------------
# Returns: subjid, trial, laser_power, pain_rating, experiment_name,
#          experiment_id   (one row per trial; trial_index added later in
#          merge_core.R once data is sorted, same as the old script did)
# =============================================================================
read_behavioural_file <- function(file_path, experiment_lookup) {
  df <- read_csv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  df <- rename_cols(df, list(
    subjid      = c("subjid", "ID", "id", "participant"),
    trial       = c("trial", "Trial", "trial_number", "Trial_Num", "trial_num"),
    laser_power = c("laser_power", "laser", "laser_level"),
    pain_rating = c("pain_rating", "pain", "rating", "Painrating")
  ))

  required <- c("subjid", "trial", "laser_power", "pain_rating")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Behavioural file missing required columns (", paste(missing, collapse = ", "),
         "): ", file_path)
  }

  exp_name <- extract_experiment_name_from_path(file_path, experiment_lookup$experiment_name)
  exp_id   <- experiment_lookup$experiment_id[experiment_lookup$experiment_name == exp_name]

  df %>%
    transmute(
      subjid          = as.integer(subjid),
      trial           = as.integer(trial),
      laser_power     = suppressWarnings(as.numeric(laser_power)),
      pain_rating     = suppressWarnings(as.numeric(pain_rating)),
      experiment_name = exp_name,
      experiment_id   = exp_id
    )
}

# =============================================================================
# read_participants_file()
# -----------------------------------------------------------------------------
# Returns: subjid, age, sex, experiment_name, experiment_id
# =============================================================================
read_participants_file <- function(file_path, experiment_name, experiment_id) {
  df <- read_tsv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  df <- rename_cols(df, list(
    participant_id = c("participant_id", "participant"),
    age            = c("age", "Age"),
    sex            = c("sex", "Sex", "gender", "Gender")
  ))

  if (!"participant_id" %in% names(df)) {
    warning("participant_id column missing in: ", file_path)
    return(NULL)
  }
  if (!"age" %in% names(df)) df$age <- NA_real_
  if (!"sex" %in% names(df)) df$sex <- NA_character_

  df %>%
    transmute(
      subjid          = as.integer(str_remove(participant_id, "^sub-")),
      age             = suppressWarnings(as.numeric(age)),
      sex             = harmonize_sex(sex),
      experiment_name = experiment_name,
      experiment_id   = experiment_id
    )
}

# =============================================================================
# read_cap_size_file()
# -----------------------------------------------------------------------------
# Returns: experiment_name, experiment_id, cap_size
# =============================================================================
read_cap_size_file <- function(file_path) {
  df <- read_csv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  required <- c("experiment_name", "experiment_id", "cap_size")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("cap_size file missing required columns: ", paste(missing, collapse = ", "))
  }

  df %>%
    transmute(
      experiment_name = as.character(experiment_name),
      experiment_id   = as.integer(experiment_id),
      cap_size        = as.character(cap_size)
    ) %>%
    distinct(experiment_name, experiment_id, .keep_all = TRUE)
}

# =============================================================================
# read_spectral_file()
# -----------------------------------------------------------------------------
# Reads one sub-XXX_spectral_chan_by_trial.csv, pivots wide metric columns
# into long format (metric_name, metric_value), lowercases metric names
# (already lowercase in this file per the real sample, but normalized
# explicitly rather than assumed, in case a future pipeline revision
# changes casing).
#
# Returns long format: subjid, trial, channel, metric_name, metric_value,
#                       experiment_name, experiment_id
# =============================================================================
read_spectral_file <- function(file_path, experiment_lookup) {
  df <- read_csv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  if (!"subjid" %in% names(df) && "subject" %in% names(df)) {
    df <- rename(df, subjid = subject)
  }
  if (!"chan_label" %in% names(df)) {
    stop("Spectral file missing chan_label column: ", file_path)
  }

  required <- c("subjid", "trial", "chan_label")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Spectral file missing required columns (", paste(missing, collapse = ", "),
         "): ", file_path)
  }

  exp_name <- extract_experiment_name_from_path(file_path, experiment_lookup$experiment_name)
  exp_id   <- experiment_lookup$experiment_id[experiment_lookup$experiment_name == exp_name]

  id_cols <- c("subjid", "trial", "chan_idx", "chan_label")
  metric_cols <- setdiff(names(df), id_cols)

  df %>%
    select(subjid, trial, channel = chan_label, all_of(metric_cols)) %>%
    pivot_longer(
      cols      = all_of(metric_cols),
      names_to  = "metric_name",
      values_to = "metric_value"
    ) %>%
    mutate(
      metric_name = str_to_lower(metric_name),
      subjid      = as.integer(subjid),
      trial       = as.integer(trial),
      channel     = as.character(channel)
    ) %>%
    filter(!is.na(metric_value)) %>%
    mutate(experiment_name = exp_name, experiment_id = exp_id)
}

# =============================================================================
# read_source_file()
# -----------------------------------------------------------------------------
# Reads one sub-XXX_source_trial.csv, pivots wide metric columns into long
# format. CRITICAL: lowercases metric names to match spectral's convention
# (BI_pre -> bi_pre, CoG_pre -> cog_pre, etc.) — confirmed necessary against
# real sample files, where source uses mixed case and spectral does not.
#
# Returns long format: subjid, trial, roi, metric_name, metric_value,
#                       experiment_name, experiment_id
# =============================================================================
read_source_file <- function(file_path, experiment_lookup) {
  df <- read_csv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  if ("subject" %in% names(df) && !"subjid" %in% names(df)) {
    df <- rename(df, subjid = subject)
  }

  required <- c("subjid", "trial", "roi")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Source file missing required columns (", paste(missing, collapse = ", "),
         "): ", file_path)
  }

  exp_name <- extract_experiment_name_from_path(file_path, experiment_lookup$experiment_name)
  exp_id   <- experiment_lookup$experiment_id[experiment_lookup$experiment_name == exp_name]

  id_cols <- c("subjid", "trial", "roi")
  metric_cols <- setdiff(names(df), id_cols)

  df %>%
    select(subjid, trial, roi, all_of(metric_cols)) %>%
    pivot_longer(
      cols      = all_of(metric_cols),
      names_to  = "metric_name",
      values_to = "metric_value"
    ) %>%
    mutate(
      metric_name = str_to_lower(metric_name),   # <- the critical harmonization step
      subjid      = as.integer(subjid),
      trial       = as.integer(trial),
      roi         = as.character(roi),
      metric_value = suppressWarnings(as.numeric(metric_value))
    ) %>%
    filter(!is.na(metric_value)) %>%
    mutate(experiment_name = exp_name, experiment_id = exp_id)
}