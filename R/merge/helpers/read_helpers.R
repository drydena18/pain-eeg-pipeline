# =============================================================================
# read_helpers.R
# -----------------------------------------------------------------------------
# Readers for the small per-experiment inputs (behaviour, participants,
# cap size, epoch metadata) and discovery of per-subject pipeline files.
#
# V2 (Oct 2026):
#   - read_spectral_file() / read_source_file() removed. Trial- and
#     subject-level feature CSVs are now ingested directly by DuckDB
#     (db_write_helpers.R: ingest_feature_csv), which unpivots in SQL instead
#     of holding every subject's long table in R memory (spectral_metrics
#     alone is of the order of 10^8 rows at full N with the v3.1 set).
#   - discover_pipeline_files(): lists files per experiment folder and takes
#     the subject ID from the file name, so no experiment-name guessing
#     from paths.
#   - read_epoch_metadata(): 07_epoch/sub-XXX_epoch_metadata.csv
#     (epoch_index, event_code, condition, intensity), written by
#     preproc_core.m >= 2.2.0.
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

# =============================================================================
# read_behavioural_file()
# -----------------------------------------------------------------------------
# File name must start with the experiment name: <experiment_name>_behaviour.csv
# Returns: subjid, trial, laser_power, pain_rating, experiment_name, experiment_id
# =============================================================================
read_behavioural_file <- function(file_path, experiment_lookup) {
  exp_name <- str_remove(basename(file_path), "_behaviour\\.csv$")
  if (!exp_name %in% experiment_lookup$experiment_name) {
    stop("Behavioural file name does not match any experiment_lookup name: ", basename(file_path))
  }
  exp_id <- experiment_lookup$experiment_id[experiment_lookup$experiment_name == exp_name]

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
    stop("Behavioural file missing required columns (", paste(missing, collapse = ", "), "): ", file_path)
  }

  df %>%
    transmute(
      subjid          = as.integer(subjid),
      trial           = as.integer(trial),
      laser_power     = suppressWarnings(as.numeric(laser_power)),
      pain_rating     = suppressWarnings(as.numeric(pain_rating)),
      experiment_name = exp_name,
      experiment_id   = as.integer(exp_id)
    )
}

# =============================================================================
# read_participants_file()
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
      experiment_id   = as.integer(experiment_id)
    )
}

# =============================================================================
# read_cap_size_file()
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
# discover_pipeline_files()
# -----------------------------------------------------------------------------
# For every experiment and every file kind in merge_cfg$files, lists
# <data_root>/<experiment_name>/<subdir>/**/<pattern> and parses the subject
# ID from the leading "sub-NNN_".
#
# Returns: experiment_name, experiment_id, subjid, kind, path
# Errors if one subject has two files of the same kind (ambiguous).
# =============================================================================
discover_pipeline_files <- function(data_root, experiment_lookup, file_specs) {
  rows <- list()
  for (i in seq_len(nrow(experiment_lookup))) {
    exp_name <- experiment_lookup$experiment_name[i]
    exp_id   <- experiment_lookup$experiment_id[i]
    for (kind in names(file_specs)) {
      root <- file.path(data_root, exp_name, file_specs[[kind]]$subdir)
      if (!dir.exists(root)) next
      paths <- list.files(root, pattern = file_specs[[kind]]$pattern,
                          recursive = TRUE, full.names = TRUE)
      if (length(paths) == 0) next
      rows[[length(rows) + 1]] <- tibble(
        experiment_name = exp_name,
        experiment_id   = as.integer(exp_id),
        subjid          = as.integer(str_match(basename(paths), "^sub-(\\d+)_")[, 2]),
        kind            = kind,
        path            = paths
      )
    }
  }
  out <- bind_rows(rows)
  if (nrow(out) == 0) {
    return(tibble(experiment_name = character(), experiment_id = integer(),
                  subjid = integer(), kind = character(), path = character()))
  }

  dup <- out %>% count(experiment_id, subjid, kind) %>% filter(n > 1)
  if (nrow(dup) > 0) {
    ex <- out %>% semi_join(dup, by = c("experiment_id", "subjid", "kind")) %>% pull(path)
    stop("Multiple files of the same kind for one subject (remove the stale copy):\n  ",
         paste(head(ex, 10), collapse = "\n  "))
  }
  out
}

# =============================================================================
# read_epoch_metadata()
# Returns: trial (= epoch_index), event_code, condition, intensity
# =============================================================================
read_epoch_metadata <- function(file_path) {
  df <- read_csv(file_path, show_col_types = FALSE,
                 col_types = cols(.default = col_character()))
  names(df) <- clean_names_local(names(df))
  required <- c("epoch_index", "event_code", "condition", "intensity")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Epoch metadata missing columns (", paste(missing, collapse = ", "), "): ", file_path)
  }
  df %>%
    transmute(
      trial      = as.integer(epoch_index),
      event_code = as.character(event_code),
      condition  = str_to_lower(coalesce(condition, "unknown")),
      intensity  = str_to_lower(coalesce(intensity, "unknown"))
    )
}