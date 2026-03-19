# =========================================================
# merge_behavioural.R
# ---------------------------------------------------------
# Merges all behavioural CSVs across experiments into a
# single standardized behavioural dataset.
# =========================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# =========================================================
# USER SETTINGS
# =========================================================
behav_dir <- "pain-eeg-pipeline/R/analysis/experiment/"
output_file <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/behavioural_master.csv"

# =========================================================
# FIXED EXPERIMENT LOOKUP
# =========================================================
experiment_lookup <- tibble::tribble(
  ~experiment_name,   ~experiment_id,
  "26ByBiosemi",      1,
  "29ByANT",          2,
  "39ByBP",           3,
  "30ByANT",          4,
  "65ByANT",          5,
  "95ByBP",           6,
  "142ByBiosemi",     7,
  "223ByBP",          8,
  "29ByBP",           9
)

# =========================================================
# HELPERS
# =========================================================
clean_names_local <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\.", "_") %>%
    str_replace_all("\\^2", "r2")
}

extract_experiment_name_from_file <- function(file_path, valid_names) {
  fname <- basename(file_path)
  matches <- valid_names[str_detect(fname, fixed(valid_names))]

  if (length(matches) == 1) return(matches)
  if (length(matches) == 0) stop("Could not match behavioural file to experiment name: ", fname)
  stop("Multiple experiment name matches found in behavioural file: ", fname)
}

standardize_behaviour <- function(file_path, experiment_lookup) {
  message("Reading: ", basename(file_path))

  df <- read_csv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  df <- df %>%
    rename(
      subjid      = any_of(c("subjid", "ID", "id", "participant")),
      trial       = any_of(c("trial", "trial_num", "trial_number")),
      laser_power = any_of(c("laser_power", "laser", "stim_intensity", "intensity")),
      pain_rating = any_of(c("pain_rating", "pain", "rating"))
    )

  if (!all(c("subjid", "trial", "laser_power", "pain_rating") %in% names(df))) {
    stop("Missing required columns in behavioural file: ", file_path)
  }

  experiment_name <- extract_experiment_name_from_file(
    file_path,
    experiment_lookup$experiment_name
  )

  experiment_id <- experiment_lookup %>%
    filter(experiment_name == !!experiment_name) %>%
    pull(experiment_id)

  df %>%
    mutate(
      experiment_name = experiment_name,
      experiment_id   = experiment_id,
      subjid          = as.integer(subjid),
      trial           = as.integer(trial)
    ) %>%
    arrange(experiment_id, subjid, trial) %>%
    group_by(experiment_id, subjid) %>%
    mutate(trial_index = row_number()) %>%
    ungroup() %>%
    mutate(
      subjid_uid = sprintf("E%02d_S%03d", experiment_id, subjid)
    )
}

# =========================================================
# READ + MERGE
# =========================================================
files <- list.files(
  behav_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No behavioural CSVs found in directory: ", behav_dir)
}

behaviour_list <- map(files, ~standardize_behaviour(.x, experiment_lookup))
behaviour_master <- bind_rows(behaviour_list)

# =========================================================
# ADD GLOBAL SUBJECT ID
# =========================================================
subject_key <- behaviour_master %>%
  distinct(experiment_id, experiment_name, subjid, subjid_uid) %>%
  arrange(experiment_id, subjid) %>%
  mutate(global_subjid = row_number())

behaviour_master <- behaviour_master %>%
  left_join(
    subject_key,
    by = c("experiment_id", "experiment_name", "subjid", "subjid_uid")
  ) %>%
  select(
    experiment_name,
    experiment_id,
    subjid,
    subjid_uid,
    global_subjid,
    trial,
    trial_index,
    laser_power,
    pain_rating
  ) %>%
  arrange(experiment_id, subjid, trial)

write_csv(behaviour_master, output_file)

message("Behavioural master saved: ", output_file)
message("Total rows: ", nrow(behaviour_master))
message("Total subjects: ", n_distinct(behaviour_master$subjid_uid))