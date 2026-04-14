# ============================================================
# merge_behavioural.R
# ============================================================
# Merges all behavioural CSVs across experiments into a
# single standardized master dataset.
# ============================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# ============================================================
# USER SETTINGS
# ============================================================
setwd("/home/UWO/darsenea/Documents/GitHub/pain-alpha-dynamics")

behav_dir = "/home/UWO/darsenea/Documents/GitHub/pain-alpha-dynamics/R/analysis/experiment"
output_file = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/behavioural_master.csv"

# ============================================================
# EXPERIMENT LOOKUP
# ============================================================
experiment_lookup <- tibble::tibble(
  experiment_name = c(
    "26ByBiosemi", "29ByANT", "39ByBP",
    "30ByANT", "65ByANT", "95ByBP",
    "142ByBiosemi", "223ByBP", "29ByBP"
  ),
  experiment_id = 1:9
)

# ============================================================
# HELPERS
# ============================================================
clean_names_local <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\.", "_") %>%
    str_replace_all("\\^2", "r2")
}

# Rename
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

extract_experiment_name <- function(file_path, valid_names) {
  fname <- basename(file_path)
  matches <- valid_names[str_detect(fname, fixed(valid_names))]
  if (length(matches) == 1L) return(matches)
  if (length(matches) == 0L) stop("Could not match behavioural file to experiment name: ", fname)
  stop("Multiple experiment name matches found in behavioural file: ", fname)
}

standardize_behaviour <- function(file_path, exp_lookup) {
  message("Reading: ", basename(file_path))

  df <- read_csv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  df <- rename_cols(df, list(
    subjid = c("subjid", "ID", "id", "participant"),
    trial = c("trial", "Trial", "trial_number", "Trial_Num", "trial_num"),
    laser_power = c("laser_power", "laser", "laser_level"),
    pain_rating = c("pain_rating", "pain", "rating", "Painrating")
  ))

  required <- c("subjid", "trial", "laser_power", "pain_rating")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns (", paste0(missing, collapse = ", "), ") in: ", file_path)
  }

  exp_name <- extract_experiment_name(file_path, exp_lookup$experiment_name)
  exp_id <- exp_lookup$experiment_id[exp_lookup$experiment_name == exp_name]

  df %>%
    mutate(
      experiment_name = exp_name, 
      experiment_id = exp_id,
      subjid = as.integer(subjid),
      trial = as.integer(trial),
      laser_power = suppressWarnings(as.numeric(laser_power)),
      pain_rating = suppressWarnings(as.numeric(pain_rating))
      ) %>%
      arrange(experiment_id, subjid, trial) %>%
      group_by(experiment_id, subjid) %>%
      mutate(trial_index = row_number()) %>%
      ungroup() %>%
      mutate(subjid_uid = sprintf("E%02d_S%03d", experiment_id, subjid))
}

# ============================================================
# READ + MERGE
# ============================================================
message("Working directory: ", getwd())
message("behav_dir exists: ", file.exists(behav_dir))

files <- list.files(
  behav_dir,
  pattern = "_behaviour\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
if (length(files) == 0) stop("No behavioural CSVs found in: ", behav_dir)

behaviour_master <- map(files, standardize_behaviour, exp_lookup = experiment_lookup) %>%
  bind_rows()

# ============================================================
# ADD GLOBAL SUBJID
# ============================================================
subject_key <- behaviour_master |>
  distinct(experiment_id, experiment_name, subjid, subjid_uid) |>
  arrange(experiment_id, subjid) |>
  mutate(global_subjid = row_number())

behaviour_master <- behaviour_master |>
  left_join(subject_key,
            by = c("experiment_id", "experiment_name", "subjid", "subjid_uid")) |>
            select(
              experiment_name, experiment_id,
              subjid, subjid_uid, global_subjid,
              trial, trial_index,
              laser_power, pain_rating
            ) |>
            arrange(experiment_id, subjid, trial)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(behaviour_master, output_file)

message("Behavioural master saved: ", output_file)
message(" Rows :", nrow(behaviour_master))
message(" Subjects: ", n_distinct(behaviour_master$subjid_uid))