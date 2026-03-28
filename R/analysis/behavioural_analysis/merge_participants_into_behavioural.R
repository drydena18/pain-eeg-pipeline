# =========================================================
# merge_participants_into_behavioural.R
# ---------------------------------------------------------
# Adds participant demographics (age, sex) from participants.tsv
# and cap size from cap_size.csv into behavioural_master.csv,
# producing behavioural_demo_master.csv.
#
# FIX LOG:
#   - Removed duplicate library(readr).
#   - rename(col = any_of(...)) is invalid dplyr syntax;
#     replaced with rename_cols() helper (same as
#     merge_behavioural.R).
#   - filter(experiment_name == !!experiment_name) had a
#     name collision; inner variable renamed to exp_name_val.
# =========================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# =========================================================
# USER SETTINGS
# =========================================================
participants_root <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/"
behav_file   <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/behavioural_master.csv"
capsize_file <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/cap_size.csv"
output_file  <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/behavioural_demo_master.csv"

# =========================================================
# FIXED EXPERIMENT LOOKUP
# =========================================================
experiment_lookup <- tibble::tribble(
  ~experiment_name,  ~experiment_id,
  "26ByBiosemi",     1L,
  "29ByANT",         2L,
  "39ByBP",          3L,
  "30ByANT",         4L,
  "65ByANT",         5L,
  "95ByBP",          6L,
  "142ByBiosemi",    7L,
  "223ByBP",         8L,
  "29ByBP",          9L
)

# =========================================================
# HELPERS
# =========================================================

clean_names_local <- function(x) {
  x |>
    str_trim() |>
    str_replace_all("\\s+", "_") |>
    str_replace_all("\\.", "_") |>
    str_replace_all("\\^2", "r2")
}

rename_cols <- function(df, mapping) {
  cur <- names(df)
  for (new_nm in names(mapping)) {
    candidates <- mapping[[new_nm]]
    hit <- intersect(candidates, cur)
    if (length(hit) > 0 && !(new_nm %in% cur)) {
      df  <- rename(df, !!new_nm := !!hit[1])
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

read_participants <- function(exp_name_val, exp_id_val, root_dir) {
  file_path <- file.path(root_dir, exp_name_val, "participants.tsv")

  if (!file.exists(file_path)) {
    warning("participants.tsv not found for: ", exp_name_val)
    return(NULL)
  }

  message("Reading participants: ", file_path)

  df <- read_tsv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  # FIX: rename(col = any_of(...)) is invalid; use helper
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

  df |>
    mutate(
      experiment_name = exp_name_val,
      experiment_id   = exp_id_val,
      subjid          = as.integer(str_remove(participant_id, "^sub-")),
      age             = suppressWarnings(as.numeric(age)),
      sex             = harmonize_sex(sex)
    ) |>
    select(experiment_name, experiment_id, subjid, age, sex)
}

# =========================================================
# LOAD BEHAVIOURAL MASTER
# =========================================================
if (!file.exists(behav_file)) stop("Behavioural master not found: ", behav_file)

behaviour_master <- read_csv(behav_file, show_col_types = FALSE)
names(behaviour_master) <- clean_names_local(names(behaviour_master))

# =========================================================
# READ PARTICIPANTS
# =========================================================
participants_master <- pmap(
  experiment_lookup,
  function(experiment_name, experiment_id) {
    # FIX: renamed args to avoid collision with tibble column names
    read_participants(experiment_name, experiment_id, participants_root)
  }
) |>
  compact() |>
  bind_rows() |>
  distinct(experiment_name, experiment_id, subjid, .keep_all = TRUE)

# =========================================================
# READ CAP SIZE
# =========================================================
if (!file.exists(capsize_file)) stop("cap_size.csv not found: ", capsize_file)

cap_df <- read_csv(capsize_file, show_col_types = FALSE)
names(cap_df) <- clean_names_local(names(cap_df))

required_cap_cols <- c("experiment_name", "experiment_id", "cap_size")
missing_cap <- setdiff(required_cap_cols, names(cap_df))
if (length(missing_cap) > 0) {
  stop("cap_size.csv missing required columns: ", paste(missing_cap, collapse = ", "))
}

cap_df <- cap_df |>
  mutate(
    experiment_name = as.character(experiment_name),
    experiment_id   = as.integer(experiment_id),
    cap_size        = as.factor(cap_size)
  ) |>
  distinct(experiment_name, experiment_id, .keep_all = TRUE)

# =========================================================
# MERGE EVERYTHING
# =========================================================
behaviour_demo_master <- behaviour_master |>
  left_join(participants_master,
            by = c("experiment_name", "experiment_id", "subjid")) |>
  left_join(cap_df,
            by = c("experiment_name", "experiment_id")) |>
  select(
    experiment_name, experiment_id,
    subjid, subjid_uid, global_subjid,
    age, sex, cap_size,
    trial, trial_index,
    laser_power, pain_rating
  ) |>
  arrange(experiment_id, subjid, trial)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(behaviour_demo_master, output_file)

message("Saved: ", output_file)
message("  Rows missing age      : ", sum(is.na(behaviour_demo_master$age)))
message("  Rows missing sex      : ", sum(is.na(behaviour_demo_master$sex)))
message("  Rows missing cap_size : ", sum(is.na(behaviour_demo_master$cap_size)))