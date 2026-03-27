# =========================================================
# merge_participants_into_behavioural.R
# ---------------------------------------------------------
# Adds participant demographics (age, sex) from participants.tsv
# and cap size from cap_size.csv into behavioural_master.csv,
# producing behavioural_demo_master.csv
# =========================================================

library(readr)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# =========================================================
# USER SETTINGS
# =========================================================
participants_root <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/"
behav_file <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/behavioural_master.csv"
capsize_file <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/cap_size.csv"
output_file <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/behavioural_demo_master.csv"

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

harmonize_sex <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  case_when(
    x %in% c("m", "male") ~ "M",
    x %in% c("f", "female") ~ "F",
    TRUE ~ NA_character_
  )
}

read_participants <- function(experiment_name, experiment_id, root_dir) {
  file_path <- file.path(root_dir, experiment_name, "participants.tsv")

  if (!file.exists(file_path)) {
    warning("participants.tsv not found for: ", experiment_name)
    return(NULL)
  }

  message("Reading participants: ", file_path)

  df <- read_tsv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  df <- df %>%
    rename(
      participant_id = any_of(c("participant_id", "participant")),
      age            = any_of(c("age", "Age")),
      sex            = any_of(c("sex", "Sex", "gender", "Gender"))
    )

  if (!"participant_id" %in% names(df)) {
    warning("participant_id missing in: ", file_path)
    return(NULL)
  }

  if (!"age" %in% names(df)) df$age <- NA_real_
  if (!"sex" %in% names(df)) df$sex <- NA_character_

  df %>%
    mutate(
      experiment_name = experiment_name,
      experiment_id   = experiment_id,
      subjid          = as.integer(str_remove(participant_id, "^sub-")),
      age             = suppressWarnings(as.numeric(age)),
      sex             = harmonize_sex(sex)
    ) %>%
    select(experiment_name, experiment_id, subjid, age, sex)
}

# =========================================================
# LOAD BEHAVIOURAL MASTER
# =========================================================
if (!file.exists(behav_file)) {
  stop("Behavioural master not found: ", behav_file)
}

behaviour_master <- read_csv(behav_file, show_col_types = FALSE)
names(behaviour_master) <- clean_names_local(names(behaviour_master))

# =========================================================
# READ PARTICIPANTS
# =========================================================
participants_list <- pmap(
  experiment_lookup,
  function(experiment_name, experiment_id) {
    read_participants(experiment_name, experiment_id, participants_root)
  }
) %>%
  compact()

if (length(participants_list) == 0) {
  stop("No participants.tsv files were successfully read.")
}

participants_master <- bind_rows(participants_list) %>%
  distinct(experiment_name, experiment_id, subjid, .keep_all = TRUE)

# =========================================================
# READ CAP SIZE
# =========================================================
if (!file.exists(capsize_file)) {
  stop("cap_size.csv not found: ", capsize_file)
}

cap_df <- read_csv(capsize_file, show_col_types = FALSE)
names(cap_df) <- clean_names_local(names(cap_df))

required_cap_cols <- c("experiment_name", "experiment_id", "cap_size")
missing_cap_cols <- setdiff(required_cap_cols, names(cap_df))
if (length(missing_cap_cols) > 0) {
  stop("cap_size.csv missing required columns: ", paste(missing_cap_cols, collapse = ", "))
}

cap_df <- cap_df %>%
  mutate(
    experiment_name = as.character(experiment_name),
    experiment_id   = as.integer(experiment_id),
    cap_size        = as.factor(cap_size)
  ) %>%
  distinct(experiment_name, experiment_id, .keep_all = TRUE)

# =========================================================
# MERGE EVERYTHING
# =========================================================
behaviour_demo_master <- behaviour_master %>%
  left_join(
    participants_master,
    by = c("experiment_name", "experiment_id", "subjid")
  ) %>%
  left_join(
    cap_df,
    by = c("experiment_name", "experiment_id")
  ) %>%
  select(
    experiment_name,
    experiment_id,
    subjid,
    subjid_uid,
    global_subjid,
    age,
    sex,
    cap_size,
    trial,
    trial_index,
    laser_power,
    pain_rating
  ) %>%
  arrange(experiment_id, subjid, trial)

write_csv(behaviour_demo_master, output_file)

message("Behaviour + demographics + cap size master saved: ", output_file)
message("Rows with missing age: ", sum(is.na(behaviour_demo_master$age)))
message("Rows with missing sex: ", sum(is.na(behaviour_demo_master$sex)))
message("Rows with missing cap_size: ", sum(is.na(behaviour_demo_master$cap_size)))