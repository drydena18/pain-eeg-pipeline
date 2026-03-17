### --------------------
# Packages / libraries
# ----------------------
library(readr)
library(dplyr)
library(purrr)
library(stringer)

# =============================================
# USER SETTINGS
# =============================================

# Directory containing all behvaioural CSVs
behav_dir <- "pain-eeg-pipeline/R/analysis/experiment/"

# Output file
output_file <- "pain-eeg-pipeline/R/analysis/behavioural-analysis/behavioural_data.csv"

# =============================================
# HELPER FUNCTION
# =============================================
standardize_behaviour <- function(file_path) {
  message('Reading: ', basename(file_path))

  df <- read_csv(file_path, show_col_types = FALSE)

  # Standardize column names
  df <- df %>%
    rename_with(~ str_trim(.x)) %>%
    rename(
      subjid = any_oc(c("subjid", "ID", "id", "participant")),
      trial  = any_of(c("trial", "trial_num", "trial_number", "Trialnum")),
      laser_power = any_of(c("laser_power", "laser", "stim_intensity", "intensity")),
      pain_rating = any_of(c("pain_rating", "pain", "rating", "Painrating"))
    )

    # Extract experiment ID from filename
    exp_id <- str_extract(basename(file_path), "exp\\d+")

    if(is.na(exp_id)) {
      exp_id <- basename(file_path) %>%
        str_remove("\\.csv$")
    }

    # Add metadata
    df <- df %>%
      mutate(
        experiment = exp_id
      )

    # Create trial index (within subject)
    df <- df %>%
      arrange(subjid, trial) %>%
      group_by(subjid) %>%
      mutate(trial_index = row_number()) %>%
      ungroup()

    # Sanity checks
    if (!all(c("subjid", "trial", "laser_power", "pain_rating") %in% names(df))) {
      stop("Missing required columns in file: ", file_path)
    }

    df
}

# =============================================
# READ ALL FILES
# =============================================
files <- list.files(
  behav_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop('No behavioural files found in directory.')
}

# =============================================
# PROCESS + MERGE
# =============================================
behaviour_list <- map(files, standardize_behaviour)

behaviour_master <- bind_rows(behaviour_list)

# =============================================
# FINAL CLEANING
# =============================================
behaviour_master <- behaviour_master %>%
  mutate(
    subjid = as.factor(subjid),
    experiment_id = as.factor(experiment_id)
  ) %>%
  select(
    experiment_id,
    subjid,
    trial,
    trial_index,
    laser_power,
    pain_rating
  )

# =============================================
# SAVE OUTPUT
# =============================================
write_csv(behaviour_master, output_file)

message("Behavioural data merged and saved to: ", output_file)

# =============================================
# Quick QC
# =============================================
message("Total rows: ", nrow(behaviour_master))
message("Total subjects: ", n_distinct(behaviour_master$subjid))

qc <- behaviour_master %>%
  group_by(experiment_id, subjid) %>%
  summarise(
    n_trials = n(),
    missing_pain = sum(is.na(pain_rating)),
    missing_laser = sum(is.na(laser_power)),
    .groups = "drop"
  )

print(qc)