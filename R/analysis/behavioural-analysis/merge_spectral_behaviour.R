# =============================================
# merge_spectral_behaviour.R
# =============================================
# Merges spectral GA trail-level CSV with behavioural master
# data and incrementally updates alpha_pain_master.csv
# without duplicating existing entries.
# =============================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)

# =============================================
# USER SETTINGS
# =============================================

# Root directory containing spectral files (searched recursively)
spectral_root <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/"

# Behavioural master file from merge_behavioural.R
behav_file <- "pain-eeg-pipeline/R/analysis/behavioural-analysis/behavioural_data.csv"

# Output master file for GAMMs
output_file <- "pain-eeg-pipeline/R/analysis/behavioural-analysis/alpha_pain_master.csv"

# Only merge GA trial-level spectral files
spectral_pattern <- "_spectral_ga_by_trial\\.csv$"

# =============================================
# HELPER: extract dataset/experiment ID from path
# =============================================
extract_experiment_id <- function(file_path) {
    # captures folder name immediately after da-analysis
    exp_id <- str_match(
        file_path,
        "da-analysis/([^/]+)/preproc/"
        )[, 2]

        if (is.na(exp_id)) {
            stop("Could not extract experiment ID from path: ", file_path)
        }

        exp_id
}

# =============================================
# HELPER: standardize one spectral file
# =============================================
read_one_spectral <- function(file_path) {
    message("Reading spectral: ", basename(file_path))

    df <- read_csv(file_path, show_col_types = FALSE)

    # Clean names
    names(df) <- names(df) %>%
        str_trim() %>%
        str_replace_all("\\s+", "_") %>%
        str_replace_all("\\^2", "r2") %>%
        str_replace_all("\\.", "_")
    
    # Standardize core identifiers
    df <- df %>%
        rename(
            subjid = any_of(c("subjid", "ID", "id", "participant")),
            trial  = any_of(c("trial", "trial_num", "trial_number")),
            slow_frac = any_of(c("slow_frac", "sf_slow_frac")),
            fooof_r2 = any_of(c("fooof_r2", "fooof_r^2", "r2"))
        )
    
    # Extract experiment ID from folder structure
    exp_id <- extract_experiment_id(file_path)

    # Basic checks
    if (!all(c("subjid", "trial") %in% names(df))) {
        stop("Spectral file missing required columns (subjid, trial): ", file_path)
    }

    df %>%
        mutate(
            experiment = exp_id
            )
}

# =============================================
# HELPER: merge one spectral df with behavioural master
# =============================================
merge_one_subject <- function(spec_df, behav_master) {
    this_subjid <- unique(spec_df$subjid)
    this_exp <- unique(spec_df$experiment_id)

    if (!length(this_subjid) != 1) {
        warning("Spectral file contains multiple subject IDs. Skipping.")
        return(NULL)
    }

    if (length(this_exp) != 1) {
        warning("Spectral file has ambiguous experiment ID. Skipping.")
        return(NULL)
    }

    behav_sub <- behav_master %>%
        filter(subjid == this_subjid, experiment_id == this_exp)
    
    if (nrow(behav_sub) == 0) {
        warning("No behavioural rows found for ", this_subjid, " in ", this_exp, ". Skipping.")
        return(NULL)
    }

    merged <- spec_df %>%
        left_join(
            behav_sub,
            by = c("experiment_id", "subjid", "trial")
        )

    n_missing_pain = sum(is.na(merged$pain_rating))
    n_missing_laser = sum(is.na(merged$laser_power))

    if (n_missing_pain > 0 || n_missing_laser > 0) {
        warning(
            paste0("Merge warning for ", this_subjid, " (", this_exp, "): ",
            "missing pain_rating = ", n_missing_pain,
            ", missing laser_power = ", n_missing_laser
            )
        )   
    }

    merged
}

# =============================================
# Load Behavioural Master
# =============================================
