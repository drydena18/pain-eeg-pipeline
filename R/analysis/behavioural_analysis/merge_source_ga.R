# =============================================================================
# merge_source_ga.R
# -----------------------------------------------------------------------------
# Collates per-subject grand-average (GA) source CSVs into a single group-level
# GA table.  Covers both pre+post spectral / LEP metrics and FOOOF metrics.
#
# Inputs (auto-discovered):
#   <source_root>/**/sub-XXX/csv/sub-XXX_source_ga.csv
#   <source_root>/**/sub-XXX/csv/sub-XXX_source_ga_fooof.csv  (optional)
#
# Outputs:
#   source_ga_master.csv        one row per (experiment, subject, ROI)
#   source_ga_fooof_master.csv  one row per (experiment, subject, ROI) — FOOOF only
#
# These GA masters are used for:
#   - Subject-level analyses of TVI_alpha, ITC, FOOOF aperiodic parameters
#   - Seeding run_gamm_source.R with GA predictors where appropriate
#   - run_classical_tests.R Tests 3 and 4
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# =============================================================================
# USER SETTINGS
# =============================================================================
source_root    <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/"
behav_file     <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/behavioural_demo_master.csv"
out_ga         <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/source_ga_master.csv"
out_ga_fooof   <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/source_ga_fooof_master.csv"

ga_pattern     <- "_source_ga\\.csv$"
fooof_pattern  <- "_source_ga_fooof\\.csv$"

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
# HELPERS
# =============================================================================
clean_names_local <- function(x) {
  x %>% str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\^2",  "r2") %>%
    str_replace_all("\\.",   "_")
}

extract_experiment_name <- function(file_path) {
  # Path: .../da-analysis/<exp_name>/source/sub-XXX/csv/...
  m <- str_match(file_path, "da-analysis/([^/]+)/source/")
  if (is.na(m[, 2])) stop("Cannot extract experiment_name: ", file_path)
  m[, 2]
}

read_one_ga <- function(file_path, experiment_lookup) {
  message("  Reading: ", basename(file_path))
  df <- tryCatch(read_csv(file_path, show_col_types = FALSE),
                 error = function(e) { warning(file_path, ": ", conditionMessage(e)); NULL })
  if (is.null(df)) return(NULL)
  names(df) <- clean_names_local(names(df))
  
  # Normalise subject column
  if ("subject" %in% names(df) && !"subjid" %in% names(df)) df <- rename(df, subjid = subject)
  if (!"subjid" %in% names(df)) { warning("No subjid column: ", file_path); return(NULL) }
  if (!"roi"    %in% names(df)) { warning("No roi column: ",    file_path); return(NULL) }
  
  exp_name <- tryCatch(extract_experiment_name(file_path),
                       error = function(e) { warning(conditionMessage(e)); NA_character_ })
  if (is.na(exp_name)) return(NULL)
  
  exp_id <- experiment_lookup %>% filter(experiment_name == exp_name) %>% pull(experiment_id)
  if (length(exp_id) != 1L) { warning("Experiment not in lookup: ", exp_name); return(NULL) }
  
  df %>% mutate(
    experiment_name = exp_name,
    experiment_id   = exp_id,
    subjid          = as.integer(subjid),
    subjid_uid      = sprintf("E%02d_S%03d", exp_id, subjid),
    roi             = as.character(roi)
  )
}

# =============================================================================
# LOAD BEHAVIOURAL DEMOGRAPHICS  (for age / sex join)
# =============================================================================
behav_demo <- NULL
if (file.exists(behav_file)) {
  behav_demo <- read_csv(behav_file, show_col_types = FALSE)
  names(behav_demo) <- clean_names_local(names(behav_demo))
  behav_demo <- behav_demo %>%
    distinct(experiment_id, subjid, .keep_all = TRUE) %>%
    select(experiment_id, subjid, subjid_uid, global_subjid, age, sex, cap_size) %>%
    mutate(experiment_id = as.integer(experiment_id),
           subjid        = as.integer(subjid))
} else {
  message("[WARN] behav_demo_master not found — demographic columns will be absent.")
}

# =============================================================================
# DISCOVER + READ GA FILES
# =============================================================================
ga_files    <- list.files(source_root, pattern = ga_pattern,    recursive = TRUE, full.names = TRUE)
fooof_files <- list.files(source_root, pattern = fooof_pattern, recursive = TRUE, full.names = TRUE)

message("Found ", length(ga_files),    " GA files.")
message("Found ", length(fooof_files), " FOOOF GA files.")

# =============================================================================
# GA MASTER
# =============================================================================
if (length(ga_files) > 0L) {
  ga_list <- map(ga_files, ~read_one_ga(.x, experiment_lookup)) %>% compact()
  
  if (length(ga_list) > 0L) {
    ga_master <- bind_rows(ga_list) %>%
      mutate(across(where(is.numeric), ~suppressWarnings(as.numeric(.x)))) %>%
      arrange(experiment_id, subjid, roi)
    
    # Join demographics
    if (!is.null(behav_demo)) {
      ga_master <- ga_master %>%
        left_join(behav_demo, by = c("experiment_id", "subjid", "subjid_uid"))
    }
    
    dir.create(dirname(out_ga), recursive = TRUE, showWarnings = FALSE)
    write_csv(ga_master, out_ga)
    message("Saved source_ga_master: ", out_ga)
    message("  Rows      : ", nrow(ga_master))
    message("  Subjects  : ", n_distinct(ga_master$subjid_uid))
    message("  ROIs      : ", paste(sort(unique(ga_master$roi)), collapse = ", "))
    message("  Columns   : ", paste(names(ga_master), collapse = ", "))
  } else {
    message("[WARN] No GA files read successfully.")
  }
} else {
  message("[WARN] No GA files found.")
}

# =============================================================================
# FOOOF GA MASTER
# =============================================================================
if (length(fooof_files) > 0L) {
  fooof_list <- map(fooof_files, ~read_one_ga(.x, experiment_lookup)) %>% compact()
  
  if (length(fooof_list) > 0L) {
    fooof_master <- bind_rows(fooof_list) %>%
      mutate(across(where(is.numeric), ~suppressWarnings(as.numeric(.x)))) %>%
      arrange(experiment_id, subjid, roi)
    
    if (!is.null(behav_demo)) {
      fooof_master <- fooof_master %>%
        left_join(behav_demo, by = c("experiment_id", "subjid", "subjid_uid"))
    }
    
    write_csv(fooof_master, out_ga_fooof)
    message("Saved source_ga_fooof_master: ", out_ga_fooof)
    message("  Rows    : ", nrow(fooof_master))
    message("  Columns : ", paste(names(fooof_master), collapse = ", "))
  } else {
    message("[WARN] No FOOOF GA files read successfully.")
  }
} else {
  message("[WARN] No FOOOF GA files found — run source pipeline with fooof.enabled = true first.")
}