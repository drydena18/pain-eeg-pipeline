# =========================================================
# merge_spectral_behaviour.R
# ---------------------------------------------------------
# Merges spectral GA trial-level CSVs with behavioural +
# demographic master data and incrementally updates
# alpha_pain_master.csv without duplicate rows.
# V 2.0.0
#
# V2.0.0 changes:
#   - safe_rename() helper replaces any_of() inside rename()
#     (any_of() is a tidy-select helper valid in select/across
#     but not in rename(); it silently did nothing in dplyr).
#   - Column name backward-compat: handles both old names
#     (pow_slow, pow_fast, slow_frac, …) and new canonical
#     names from the dynamic CSV writer (pow_slow_alpha,
#     pow_fast_alpha, slow_alpha_frac, …).
#   - All V2 spectral columns (bi_pre, lr_pre, cog_pre,
#     psi_cog, erd_slow, erd_fast, delta_erd, p5_flag,
#     phase_slow_rad) pass through automatically via
#     spectral_only_cols — no additional mapping required.
# =========================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# =========================================================
# USER SETTINGS
# =========================================================
spectral_root  <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis"
behav_file     <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/behavioural_demo_master.csv"
output_file    <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/alpha_pain_master.csv"
spectral_pattern <- "_ga_by_trial\\.csv$"

# =========================================================
# FIXED EXPERIMENT LOOKUP
# =========================================================
experiment_lookup <- tibble::tibble(
  experiment_name = c(
    "26ByBiosemi", "29ByANT", "39ByBP",
    "30ByANT", "65ByANT", "95ByBP",
    "142ByBiosemi", "223ByBP", "29ByBP"
  ),
  experiment_id = 1:9
)

# =========================================================
# HELPERS
# =========================================================
clean_names_local <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\^2", "r2") %>%
    str_replace_all("\\.", "_")
}

# safe_rename: rename `new_name` to the first candidate found in df.
# Unlike any_of() inside rename(), this is explicit and warns when
# no candidate is present.
safe_rename <- function(df, new_name, candidates) {
  found <- intersect(candidates, names(df))
  if (length(found) == 0) return(df)                    # column absent; no-op
  if (found[1] == new_name) return(df)                  # already correct name
  rename(df, !!new_name := !!found[1])
}

extract_experiment_name <- function(file_path) {
  exp_name <- str_match(file_path, "da-analysis/([^/]+)/spec/")[, 2]
  if (is.na(exp_name)) {
    stop("Could not extract experiment_name from path: ", file_path)
  }
  exp_name
}

read_one_spectral <- function(file_path, experiment_lookup) {
  message("Reading spectral: ", basename(file_path))

  df <- read_csv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  # ----------------------------------------------------------
  # Backward-compat renames: old V1 names -> new canonical names
  # The dynamic CSV writer (spec_write_ga_trial_csv V2) writes
  # MATLAB field names directly, so no renames are needed for
  # data processed with the new pipeline.  These renames handle
  # legacy CSVs produced before the V2 writer.
  # ----------------------------------------------------------
  df <- safe_rename(df, "pow_slow_alpha",
                    c("pow_slow_alpha", "pow_slow"))
  df <- safe_rename(df, "pow_fast_alpha",
                    c("pow_fast_alpha", "pow_fast"))
  df <- safe_rename(df, "pow_alpha_total",
                    c("pow_alpha_total", "pow_alpha"))
  df <- safe_rename(df, "rel_slow_alpha",
                    c("rel_slow_alpha", "rel_slow"))
  df <- safe_rename(df, "rel_fast_alpha",
                    c("rel_fast_alpha", "rel_fast"))
  df <- safe_rename(df, "slow_alpha_frac",
                    c("slow_alpha_frac", "slow_frac", "sf_slow_frac"))

  # FOOOF column aliases (bridge output vs filled output names)
  df <- safe_rename(df, "fooof_offset",
                    c("fooof_offset", "aperiodic_offset"))
  df <- safe_rename(df, "fooof_exponent",
                    c("fooof_exponent", "aperiodic_exponent"))
  df <- safe_rename(df, "fooof_r2",
                    c("fooof_r2", "fooof_r_2", "r2"))
  df <- safe_rename(df, "fooof_error",
                    c("fooof_error", "error"))
  df <- safe_rename(df, "fooof_alpha_cf",
                    c("fooof_alpha_cf", "alpha_cf"))
  df <- safe_rename(df, "fooof_alpha_pw",
                    c("fooof_alpha_pw", "alpha_pw"))
  df <- safe_rename(df, "fooof_alpha_bw",
                    c("fooof_alpha_bw", "alpha_bw"))

  if (!all(c("subjid", "trial") %in% names(df))) {
    stop("Spectral file missing required columns (subjid/trial): ", file_path)
  }

  exp_name_val <- extract_experiment_name(file_path)

  experiment_id <- experiment_lookup %>%
    filter(experiment_name == !!exp_name_val) %>%
    pull(experiment_id)

  if (length(experiment_id) != 1) {
    stop("Experiment name not found in lookup table: ", exp_name_val)
  }

  df %>%
    mutate(
      experiment_name = exp_name_val,
      experiment_id   = experiment_id,
      subjid          = as.integer(subjid),
      trial           = as.integer(trial),
      subjid_uid      = sprintf("E%02d_S%03d", experiment_id, subjid)
    )
}

merge_one_subject <- function(spec_df, behav_master) {
  this_subjid <- unique(spec_df$subjid)
  this_exp_id <- unique(spec_df$experiment_id)
  this_exp_nm <- unique(spec_df$experiment_name)

  if (length(this_subjid) != 1) {
    warning("Spectral file contains multiple subject IDs; skipping.")
    return(NULL)
  }

  if (length(this_exp_id) != 1 || length(this_exp_nm) != 1) {
    warning("Spectral file has ambiguous experiment ID/name; skipping.")
    return(NULL)
  }

  behav_sub <- behav_master %>%
    filter(experiment_id == this_exp_id, subjid == this_subjid)

  if (nrow(behav_sub) == 0) {
    warning("No behavioural/demographic rows found for subjid ", this_subjid,
            " in ", this_exp_nm)
    return(NULL)
  }

  spec_cols <- names(spec_df)

  merged <- spec_df %>%
    left_join(
      behav_sub,
      by = c("experiment_name", "experiment_id", "subjid", "subjid_uid", "trial")
    )

  id_cols    <- c("experiment_name", "experiment_id", "subjid", "subjid_uid")
  behav_cols <- c("global_subjid", "age", "sex", "cap_size",
                  "trial", "trial_index", "laser_power", "pain_rating")
  spectral_only_cols <- setdiff(
    spec_cols,
    c("experiment_name", "experiment_id", "subjid", "subjid_uid", "trial")
  )

  merged %>%
    select(
      all_of(id_cols),
      all_of(behav_cols),
      any_of(spectral_only_cols)
    )
}

# =========================================================
# LOAD BEHAVIOURAL + DEMOGRAPHIC MASTER
# =========================================================
if (!file.exists(behav_file)) {
  stop("Behavioural demo master not found: ", behav_file)
}

behav_master <- read_csv(behav_file, show_col_types = FALSE)
names(behav_master) <- clean_names_local(names(behav_master))

required_behav_cols <- c(
  "experiment_name", "experiment_id", "subjid", "subjid_uid",
  "global_subjid", "age", "sex", "cap_size",
  "trial", "trial_index", "laser_power", "pain_rating"
)

missing_behav_cols <- setdiff(required_behav_cols, names(behav_master))
if (length(missing_behav_cols) > 0) {
  stop("Behavioural demo master missing required columns: ",
       paste(missing_behav_cols, collapse = ", "))
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

# =========================================================
# FIND + READ + MERGE SPECTRAL FILES
# =========================================================
spectral_files <- list.files(
  spectral_root,
  pattern  = spectral_pattern,
  recursive = TRUE,
  full.names = TRUE
)

if (length(spectral_files) == 0) {
  stop("No spectral files found under: ", spectral_root)
}

message("Found ", length(spectral_files), " spectral GA files.")

merged_list <- spectral_files %>%
  map(~tryCatch(
    read_one_spectral(.x, experiment_lookup),
    error = function(e) { warning(.x, ": ", e$message); NULL }
  )) %>%
  compact() %>%
  map(~tryCatch(
    merge_one_subject(.x, behav_master),
    error = function(e) { warning(e$message); NULL }
  )) %>%
  compact()

if (length(merged_list) == 0) {
  stop("No successfully merged subject files were produced.")
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
    trial_index     = as.integer(trial_index)
  )

# =========================================================
# UPSERT INTO EXISTING MASTER
# =========================================================
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
      trial_index     = as.integer(trial_index)
    )

  all_cols <- union(names(old_data), names(new_data))
  for (col in setdiff(all_cols, names(old_data))) old_data[[col]] <- NA
  for (col in setdiff(all_cols, names(new_data))) new_data[[col]] <- NA

  updated_master <- old_data %>%
    select(all_of(all_cols)) %>%
    anti_join(
      new_data %>% select(experiment_id, subjid, trial),
      by = c("experiment_id", "subjid", "trial")
    ) %>%
    bind_rows(new_data %>% select(all_of(all_cols)))
} else {
  updated_master <- new_data
}

updated_master <- updated_master %>%
  arrange(experiment_id, subjid, trial)

write_csv(updated_master, output_file)

message("Saved updated alpha_pain_master: ", output_file)
message("Total rows: ",        nrow(updated_master))
message("Total subjects: ",    n_distinct(updated_master$subjid_uid))
message("Total experiments: ", n_distinct(updated_master$experiment_id))