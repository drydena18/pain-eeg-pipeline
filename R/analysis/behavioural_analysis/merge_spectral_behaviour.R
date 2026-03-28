# =========================================================
# merge_spectral_behaviour.R
# ---------------------------------------------------------
# Merges spectral GA trial-level CSVs with behavioural +
# demographic master data and incrementally upserts into
# alpha_pain_master.csv without duplicate rows.
#
# V 2.0.0 changes
# ---------------
# * Type enforcement for new pre-stim interaction metric
#   columns: bi_pre, lr_pre, cog_pre, psi_cog, delta_erd,
#   phase_slow_rad (numeric) and p5_flag (integer).
#   These columns pass through read_one_spectral() via the
#   spectral_only_cols passthrough — no rename alias needed.
# * Mixed-vintage guard: check_spectral_columns() warns when
#   a pre-V2 CSV is missing the new columns so mixed batches
#   degrade gracefully (old rows get NA for new columns).
# * Removed dead slow_frac rename alias (the GA CSV never
#   wrote a column named "slow_frac"; it was always written
#   as "slow_alpha_frac" and passed through automatically).
# * Removed "r2" from the fooof_r2 alias list (too broad;
#   would clobber any column literally named r2).
# * Unified pipe style to base R |> (R >= 4.1).
# =========================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

# =========================================================
# USER SETTINGS
# =========================================================
spectral_root   <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/"
behav_file      <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/behavioural_demo_master.csv"
output_file     <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/alpha_pain_master.csv"
spectral_pattern <- "_ga_by_trial\\.csv$"

# =========================================================
# FIXED EXPERIMENT LOOKUP
# =========================================================
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

# =========================================================
# HELPERS
# =========================================================
clean_names_local <- function(x) {
  x |>
    str_trim() |>
    str_replace_all("\\s+", "_") |>
    str_replace_all("\\^2",  "r2") |>
    str_replace_all("\\.",   "_")
}

extract_experiment_name <- function(file_path) {
  exp_name <- str_match(file_path, "da-analysis/([^/]+)/preproc/")[, 2]
  if (is.na(exp_name)) {
    stop("Could not extract experiment_name from path: ", file_path)
  }
  exp_name
}

# Columns introduced in V2 of the spectral pipeline.
# Used only to emit informative warnings on mixed-vintage batches;
# their absence is not an error — old rows simply get NA.
v2_spectral_cols <- c(
  "bi_pre", "lr_pre", "cog_pre", "psi_cog",
  "delta_erd", "phase_slow_rad", "p5_flag"
)

check_spectral_columns <- function(df, file_path) {
  missing_v2 <- setdiff(v2_spectral_cols, names(df))
  if (length(missing_v2) > 0) {
    message(
      "[WARN] Pre-V2 spectral file (missing new columns): ", basename(file_path),
      "\n       Missing: ", paste(missing_v2, collapse = ", "),
      "\n       These rows will receive NA for those columns in the master CSV."
    )
  }
}

# Enforce canonical types on all spectral metric columns.
# Uses mutate(across(any_of(...))) so it is safe when columns are absent.
normalize_spectral_types <- function(df) {
  numeric_cols <- c(
    # Full-epoch features (V1)
    "paf_cog_hz",
    "pow_slow_alpha", "pow_fast_alpha", "pow_alpha_total",
    "rel_slow_alpha",  "rel_fast_alpha",
    "sf_ratio", "sf_logratio", "sf_balance", "slow_alpha_frac",
    # FOOOF
    "fooof_offset", "fooof_exponent", "fooof_r2", "fooof_error",
    "fooof_alpha_cf", "fooof_alpha_pw", "fooof_alpha_bw",
    "fooof_alpha_cf_filled", "fooof_alpha_pw_filled", "fooof_alpha_bw_filled",
    "fooof_alpha_found",
    # Pre-stim interaction metrics (V2)
    "bi_pre", "lr_pre", "cog_pre", "psi_cog",
    "erd_slow", "erd_fast", "delta_erd",
    "phase_slow_rad"
  )
  integer_cols <- c("p5_flag")

  df |>
    mutate(across(any_of(numeric_cols),  ~ suppressWarnings(as.numeric(.x)))) |>
    mutate(across(any_of(integer_cols),  ~ suppressWarnings(as.integer(.x))))
}

# ---------------------------------------------------------------
read_one_spectral <- function(file_path, experiment_lookup) {
  message("Reading spectral: ", basename(file_path))

  df <- read_csv(file_path, show_col_types = FALSE)
  names(df) <- clean_names_local(names(df))

  # Rename only columns whose canonical name may vary across legacy files.
  # New V2 columns (bi_pre, lr_pre, cog_pre, …) are already written with
  # their canonical names by spec_write_ga_trial_csv.m and pass through
  # the spectral_only_cols mechanism without renaming.
  df <- df |>
    rename(
      subjid         = any_of(c("subjid",         "ID",         "id",         "participant")),
      trial          = any_of(c("trial",           "trial_num",  "trial_number")),
      fooof_offset   = any_of(c("fooof_offset",   "aperiodic_offset")),
      fooof_exponent = any_of(c("fooof_exponent", "aperiodic_exponent")),
      # fooof_r2 alias: "r2" deliberately omitted — too broad, would clobber
      # any column literally named r2. FOOOF writes fooof_r2 directly.
      fooof_r2       = any_of(c("fooof_r2",       "fooof_r_2")),
      fooof_error    = any_of(c("fooof_error",    "error")),
      fooof_alpha_cf = any_of(c("fooof_alpha_cf", "alpha_cf")),
      fooof_alpha_pw = any_of(c("fooof_alpha_pw", "alpha_pw")),
      fooof_alpha_bw = any_of(c("fooof_alpha_bw", "alpha_bw"))
    )

  if (!all(c("subjid", "trial") %in% names(df))) {
    stop("Spectral file missing required columns (subjid/trial): ", file_path)
  }

  # Warn if this is a pre-V2 file missing new interaction-metric columns
  check_spectral_columns(df, file_path)

  # Enforce types before any join
  df <- normalize_spectral_types(df)

  experiment_name <- extract_experiment_name(file_path)

  experiment_id <- experiment_lookup |>
    filter(experiment_name == !!experiment_name) |>
    pull(experiment_id)

  if (length(experiment_id) != 1) {
    stop("Experiment name not found in lookup table: ", experiment_name)
  }

  df |>
    mutate(
      experiment_name = experiment_name,
      experiment_id   = experiment_id,
      subjid          = as.integer(subjid),
      trial           = as.integer(trial),
      subjid_uid      = sprintf("E%02d_S%03d", experiment_id, subjid)
    )
}

# ---------------------------------------------------------------
merge_one_subject <- function(spec_df, behav_master) {
  this_subjid <- unique(spec_df$subjid)
  this_exp_id <- unique(spec_df$experiment_id)
  this_exp_nm <- unique(spec_df$experiment_name)

  if (length(this_subjid) != 1) {
    warning("Spectral file contains multiple subject IDs; skipping.")
    return(NULL)
  }

  if (length(this_exp_id) != 1 || length(this_exp_nm) != 1) {
    warning("Spectral file has ambiguous experiment ID / name; skipping.")
    return(NULL)
  }

  behav_sub <- behav_master |>
    filter(
      experiment_id == this_exp_id,
      subjid        == this_subjid
    )

  if (nrow(behav_sub) == 0) {
    warning("No behavioural/demographic rows for subjid ", this_subjid,
            " in ", this_exp_nm)
    return(NULL)
  }

  spec_cols <- names(spec_df)

  merged <- spec_df |>
    left_join(
      behav_sub,
      by = c("experiment_name", "experiment_id", "subjid", "subjid_uid", "trial")
    )

  id_cols           <- c("experiment_name", "experiment_id", "subjid", "subjid_uid")
  behav_cols        <- c("global_subjid", "age", "sex", "cap_size",
                         "trial", "trial_index", "laser_power", "pain_rating")
  # All spectral columns from this file (including V2 columns if present)
  spectral_only_cols <- setdiff(
    spec_cols,
    c("experiment_name", "experiment_id", "subjid", "subjid_uid", "trial")
  )

  merged |>
    select(
      all_of(id_cols),
      all_of(behav_cols),
      all_of(spectral_only_cols)
    )
}

# =========================================================
# LOAD BEHAVIOURAL + DEMOGRAPHIC MASTER
# =========================================================
if (!file.exists(behav_file)) {
  stop("Behavioural demo master not found: ", behav_file)
}

behav_master        <- read_csv(behav_file, show_col_types = FALSE)
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

behav_master <- behav_master |>
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
# FIND SPECTRAL FILES
# =========================================================
spectral_files <- list.files(
  spectral_root,
  pattern   = spectral_pattern,
  recursive = TRUE,
  full.names = TRUE
)

if (length(spectral_files) == 0) {
  stop("No spectral GA CSV files found under: ", spectral_root)
}

message("Found ", length(spectral_files), " spectral GA CSV files.")

# =========================================================
# READ + MERGE
# =========================================================
merged_list <- spectral_files |>
  map(\(fp) read_one_spectral(fp, experiment_lookup)) |>
  map(\(df) merge_one_subject(df, behav_master)) |>
  compact()

if (length(merged_list) == 0) {
  stop("No successfully merged subject files were produced.")
}

# ---------------------------------------------------------------
# Bind rows across subjects.
# bind_rows() fills absent columns with NA, so mixed-vintage
# files (some with V2 columns, some without) combine cleanly.
# ---------------------------------------------------------------
new_data <- bind_rows(merged_list) |>
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
  ) |>
  # Re-apply type enforcement after bind_rows() in case any column
  # was promoted to character during the row-bind of mixed files
  normalize_spectral_types()

message("New data rows: ", nrow(new_data),
        " | Subjects: ", n_distinct(new_data$subjid_uid))

# Report V2 column presence in the incoming batch
v2_present <- intersect(v2_spectral_cols, names(new_data))
v2_absent  <- setdiff(v2_spectral_cols, names(new_data))
if (length(v2_present) > 0) {
  message("V2 interaction-metric columns present in batch: ",
          paste(v2_present, collapse = ", "))
}
if (length(v2_absent) > 0) {
  message("[WARN] V2 columns absent from entire batch (all files are pre-V2): ",
          paste(v2_absent, collapse = ", "))
}

# =========================================================
# UPSERT INTO EXISTING MASTER
# =========================================================
if (file.exists(output_file)) {
  message("Existing master found; upserting...")

  old_data        <- read_csv(output_file, show_col_types = FALSE)
  names(old_data) <- clean_names_local(names(old_data))

  old_data <- old_data |>
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
    ) |>
    normalize_spectral_types()

  # Union of all columns; pad whichever side is missing a column with NA
  all_cols <- union(names(old_data), names(new_data))
  for (col in setdiff(all_cols, names(old_data))) old_data[[col]] <- NA
  for (col in setdiff(all_cols, names(new_data))) new_data[[col]] <- NA

  old_data <- old_data |> select(all_of(all_cols))
  new_data <- new_data |> select(all_of(all_cols))

  # Drop old rows for any (experiment_id, subjid, trial) appearing in the
  # incoming batch, then bind the new versions.  This is a whole-row replace;
  # no column-level merging is attempted.
  n_old <- nrow(old_data)
  updated_master <- old_data |>
    anti_join(
      new_data |> select(experiment_id, subjid, trial),
      by = c("experiment_id", "subjid", "trial")
    ) |>
    bind_rows(new_data)

  message("Upsert: replaced ",
          n_old - (nrow(updated_master) - nrow(new_data)),
          " old rows with ", nrow(new_data), " new rows.")

} else {
  message("No existing master found; creating from scratch.")
  updated_master <- new_data
}

updated_master <- updated_master |>
  arrange(experiment_id, subjid, trial)

write_csv(updated_master, output_file)

message("Saved updated alpha_pain_master.csv -> ", output_file)
message("Total rows     : ", nrow(updated_master))
message("Total subjects : ", n_distinct(updated_master$subjid_uid))
message("Total experiments: ", n_distinct(updated_master$experiment_id))