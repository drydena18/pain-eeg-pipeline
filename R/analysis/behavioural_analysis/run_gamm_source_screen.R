# =============================================================================
# run_gamm_source_screen.R
# =============================================================================
# Combo-capable replacement for run_gamm_core.R
#
# Screens every candidate alpha metrics (Balance Index, Log-Ratio, Centre of
# Gravity, Temporal Variability Index, delta-ERD, sf_balance, FOOOF aperiodic
# params, ...) against every ROI in the six-region pain neuromatrix
# (S1/postcentral, M1/precentral, ACC, insula, SII/supramarginal,
# dlPFC/rostralmiddlefrontal - bilaterl), to find which metric x region
# combinations show real signal. This is Aim 1b's exploratory screen.
#
# Output:
#   DuckDB table gamm_fitted in <out_dir>/gamm_results.duckdb
#       (one row per roi x model_name, written incrementally as each fit
#        completes - see gamm_combo_core.R for the upsert/crash-safety design)
#   <out_dir>/<roi>/<model_name>.rds        every fitted model object
#   <out_dir>/<roi>/model_input_<roi>.csv   z-scored model input (QC trail)
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)

source("gamm_combo_core.R")

# =============================================================================
# USER SETTINGS
# =============================================================================
data_file <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/R/source_pain_master.csv"
out_dir <- "/cofs/seminowicz/eegPainDatasets/CNED/da-analysis/R/gamm_outputs_source_screen"
db_path <- file.path(out_dir, "gamm_results.duckdb")

min_fooof_r2 <- 0.80
min_trials_roi <- 20L

COMBO_SIZES <- 1L # <- flip to c(1L, 2L) or c(1L, 2L, 3L)
TENSOR_PAIRS <- list(c("pow_slow_z", "pow_fast_z")) # fit as te() not s()+s()

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================
clean_names_local <- function(x) {
  x %>% str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("\\^2",  "r2") %>%
    str_replace_all("\\.",   "_")
}
zscore      <- function(x) as.numeric(scale(x))
safe_factor <- function(x) as.factor(as.character(x))

# =============================================================================
# LOAD DATA
# =============================================================================
if (!file.exists(data_file)) {
    stop("source_pain_master.csv not found: ", data_file, "\nRun merge_source_spectral.R first.")
}

df_raw <- read_csv(data_file, show_col_types = FALSE)
names(df_raw) <- clean_names_local(names(df_raw))

required_cols <- c("experiment_name", "experiment_id")