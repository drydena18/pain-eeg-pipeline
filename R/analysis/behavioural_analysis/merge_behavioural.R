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
behav_dir = "pain-eeg-pipeline/R/analysis/behavioural_analysis/experiment"
output_file = "pain-eeg-pipeline/R/analysis/behavioural_analysis/behavioural_master.csv"

# ============================================================
# EXPERIMENT LOOKUP
# ============================================================
experiment_lookup = tibble::tibble(
  ~experiment_name, ~experiment_id,
  "26ByBiosemi", 1L,
  "29ByANT", 2L,
  "39ByBP", 3L,
  "30ByANT", 4L,
  "65ByANT", 5L,
  "95ByBP", 6L,
  "142ByBiosemi", 7L,
  "223ByBP", 8L,
  "29ByBP", 9L
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

extract_experiment_name <- function()