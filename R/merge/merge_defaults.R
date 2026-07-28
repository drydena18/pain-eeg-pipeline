# =============================================================================
# merge_defaults.R
# =============================================================================
# Configuration for the merge stage (merge_core.R + assign_split.R).
# Replaces the per-script harcoded paths/lookups duplicated across the old
# 5 merge_*.R scripts.
# =============================================================================

merge_cfg <- list(

    # -- Source Roots ----
    # Root containing per-experiment subfolders, each with resource/,
    # spec/sub-XXX/csv/, source/sub-XXX/csv
    data_root = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis",

    # Behavioural CSVs live in a separate tree from the per-subject pipeline
    # output
    behav_dir = "/home/UWO/darsenea/Documents/GitHub/pain-alpha-dynamics/R/analysis/experiment",

    cap_size_file = "/home/UWO/darsenea/Documents/GitHub/pain-alpha-dynamics/R/analysis/behavioural_analysis/cap_size.csv",

    # -- Output Database ----
    data_db_path = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/db/data.duckdb",

    # -- File Patterns (per-subject pipeline outputs) ----
    behav_pattern = "_behaviour\\.csv$",
    spectral_pattern = "_spectral_chan_by_trial\\.csv$",
    source_pattern = "_source_trial\\.csv$",
    participants_filename = "participants.tsv",

    # -- Experiment Lookup ----
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
    ),

    # -- Train / holdout split ----
    # Incremental, stratified by cap_size, never reassigned once set (see
    # helpers/split_helpers.R). Fixed seed for reproducibility; do not change
    # once any subject has been assigned, or existing assignments won't be
    # reproducible from a fresh run
    split_siid = 42L,
    split_train_frac = 0.75
)