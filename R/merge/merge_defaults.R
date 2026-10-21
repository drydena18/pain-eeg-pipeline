# =============================================================================
# merge_defaults.R
# =============================================================================
# Configuration for the merge stage (merge_core.R). Run merge_core.R from
# R/merge/ (see R_README.md, Run Sequence).
#
# V2 (Oct 2026) — matches pipeline feature set v3.1 (spectral_core.m 3.1.0,
# source_core.py 4.1.0):
#   - experiment_lookup used `<-` inside list(), which made
#     merge_cfg$experiment_lookup NULL; now `=`.
#   - Per-experiment output folders are listed explicitly (no recursive
#     search of the whole da-analysis tree, no file-name -> experiment guess).
#   - Adds epoch metadata, subject-level summaries, trial-alignment QC and a
#     rebuild switch for derived tables (subject_split is never dropped).
# =============================================================================

merge_cfg <- list(

    # -- Roots ----
    # <data_root>/<experiment_name>/{preproc, spec, source, resource}
    data_root = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis",

    # Behavioural CSVs (one per experiment: <experiment_name>_behaviour.csv)
    behav_dir = "/home/UWO/darsenea/Documents/GitHub/pain-alpha-dynamics/R/data/behaviour",

    cap_size_file = "/home/UWO/darsenea/Documents/GitHub/pain-alpha-dynamics/R/data/cap_size.csv",

    # -- Output database ----
    data_db_path = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/db/data.duckdb",

    # Pipeline feature set these tables are built from. Written to ingest_log
    # and checked by the GAMM stage.
    feature_set = "v3.1",

    # -- Rebuild ----
    # TRUE: drop and recreate every DERIVED table (subjects, trials, epochs,
    # alignment_qc, *_metrics, ingest_log) at the start of each run, so the
    # database always mirrors the files on disk. subject_split is NEVER
    # dropped: it is the only state that must persist across runs.
    rebuild_derived_tables = TRUE,

    # -- File layout (per experiment) ----
    behav_pattern = "_behaviour\\.csv$",
    participants_filename = "participants.tsv",
    files = list(
        epoch_meta       = list(subdir = "preproc", pattern = "^sub-\\d+_epoch_metadata\\.csv$"),
        spectral_trial   = list(subdir = "spec",    pattern = "^sub-\\d+_spectral_chan_by_trial\\.csv$"),
        spectral_subject = list(subdir = "spec",    pattern = "^sub-\\d+_spectral_chan_summary\\.csv$"),
        source_trial     = list(subdir = "source",  pattern = "^sub-\\d+_source_trial\\.csv$"),
        source_subject   = list(subdir = "source",  pattern = "^sub-\\d+_source_ga\\.csv$")
    ),

    # -- Trial alignment ----
    # Behaviour is joined to EEG by trial = epoch index. A subject is
    # EXCLUDED from every EEG table unless:
    #   - EEG epoch count == behavioural trial count (and the spectral/source
    #     files agree with the epoch count), and
    #   - for experiments whose epochs carry low/high labels, every "high"
    #     epoch has the subject's higher laser_power and every "low" epoch
    #     the lower one (exact; laser order is randomized, so a one-trial
    #     shift breaks ~half the trials).
    # Status per subject is written to alignment_qc either way.
    exclude_misaligned = TRUE,

    # -- Rating range ----
    # Subjects with ANY rating outside this range (data-entry errors, e.g. a
    # 34 on a 0-10 scale) get no EEG rows; listed in alignment_qc with
    # status "rating_out_of_range". NA ratings are not affected.
    rating_range = c(0, 10),

    # -- Experiment lookup ----
    experiment_lookup = tibble::tibble(
        experiment_name = c(
            "26ByBiosemi", "29ByANT", "39ByBP",
            "30ByANT", "65ByANT", "95ByBP",
            "142ByBiosemi", "223ByBP", "29ByBP"
        ),
        experiment_id = 1:9
    ),

    # -- Train / holdout split ----
    # Incremental, stratified by cap_size, never reassigned once set (see
    # helpers/split_helpers.R). Do not change the seed once any subject has
    # been assigned.
    split_seed = 42L,
    split_train_frac = 0.75
)